from mpi_wrapper import *
import numpy.random as nr
import conifertools as cp
import time
import argparse
import re

def get_data(conifer_file, samples = None, chromosomes=None):
    
    p = cp.ConiferPipeline(conifer_file)
    if samples == None:
        samples = p.samples
    familyIDs = sorted(list(set(map(lambda x: x[0:5], samples))))
    if chromosomes == None:
        chromosomes = range(1,24)
    
    #for s in samples:
    #   data = p.getConiferData(s,chromosomes)
    #   yield {"data":data,"sample":s,"chrom":chromosomes}
    for f in familyIDs:
        data = []
        sampleIDs = []
        for rel in ["mo","fa","p1","s1","s2","s3"]:
            sampleID = "%s.%s" % (f, rel)
            if sampleID in p.samples:
                sampleIDs.append(sampleID)
                d = p.preprocess(p.getConiferData(sampleID, range(1,24)))
                data.append(d)
        yield {"data":data,"samples":sampleIDs,"chrom":chromosomes,"familyID": f}

def mapinit_func():
    global p
    p = cp.ConiferPipeline() # empty pipeline for methods only

def map_func(x):
    attempts = 0
    while attempts < (args.n_retry+1):
        #try:
        print "[STATUS] now running %s" % x["familyID"]
        calls = p.segment(x["data"], x["samples"])
        print "[STATUS] finished %s" % x["familyID"]
        return calls
 #       except:
 #           print "[WARNING] %s could not segment... Retrying" % x["familyID"]
 #           attempts +=1
    print "[ERROR] %s failed to segment after %d tries" % (x["familyID"], (args.n_retry+1))
    return 0


def reduceinit_func():
    global all_calls
    all_calls = cp.CallTable()


def reduce_func(x):
    if x != 0:
        print x.calls
        all_calls.appendCalls(x)

def reduceexit_func():
    try:
        all_calls.calls.to_csv(out_file,sep="\t")
    except RuntimeError:
        print "[RUN ERROR]"
        print all_calls


if __name__ == "__main__":  
    global out_file
    # parse arguments
    parser = argparse.ArgumentParser(prog="call.mpi.py", description="SGE submit script for ConiferPipeline class. Author: Nik Krumm, 2013")
    parser.add_argument('--out',action='store', required=True, metavar='/path/to/calls.txt',  nargs="?",help="output location of calls.txt file")
    parser.add_argument('--conifer_file',action='store', required=True, metavar='/path/to/analysis.hdf5',nargs="?",help="location of conifer data file")
    parser.add_argument('--chromosomes',type=int, nargs="+", required=False, metavar='2 or 1 2 3 4, etc.',help="Chromosomes to process")
    parser.add_argument('--samples',type=str, nargs="+", required=False, metavar='sampleA sampleB etc.',help="Samples to process")
    parser.add_argument('--n_retry',type=int, default=0, required=False, metavar='3',help="Max number of retries of CGHCall method before failing with no calls")
    parser.add_argument('--verbose',action="store_true", default=False, required=False,help="Inidicate verbose output")
    args = parser.parse_args()
    
    out_file = args.out
    
    start_time = time.time()
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank   = comm.Get_rank()
    
    data_iter = get_data(args.conifer_file, chromosomes=args.chromosomes, samples=args.samples)
    
    if rank==0:
        masterloop(comm = comm, rank = 0, nprocs = nprocs,\
                    data_function = data_iter,\
                    reduce_function = reduce_func,\
                    reduceinit_function = reduceinit_func,\
                    reduceexit_function = reduceexit_func)
        
        print "Total Time: ",  time.time()-start_time
        
    else:
        slaveloop(comm = comm, rank= rank,\
                    init_function=mapinit_func,\
                    map_function = map_func)