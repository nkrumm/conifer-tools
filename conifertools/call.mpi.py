from mpi_wrapper import *
import numpy.random as nr
import numpy as np
import conifertools as cp
import time
import argparse
import re

def get_data(conifer_file, samples = None, families = None, chromosomes=None):
    """
    Retrieves data from the conifer analysis.hdf5 file and yields a dictionary of data.
    With no parameters, will yield all samples and all chromosomes by sample.
    If samples and/or chromosomes is given, the data yielded will be restricted to these.
    If families is also given, the samples will be yielded in groups based on the values 
        in the families variable. For this to work, the list of samples and the list of 
        families should be *the same length*. eg:
        samples = [100.p1, 100.mo, 100.fa, 201.p1, 201.mo]
        families = [100, 100, 100, 201, 201]
    """
    p = cp.ConiferPipeline(conifer_file)
    if samples == None:
        samples = p.samples

    if chromosomes == None:
        chromosomes = range(1,24)

    if families is not None:
        if len(families) != len(samples):
            print "[FATAL ERROR] list of samples and families not the same length!"
            sys.exit(1)

        families = np.array(families)
        samples = np.array(samples)
        
        for f in np.unique(families):
            data = []
            sampleIDs = list(samples[families == f])
            for sampleID in sampleIDs: 
                d = p.preprocess(p.getConiferData(sampleID, chromosomes))
                data.append(d)
            yield {"data":data,"samples":sampleIDs,"chrom":chromosomes,"familyID": f}
    else:
        # just iterate through the samples
        for s in samples:
            data = p.getConiferData(s,chromosomes)
            yield {"data":data,"sample":s,"chrom":chromosomes}

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
        if calls != 0:
            return calls
        attempts += 1

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
    parser.add_argument('--families',type=str, nargs="+", required=False, metavar='familyA familyB etc.',help="Family IDs corresponding to samples in --samples list")
    parser.add_argument('--n_retry',type=int, default=0, required=False, metavar='3',help="Max number of retries of CGHCall method before failing with no calls")
    parser.add_argument('--verbose',action="store_true", default=False, required=False,help="Inidicate verbose output")
    args = parser.parse_args()
    
    out_file = args.out
    
    start_time = time.time()
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank   = comm.Get_rank()
    
    data_iter = get_data(args.conifer_file, chromosomes=args.chromosomes, families=args.families, samples=args.samples)
    
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