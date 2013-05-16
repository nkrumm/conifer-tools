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
	if chromosomes == None:
		chromosomes = range(1,24)
	
	for s in samples:
		for chrom in chromosomes:
			data = p.getConiferData(s,chrom)
			yield {"data":data,"sample":s,"chrom":chrom}

def mapinit_func():
	global p
	p = cp.ConiferPipeline() # empty pipeline for methods only

def map_func(x):
	data = p.preprocess(x["data"])
	attempts = 0
	while attempts < (args.n_retry+1):
		try:
			print "[STATUS] now running %s:chr%s" % (x["sample"], x["chrom"])
			calls = p.segment(data,x["sample"],x["chrom"])
			print "[STATUS] finished %s:chr%s" % (x["sample"], x["chrom"])
			return calls
		except:
			print "[WARNING] %s:chr%d could not segment... Retrying" % (x["sample"], x["chrom"])
			attempts +=1
	print "[ERROR] %s:chr%d failed to segment after %d tries" % (x["sample"], x["chrom"], (args.n_retry+1))
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