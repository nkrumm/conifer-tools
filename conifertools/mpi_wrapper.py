from mpi4py import MPI
import sys

TAG_WORKER_FINISHED = 20
TAG_RAW_DATA = 10
TAG_SETUP = 0
TAG_KILL = 5

def masterloop(comm, rank, nprocs, data_function, reduce_function, reduceinit_function, reduceexit_function):
	completedRuns = 0
	num_workers = nprocs - 1
	workerQueue = []
	
	reduceinit_function()
	
	for i in range(1, nprocs):
		comm.send("START",dest = i, tag = TAG_SETUP)
		print "master sending START to ", str(i)
	while len(workerQueue) < num_workers:
		worker_rank = comm.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG)
		print "master received OK from ", worker_rank
		workerQueue.append(worker_rank)
		sys.stdout.flush()
	
	stillWorking = True
	currentJobs = {}
	
	while stillWorking or (len(currentJobs) > 0):
		msg = comm.Iprobe(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG)
		if msg:
			status = MPI.Status()
			data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
			if status.tag == TAG_WORKER_FINISHED:
				# a worker has finished, do something
				reduce_function(data)
				
				del currentJobs[status.source]
				# this worker is ready for another batch!
				workerQueue.append(status.source)	
				sys.stdout.flush()
		
		if len(workerQueue) > 0:
			# get next worker ID
			worker_rank = workerQueue.pop()
			# get data from file
			try:
				
				out_data = data_function.next()
				
				currentJobs[worker_rank] = 1
				# send the data
				comm.send(out_data, dest = worker_rank, tag = TAG_RAW_DATA)
			except StopIteration:
				stillWorking = False
			
			
	
	### KILL ALL PROCESSES WHEN DONE
	for i in range(1, nprocs):
		comm.send("KILL",dest = i, tag = TAG_KILL)
		print "master sending KILL to ", str(i)
	
	reduceexit_function()


def slaveloop(comm, rank, init_function, map_function):
	data = comm.recv(source=0, tag=TAG_SETUP)
	#print rank, data
	
	init_function()
	
	comm.send(rank) #send startup message
	stillWorking = True
	while stillWorking:
		msg = comm.Iprobe(source=0,tag=MPI.ANY_TAG)
		if msg:
			status = MPI.Status()
			data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
			
			
			if status.tag == TAG_RAW_DATA:
				result = map_function(data)
				comm.send(result, dest=0, tag=TAG_WORKER_FINISHED)
			elif status.tag == TAG_KILL:
				stillWorking = False

