import drmaa
import os
import subprocess
import tempfile
import stat
import time

class LocalJobTemplate(drmaa.JobTemplate):
    queue = "all.q"
    
    def __enter__(self):
        instance = super(LocalJobTemplate, self).__enter__()
        instance.workingDirectory = os.getcwd()
        instance.joinFiles = False
        instance.nativeSpecification = "-q %s -S /bin/bash" % self.queue
        instance.outputPath=":%s" % os.environ['HOME']
        return instance
    
    
def build_module_command(module_list):
    cmd = ". /etc/profile.d/modules.sh\n"
    cmd +="if test ! -z $MODULESHOME; then \n"
    for m in module_list:
        cmd += "module load " + m + " \n"
    cmd += "fi\n"
    return cmd

def build_qsubscript(remote_command, modules=None, mpirun=False):
    
    header = "#!/bin/bash\n"
    
    if modules != None:
        module_command = build_module_command(modules)
    else:
        module_command = ""
    if mpirun:
        mpi_command = "mpirun -x PATH -x LD_LIBRARY_PATH --prefix $MPIBASE -mca plm ^rshd -mca btl ^openib "
    else:
        mpi_command = ""
    
    
    script = " ".join([header, module_command, mpi_command, remote_command, "\n"])
    
    script_f = tempfile.NamedTemporaryFile(delete=False)
    #script_f = open("tmp000.sh",'w')
    script_fname = script_f.name
    script_f.write(script)
    script_f.close()
    return script_fname
    
def submit_job(remote_command, name="mpirun", modules = None, native_specification=None, mpirun=False, wait=False):
    """
    Submit a single job to the cluster.

    The optional ``native_specification`` argument specifies parameters for
    the SGE job such as the amount of memory and CPUs to request. The format
    is a string as it would appear on the command line in a qsub
    command. For example, the following string requests 4 GB of RAM and
    10-20 CPUs:
    
        -h_vmem=4G -pe orte 10-20
    
    Takes an optional argument, ``wait``, which tells the method to wait for
    the cluster jobs to finish before returning.
    """
    job_id = None
    with drmaa.Session() as session:
        with LocalJobTemplate() as job_template:
            job_template.workingDirectory = os.getcwd()
            job_template.joinFiles = True
            job_template.nativeSpecification = "-S /bin/bash -b n"
            
            script_fname = build_qsubscript(remote_command,modules,mpirun)
            job_template.remoteCommand = script_fname
            if native_specification is not None:
                job_template.nativeSpecification = " ".join([job_template.nativeSpecification, native_specification])
            
            job_template.jobName = name
            
            
            try:
                job_id = session.runJob(job_template)
            except drmaa.errors.DeniedByDrmException, e:
                print e
                return -1
            
            if wait:
                return_value = session.wait(job_id, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        
    return job_id   