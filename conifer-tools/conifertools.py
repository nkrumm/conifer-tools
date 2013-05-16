import numpy as np
import pandas
try:
    import rpy2
    import rpy2.rlike.container as rlc
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    cghcall = importr("CGHcall")
    cghbase = importr("CGHbase")
except:
    print "Import Error: Could not import rpy2 and/or CGHcall. Installed?"
import sys
import os
import time
import types
import hcluster
try:
    import rpkm
except ImportError:
    print "Import Error: Could not import rpkm class-- is pytables installed and working?"


class CallTable(object):
    """A CallTable contains the call output from a ConiferPipeline analysis. Each row in the CallTable
    has the following fields:
    sampleID            The sample with the call
    chromosome          chromosome or contig name
    start               genomic start position
    stop                genomic stop position
    state               Call is a deletion (-1) or duplication (+1)
    start_exon          exon index of start
    stop_exon           exon index of stop
    num_probes          number of exons
    size_bp             size of call in basepairs
    probability         Probability of call (output from DNACopy algorithm)
    median_svdzrpkm     Median SVD-ZRPKM signal across probes in call
    stdev_svdzrpkm      Standard deviation of SVD-ZRPKM signal across probes in call
    """
    default_columns=["sampleID","chromosome","start","stop","state","start_exon","stop_exon","num_probes","size_bp","probability","median_svdzrpkm","stdev_svdzrpkm"]

    def __init__(self,calls=None):
        converters = {"sampleID": str}
        if isinstance(calls,pandas.DataFrame):
            self.calls = calls
        elif isinstance(calls,str):
            if os.path.exists(calls):
                self.calls = pandas.read_csv(calls,sep="\t",index_col=0, converters = converters)
        elif isinstance(calls, list):
            self.calls = pandas.read_csv(calls[0], sep="\t", index_col=0, converters = converters)
            print calls[0]
            print self.calls
            for p in calls[1:]:
                if os.path.exists(p):
                    print p
                    self.calls = self.calls.append(
                        pandas.read_csv(p, sep="\t", index_col=0, converters = converters),
                        ignore_index=True)
                
        else:
            self.calls = pandas.DataFrame(columns=self.default_columns)
    
    def __repr__(self):
        outstr= "CallTable Object\n"
        outstr+="Summary:\n"
        outstr+="  %d Call(s)" % len(self.calls)
        outstr+="  %d Samples(s)" % len(set(self.calls.sampleID))
        if "cnvrID" in self.calls:
            outstr+="  %d CNVRs" % len(set(self.calls.cnvrID))
        return outstr
    
    def __iter__(self):
        for ix, call in self.calls.iterrows():
            yield call
    
    #def __getitem__(self, name):
        # Todo, may want to return CallsTable() instance here for consistancy
    #   return self.calls[name]
    
    #def __getattr__(self, name):
        # Todo, may want to return CallsTable() instance here for consistancy
    #   return self.calls[name]
    
    def __delitem__(self,name):
        del self.calls[name]
    
    def __delattr__(self,name):
        del self.calls[name]
    
    def head(self, n=10):
        return self.calls.head(n)
    
    def _error(self, t):
        print t
        return 0
    
    def __add__(self,other):
        return CallTable(self.calls.append(other.calls,ignore_index=True))
    
    def appendCalls(self,a):
        """append another CallTable object to this one"""
        self.calls = self.calls.append(a.calls,ignore_index=True)
        
    def _add(self, a):
        try:
            assert isinstance(a, dict)
        except AssertionError, e:
            raise( AssertionError( "Call to add must be of type dict or equivalent!. %s"%e ) )
        
        self.calls = self.calls.append(pandas.Series(a),ignore_index=True)
    
    def addCall(self, a):
        if isinstance(a,list):
            for i in a:
                self._add(i)
        else:
            self._add(a)
    
    
    def filter(self, f):
        if isinstance(f, CallFilterTemplate):
            return f(self)
        elif type(f) == types.FunctionType:
            #TODO, use pyparsing or yacc to allow for simple boolean statements to be evaluated
            return CallTable(self.calls[map(lambda x: f(x[1]), self.calls.iterrows())])
        else:
            print "Filter is not of supported type (lambda or CallFilterTemplate)"
            return 0
    
    def _genColumnName(self,name, tbl):
        if name not in tbl:
            return name
        else:
            i = 1
            while "%s_%d" % (name, i) in tbl:
                i += 1
            return "%s_%d" % (name, i)
    
    def annotate(self, f, name=None):
        if isinstance(f, CallFilterTemplate):
            return f(self,filter_on=False)
        elif type(f) == types.FunctionType:
            if name != None:
                self.calls[name] = map(lambda x: f(x[1]), self.calls.iterrows())
            else:
                name,i = "annotation", 1
                while "%s_%d" % (name, i) in self.calls:
                    i += 1
                self.calls["%s_%d" % (name, i)] = map(lambda x: f(x[1]), self.calls.iterrows())
            return CallTable(self.calls)
        else:
            print "Annotater is not of supported type (lambda or CallFilterTemplate)"
            return 0
    
    
    def clusterCalls(c,gamma=0.9,cophenetic_cutoff=0.85):
        """find CNVRs in calls"""
        cnvrs = hcluster.hcluster_cnvr(c.calls,gamma=gamma,cophenetic_cutoff=cophenetic_cutoff)
        #print cnvrs
        return CallTable(cnvrs)
    
    def clusterCallsByCohort(c,cohort_field="cohort",gamma=0.9,cophenetic_cutoff=0.85):
        """find CNVRs in calls"""
        cnvrs = hcluster.hcluster_cnvr_by_cohort(c.calls,cohort_field,gamma=gamma,cophenetic_cutoff=cophenetic_cutoff)
        #print cnvrs
        return CallTable(cnvrs)

    def mergeCalls(self,merge_window):
        """merge calls within merge_window probe distance"""
        print "Not Implemented"
        pass
    
    def save(self,f,**kwargs):
        if kwargs.get("sep") == None:
            kwargs["sep"] = "\t"
        
        self.calls.to_csv(f,**kwargs)

class CallFilterTemplate():
    def __init__(self, p, file, name, type=None, func = None, coerce_chrXY_to_int=True):
        """Create a new CallFilter from a bed file and exon list"""
        try:
            bedfile = pandas.read_csv(file,sep="\t",index_col=False,names=["chr","start","stop"])
        except:
            print "Exception: Could not read bed file or load previously save filter from", file

        bedfile = bedfile.sort(["chr","start","stop"])
        bedfile_chrs = set(bedfile["chr"].values)
        self._filter = pandas.DataFrame()
        
        for probe_tbl in p.r.h5file.root.probes:
            probe_array = pandas.DataFrame(probe_tbl.read()[["start","stop"]])
            probe_array_chr = int(probe_tbl.title.replace("chr",""))
            probe_array["chr"] = probe_array_chr

            # semi-hack to convert chrX/Y to chr23/chr24
            if coerce_chrXY_to_int and (probe_array_chr not in bedfile_chrs):
                if probe_array_chr == 23:
                    chr_string = "chrX"
                elif probe_array_chr == 24:
                    chr_string = "chrY"
                else:
                    chr_string = probe_tbl.title
            else:
                chr_string = probe_tbl.title

            bed_array = bedfile[bedfile["chr"] == chr_string]
            # this line returns if the intervals of the probe_array are contained or overlapping the bed_array.
            # It works by using the trick that the searchsorted function will return differing indices if a start or stop
            # from bed_array falls between one of the probes, in that case the searchsorted result will not be the same.    
            probe_array["feature"] = np.searchsorted(np.sort(bed_array["start"].values), probe_array["stop"].values) != np.searchsorted(np.sort(bed_array["stop"].values), probe_array["start"].values)
            self._filter = self._filter.append(probe_array,ignore_index=True)
        
        if type == "overlap":
            self.type="overlap"
            self.func = func
            self.name = name
        elif type=="count":
            self.type="count"
            self.func = func
            self.name = name
        elif type=="contains":
            self.type="contains"
            self.func = func
            self.name = name
        else:
            print "filter type %s not yet implemented" % type
            return None
        
    
    def _count(self, row):
        return np.sum((self._filter["feature"] == True) & (self._filter["chr"].values == row["chromosome"]) & (self._filter["start"].values<=row["stop"]) & (self._filter["stop"].values >= row["start"]))
    
    def _frac(self,count, row):
        return float(count)/row["num_probes"]
    
    def _contains(self, row):
        return self._count(row) > 0 
    
    def _genColumnName(self,name, tbl):
        if name not in tbl:
            return name
        else:
            i = 1
            while "%s_%d" % (name, i) in tbl:
                i += 1
            return "%s_%d" % (name, i)
    
    def __call__(self,in_table,filter_on=True):
        if self.type=="overlap":
            if filter_on and (self.func is not None):
                ffunc = lambda x: self.func(self._frac(self._count(x[1]),x[1]))
                return CallTable(in_table.calls.ix[map(ffunc, in_table.calls.iterrows())])
            else:
                in_table.calls[self._genColumnName(self.name,in_table.calls)] = map(lambda x: self._frac(self._count(x[1]),x[1]), in_table.calls.iterrows())
                return in_table
        elif self.type=="count":
            if filter_on and (self.func is not None):
                return CallTable(in_table.calls.ix[map(lambda x: self.func(self._count(x[1])), in_table.calls.iterrows())])
            else:
                in_table.calls[self._genColumnName(self.name,in_table.calls)] = map(lambda x: self._count(x[1]), in_table.calls.iterrows())
                return in_table
        if self.type=="contains":
            if filter_on and (self.func is not None):
                ffunc = lambda x: self.func(self._contains(x[1]))
                return CallTable(in_table.calls.ix[map(ffunc, in_table.calls.iterrows())])
            else:
                in_table.calls[self._genColumnName(self.name,in_table.calls)] = map(lambda x: self._contains(x[1]), in_table.calls.iterrows())
                return in_table
        else:
            print "filter type %s not yet implemented" % self.type
            return in_table
    
    
    def save(self, filename):
        #f = open(filename,'wb')
        #print self.__dict__
        #tmp_dict = {}
        #pickle.dump(self.__dict__,f,2)
        #f.close()
        print "not implemented"
        pass

class ConiferPipeline:
    def __init__(self,conifer_file=None,samples = None):
        """A ConiferPipeline is designed to segment and call CNVs from a CoNIFER analysis file"""    
        if conifer_file is not None:
            self.h5file = conifer_file
            # TODO test if valid file
            try:
                self.r = rpkm.rpkm_reader(self.h5file)
            except:
                print "Error in reading or opening HDF5 analysis file"
                return None
            # get contig names
            self.contigs = self.r.getContigList()
            
            # get sample names
            if samples == None:
                self.samples = self.r.getSampleList()
            else:
                self.samples = samples
        else:
            pass
        
        rpy2.rinterface.set_writeconsole(self.Rlog)
        
        self.stderr_log = ""
        self.stdout_log = ""
        self.log_level = 0
    
    def __del__(self):
        try:
            self.r.h5file.close()
        except:
            pass
    
    @property
    def samples(self):
        #Samples to be processed in this pipeline. The default for the ConiferPipeline
        return self.__samples__
    
    @samples.setter
    def samples(self, value):
        self.__samples__ = value
    
    @property
    def contigs(self):
        #Samples to be processed in this pipeline. The default for the ConiferPipeline
        return self.__contigs__
    
    @contigs.setter
    def contigs(self, value):
        self.__contigs__ = value
        
    def Rlog(self, x):
        self.log("stdout",x)
    
    def log(self, logtype, text):
        if logtype =="stderr":
            self.stderr_log += text
        elif logtype =="stdout":
            self.stdout_log += text
        #print text
    
    def is_sequence(self, arg):
        return (not hasattr(arg, "strip") and hasattr(arg, "__getitem__") or hasattr(arg, "__iter__"))
    
    def getProbesByChrom(self, chrom):
        if isinstance(chrom, str):
            chrom = int(chrom.lstrip("chr"))
        return pandas.DataFrame(self.r.h5file.root.probes._f_getChild("probes_chr%d" % chrom).read())

    def getConiferData(self, sample, chrom):
        #import gzip
        #
        #fn = "/net/eichler/vol8/home/nkrumm/EXOMES/Quad_Analysis/calling/UW_v2/chr%(chrom)d/%(sample)s.chr%(chrom)d.svdzrpkm.txt.gz" % {"sample":sample,"chrom":chrom}
        #d = pandas.read_csv(fn,compression="gzip",names=["CHROMOSOME","START","STOP","svdzrpkm"],sep="\t")
        return self.r.getChromosomeBySample(sample,chrom)
    
    def preprocess(self, data, min_val=-3,max_val=3):
        data.rpkm = np.clip(data.rpkm,min_val,max_val)
        return data
    
    def findBreakPoints(self, x):
        #bkpts = np.where(np.logical_xor(x[1:],x[:-1],))[0]
        bkpts = np.where(np.logical_xor(np.hstack([x[0],x]),np.hstack([x,x[-1]])))[0]
        if x[0] == True:
            bkpts = np.hstack([0,bkpts])
        if x[-1] == True:
            bkpts = np.hstack([bkpts,len(x)])
        if len(bkpts) <= 2:
            if bkpts[1]-bkpts[0] == 1:
                return [tuple([bkpts[0],bkpts[0]])]
            return [tuple(bkpts)]
        else:
            return [(b[0],b[0]) if (b[1]-b[0])==1 else tuple(b) for b in bkpts.reshape(len(bkpts)/2,2)]
    
    def runCGHCall(self, data,sampleID):
        if isinstance(data,rpkm.rpkm_data):
            od = rlc.OrdDict([("NAME", robjects.vectors.IntVector(data.exons["probeID"])), \
                  ("CHROMOSOME",  robjects.vectors.IntVector(np.repeat(data.contig,len(data.exons)))), \
                  ("START_POS",   robjects.vectors.IntVector(data.exons["start"])), \
                  ("STOP_POS",    robjects.vectors.IntVector(data.exons["stop"])),
                  (str(sampleID), robjects.vectors.FloatVector(data.rpkm))])
        else:
            print "unknown data type input!"
            return 0
        
        try:
            t  = cghbase.make_cghRaw(robjects.DataFrame(od))
            t2 = cghcall.preprocess(t)
            t3 = cghcall.segmentData(t2,**{'alpha':0.02, 'undo.splits': "sdundo", 'undo.SD':2})
            t4 = cghcall.postsegnormalize(t3)
            t5 = cghcall.CGHcall(t4)
            try:
                t6 = cghcall.ExpandCGHcall(t5,t4)
            except:
                self.log("stderr", "%(sampleID)s: Failed to ExpandCGHcall() on chromosome %(chrom)d\n" % {"sampleID":sampleID,"chrom":data["CHROMOSOME"][0]})
                print "%(sampleID)s: Failed to ExpandCGHcall() on chromosome %(chrom)d\n" % {"sampleID":sampleID,"chrom":data["CHROMOSOME"][0]}
                raise
        except:
            #print "[ERROR]: Failed on R code for sample %s" % sampleID
            raise
        
        out=np.vstack([np.array(cghbase.chromosomes(t6)),
                       np.array(cghbase.bpstart(t6)),
                       np.array(cghbase.bpend(t6)),
                       np.array(cghbase.copynumber(t6)).flatten(),
                       np.array(cghbase.calls(t6)).flatten(),
                       np.array(cghbase.probdloss(t6), dtype=np.float).flatten(),
                       np.array(cghbase.probloss(t6), dtype=np.float).flatten(),
                       np.array(cghbase.probnorm(t6), dtype=np.float).flatten(),
                       np.array(cghbase.probgain(t6), dtype=np.float).flatten(),
                       np.array(cghbase.probamp(t6), dtype=np.float).flatten()]).transpose()
        
        out = pandas.DataFrame(out,columns=["chromosome","start","stop","svdzrpkm","call","p_dloss","p_loss","p_norm","p_gain","p_amp"])
        return out
    
    def segment(self, data, sample,chromosome):
        # convert "gain" and "amplification" to "gain" only
        callint2colname = {-2: "p_dloss", -1: "p_loss", 1: "p_gain", 2:"p_amp"}
        # fire up R, run CGHCall
        try:
            out = self.runCGHCall(data, sample)
        except:
            print "[ERROR] in runCGHCall() Method, %d, %s" % (len(data), sample)
            pass
        
        # post process 
        out["call"] = np.sign(out["call"])
        out["p_gain"] = np.maximum(out["p_gain"], out["p_amp"])
        out["p_loss"] = np.maximum(out["p_loss"], out["p_dloss"])
        call_states_present = set(out["call"])-set([0])
        
        c = CallTable()
        for i in call_states_present:
            calls = self.findBreakPoints(np.array(1 * (out["call"] == i)))
            for call in calls:
                prob_vals = out.ix[call[0]:call[1]-1][callint2colname[i]]
                svdzrpkm_vals = out.ix[call[0]:call[1]-1]["svdzrpkm"]
                start_bp = out.ix[call[0]]["start"]
                stop_bp = out.ix[call[1]-1]["stop"]
                start_exon= np.searchsorted(out["start"], start_bp)
                stop_exon = np.searchsorted(out["stop"], stop_bp)
                #np.median(prob_vals), np.median(svdzrpkm_vals), np.std(svdzrpkm_vals))                
                t = {"sampleID": sample,\
                     "chromosome":int(chromosome),\
                     "start":int(start_bp),\
                     "stop":int(stop_bp),\
                     "state":int(i),\
                     "start_exon":int(start_exon),\
                     "stop_exon":int(stop_exon),\
                     "num_probes":int(stop_exon-start_exon+1),\
                     "size_bp":int(stop_bp - start_bp),\
                     "probability": np.median(prob_vals),
                     "median_svdzrpkm": np.median(svdzrpkm_vals),
                     "stdev_svdzrpkm":np.std(svdzrpkm_vals)}

                
                c.addCall(t)
        
        return c
    
    def makeCalls(self,chromosomes=None,samples=None,verbose=False):
        if chromosomes == None:
            chromosomes = self.contigs
        if not self.is_sequence(chromosomes):
            chromosomes = [chromosomes]
        
        if samples == None:
            samples = self.samples
        if not self.is_sequence(samples):
            samples = [samples]
        
        calls = CallTable()
        for s in samples:
            sample_calls = CallTable()
            if self.log_level >=1:
                print s, 
            for chrom in chromosomes:
                if self.log_level >=1:
                    print chrom,
                data = self.preprocess(self.getConiferData(s,chrom))
                try:
                    chrom_calls = self.segment(data,s,chrom)
                except:
                    #skip this chromosome
                    print "Failed on chromosome, ", chrom
                    continue
                sample_calls.appendCalls(chrom_calls)
            if len(sample_calls.calls) > 0:
                calls.appendCalls(sample_calls)
                if self.log_level >=1:
                    print "\n", # newline needed
        return calls
    
    
    def makeCallsMPI(self, outfile=None, chromosomes=None,samples=None,n_cpus=10, verbose=False,n_retry=0):
        import mpi_submit as mpi
        import tempfile 
        if outfile == None:
            outf = tempfile.NamedTemporaryFile(dir=os.getcwd(), prefix="temp.calls.", delete=False)
            outfile= outf.name
        
        arg_dict = {"out": outfile, "conifer_file":self.h5file, "n_retry":n_retry}
        if verbose:
            arg_dict["verbose"] == "" # just need --verbose, nothing else

        if samples != None:
            try:
                arg_dict["samples"] = " ".join(samples)
            except TypeError:
                arg_dict["samples"] = samples
        if chromosomes != None:
            try:
                arg_dict["chromosomes"] = " ".join([str(c) for c in chromosomes])
            except TypeError:
                arg_dict["chromosomes"] = chromosomes
        
        arg_string = " ".join(["--%s %s" % (k,v) for k,v in arg_dict.items()])
        
        remote_command = "%s %s/call.mpi.py" % (sys.executable, os.path.dirname(os.path.realpath(__file__)))
        modules = ["modules", "modules-init", "modules-gs", "modules-eichler", "zlib/1.2.5", "hdf5/1.8.8", "bedtools/latest", "openmpi/1.5.3", "R/latest", "python/2.7.2", "numpy/1.6.1", "scipy/0.10.0", "lzo/2.06", "pytables/2.3.1_hdf5-1.8.8", "MySQLdb/1.2.3"]
        specs = "-S /bin/bash -cwd -j y -l h_vmem=2G -pe orte %d" % (n_cpus)
        cmd = remote_command + " " + arg_string
        print "makeCallsMPI qsub command:" , cmd
        mpi.submit_job(cmd,name="makeCallsMPI", modules = modules, native_specification=specs, mpirun=True, wait=True)
        
        outf.close()
        
        # read calls from tempfile and return
        calls = pandas.read_csv(outfile, header=0,sep="\t",index_col=0)
        os.remove(outfile)
        return CallTable(calls=calls)



if __name__ == "__main__":
    p = ConiferPipeline("/net/eichler/vol8/home/nkrumm/EXOMES/Autism_Main/SVD15_199.10bpminprobe.h5")
    import time 
    
    t1 = time.time()
    
    calls = p.makeCallsMPI(chromosomes=None, samples=None,n_cpus=351,n_retry=3)
    
    print "Done in ", time.time()-t1, "seconds"
    print calls.calls.head(100)
    calls.save("/net/eichler/vol8/home/nkrumm/EXOMES/CoNIFER/conifer_v0.3_DEV/calls_fulltest.alpha002undoSD.retry.txt")
