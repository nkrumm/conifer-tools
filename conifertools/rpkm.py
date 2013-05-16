import tables
import numpy as np


class rpkm_data:
    def __init__(self):
        self.rpkm = None
        self.samples = None
        self.exons = None
        self.isGenotype = False
        self.calls = []
        self.refined_calls = []

    def smooth(self, window=15, padded=False):  # todo, fix the padding here
        if self.isGenotype:
            print "Warning: the data in this rpkm_data container are single genotype values. Smoothing will have no effect!"
            return self.rpkm

        if window > 0:
            weightings = np.blackman(window)
            weightings = weightings/weightings.sum()
            smoothed_data = np.array([])
            for row in self.rpkm.transpose():
                smoothed = np.convolve(row, weightings)[(window-1)/2:-((window-1)/2)]
                if len(smoothed_data) == 0:
                    smoothed_data = smoothed
                else:
                    smoothed_data = np.vstack([smoothed_data, smoothed])

            self.rpkm = smoothed_data.transpose()
            return self.rpkm
        else:
            return self.rpkm

    def getSample(self, sampleIDs):
        sample_array = np.array(self.samples)
        if isinstance(sampleIDs, list):
            mask = np.zeros(len(sample_array), dtype=np.bool)
            for sampleID in sampleIDs:
                mask = np.logical_or(mask, sample_array == str(sampleID))

            return self.rpkm[:, mask]
        else:
            mask = sample_array == str(sampleID)
            return self.rpkm[:, mask]

    def getSamples(self, sampleIDs):
        return self.getSample(sampleIDs)

    @property
    def shape(self):
        if self.isGenotype:
            return [len(self.samples), 1]
        else:
            return [len(self.samples), len(self.exons)]


class rpkm_reader:
    def __init__(self, rpkm_fn=None):
        """Initialize an rpkm_reader instance. Specify the location of the data file"""

        if rpkm_fn is None:
            print "Must specify RPKM HDF5 file!"
            return 0
        # set up file access
        self.h5file = tables.openFile(rpkm_fn, mode='r')
        self.sample_table = self.h5file.root.samples.samples

    def __del__(self):
        self.h5file.close()

    def getExonValuesByRegion(self, chromosome, start=None, stop=None, exon_padding=0, sampleList=None, genotype=False, overlap=True):
        probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
        if (start is not None) and (stop is not None):
            if overlap:
                table_rows = probe_tbl.getWhereList('(stop >= %d) & (start <= %d)' % (start, stop))
            else:
                table_rows = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start, stop))
        else:
            table_rows = probe_tbl.getWhereList('(start >= 0) & (stop <= 1000000000)')

        data_tbl = self.h5file.root._f_getChild("chr" + str(chromosome))

        if sampleList is None:
            num_samples = data_tbl._v_nchildren
            samples = data_tbl
        else:
            num_samples = len(sampleList)
            samples = [data_tbl._f_getChild("sample_" + s) for s in sampleList]

        if exon_padding > 0:
            left_padding = max(table_rows[0]-exon_padding, 0)
            right_padding = min(table_rows[-1] + exon_padding, len(probe_tbl))
            table_rows = range(left_padding, right_padding)

        data = np.empty([num_samples, len(table_rows)], dtype=np.float)

        out_sample_list = []
        cnt = 0
        for sample_tbl in samples:
            d = sample_tbl.readCoordinates(table_rows, field="rpkm")
            data[cnt, :] = d
            cnt += 1
            out_sample_list.append(sample_tbl.title)

        d = rpkm_data()
        if genotype:  # return average #todo-- implement median and SD?
            d.rpkm = data.transpose().mean(axis=0)
            d.isGenotype = True
        else:  # return all data points
            d.rpkm = data.transpose()
        d.samples = out_sample_list
        d.exons = probe_tbl.readCoordinates(table_rows)
        d.contig = chromosome
        d.data_range = table_rows
        return d

    def getExonValuesByExons(self, chromosome, start_exon, stop_exon, sampleList=None, genotype=False):

        probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
        #table_rows = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start,stop))
        start_exon = max(start_exon, 0)
        stop_exon = min(stop_exon, probe_tbl.nrows)
        #print start_exon, stop_exon
        table_rows = np.arange(start_exon, stop_exon, 1)
        data_tbl = self.h5file.root._f_getChild("chr" + str(chromosome))

        if sampleList is None:
            num_samples = data_tbl._v_nchildren
            samples = data_tbl
        else:
            num_samples = len(sampleList)
            samples = [data_tbl._f_getChild("sample_" + s) for s in sampleList]

        data = np.empty([num_samples, len(table_rows)], dtype=np.float)

        out_sample_list = []
        cnt = 0
        for sample_tbl in samples:
            d = sample_tbl.readCoordinates(table_rows, field="rpkm")
            data[cnt, :] = d
            cnt += 1
            out_sample_list.append(sample_tbl.title)

        d = rpkm_data()
        if genotype:  # return average #todo-- implement median and SD?
            d.rpkm = data.transpose().mean(axis=0)
            d.isGenotype = True
        else:  # return all data points
            d.rpkm = data.transpose()
        d.samples = out_sample_list
        d.exons = probe_tbl.readCoordinates(table_rows)
        d.contig = chromosome
        return d

    def getChromosomeBySample(self, sampleID, chromosome, getexons=True):
        d = rpkm_data()
        data_tbl = self.h5file.root._f_getChild("chr" + str(chromosome))
        sample_tbl = data_tbl._f_getChild("sample_" + sampleID)
        d.rpkm = sample_tbl.read(field="rpkm")
        if getexons:
            probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
            d.exons = probe_tbl.read()
            d.contig = chromosome
        return d

    def getAllChromosomesBySample(self, sampleID):
        d = rpkm_data()
        for chr in range(1, 25):
            tbl = self.h5file.root._f_getChild("chr" + str(chr))._f_getChild("sample_" + sampleID)
            data = tbl.read(field="rpkm")
            if chr == 1:
                d.rpkm = data
            else:
                d.rpkm = np.hstack([d.rpkm, data])
        return d

    def getSampleList(self, cohort=None, sex=None, ethnicity=None, custom=None):
        """Return a list of available samples in the current data file. Specifying no arguments will return all available samples"""

        readWhereStr = ""
        if custom is not None:
            readWhereStr = custom
        else:
            if cohort is not None:
                if isinstance(cohort, list):
                    for c in cohort:
                        readWhereStr += "(cohort=='%s') | " % c
                    readWhereStr = readWhereStr.strip(" |")
                    readWhereStr += " & "
                else:
                    readWhereStr += "(cohort=='%s') " % cohort
            if sex is not None:
                if sex not in ['M', 'F']:
                    sex = sex.upper()[0]
                readWhereStr += " (sex=='%s') &" % sex
            if ethnicity is not None:
                readWhereStr += " (ethnicity=='%s') &" % ethnicity

            readWhereStr = readWhereStr.strip(" &")  # remove leading or trailing characters
        if readWhereStr != "":
            #print readWhereStr
            sampleIDs = self.sample_table.readWhere(readWhereStr, field='sampleID')
        else:
            sampleIDs = self.sample_table.read(field='sampleID')

        return sampleIDs

    def getSampleInfo(self, sampleID):
        return self.sample_table.readWhere("(sampleID=='%s')" % sampleID)

    def getExonIDs(self, chromosome, start=None, stop=None, anyIntersect=False):
        probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
        if anyIntersect:
            exons = probe_tbl.getWhereList('(stop >= %d) & (start <= %d)' % (start, stop))
        else:
            exons = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start, stop))
        return exons

    def getContigList(self):
        return self.h5file.root.probes._v_children.keys()
