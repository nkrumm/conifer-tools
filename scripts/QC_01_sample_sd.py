from conifertools import ConiferPipeline
import argparse
import numpy as np

if __name__ == "__main__":

parser = argparse.ArgumentParser()
parser.add_argument("--infile", "-i", action="store", required=True)
parser.add_argument("--outfile", "-o", action="store", required=True)
parser.add_argument("--plotfile", action="store", required=False, default=None)
parser.add_argument("--verbose", action="store_true", required=False, default=False)

p = ConiferPipeline(args.infile)

out_sd_values = []
with open(args.outfile, 'w') as out_file:
    for sampleID in p.samples:
        sample_data = []
        for c in xrange(1,24):
            sample_data.extend(p.getConiferData(sampleID, c).rpkm)
        sample_sd = np.std(np.array(sample_data))
        out_sd_values.append(sample_sd)
        if args.verbose:
            print sampleID, sample_sd
        out_file.write("%s\t%f\n" % (sampleID, sample_sd))

    if args.plotfile:
        from conifertools.plotting import QC_SampleSD_Plot
        kwargs = {"bins": 20, "color":'r'}
        QC_SampleSD_Plot(out_sd_values, title='SD values', outfile=args.plotfile, logscale=True, **kwargs)