import pandas
pandas.set_option("display.line_width", 200)
from conifertools import *
import argparse
import numpy as np
import scipy.stats


def QC_Sample_Plot(calls, n_samples=None, title=None, outfile=None, outfile_list=None):
    g = calls.calls[["sampleID"]].groupby("sampleID", as_index=False)
    counts = np.sort(g.size().values)
    if n_samples is None:
        n_samples = len(counts)
    else:
        counts = np.hstack([np.zeros(n_samples-len(counts)), counts])
    fig, ax = plt.subplots(figsize=(8, 5))
    th = scipy.stats.scoreatpercentile(counts, 97.5)
    sample_cutoff = np.min(np.where(counts >= th)[0])

    ax.plot(counts[0:sample_cutoff], c='b',
            label="%d samples pass" % sample_cutoff)
    ax.scatter(np.where(counts >= th), counts[counts >= th],
               marker='x', lw=0.5, c='r', s=20,
               label="%d samples filtered" % (len(counts) - sample_cutoff))
    ax.axhline(th, c='r')
    ax.axvline(sample_cutoff, c="r")
    ax.annotate("Cutoff = %d calls/sample" % th,
                xy=[0.05, th], xycoords=(ax.transAxes, ax.transData))
    ax.legend(loc='lower left')
    ax.set_xlim([0, 1.1*len(counts)])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.set_xlabel("Sample Index (sorted by # of calls)")
    ax.set_ylabel("# of calls [Log scale]")
    if title is None:
        ax.set_title("Number of calls per sample")
    ax.set_yscale("symlog")
    if outfile is not None:
        plt.savefig(outfile)
    if outfile_list is not None:
        df = pandas.DataFrame(g.size(), columns=["call_count"])
        df = df.sort("call_count")
        df["threshold"] = th
        df["pass_threshold"] = df["call_count"] <= th
        df.to_csv(outfile_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", "-i", "--infiles",
                        nargs="+", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    parser.add_argument("--outfile_list", action="store", required=False)
    parser.add_argument("--title", required=False, default=None)
    parser.add_argument("--n_samples", type=int, required=False, default=None)
    args = parser.parse_args()

    calls = CallTable(args.infile)

    QC_Sample_Plot(calls,
                   n_samples=args.n_samples,
                   title=args.title,
                   outfile=args.outfile,
                   outfile_list=args.outfile_list)
