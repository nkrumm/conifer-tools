from conifertools import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", "-i", nargs="+", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    parser.add_argument("--title", required=False, default=None)
    args = parser.parse_args()

    calls = CallTable(calls=args.infile)

    QC_Chromosome_Plot(calls, title=args.title, outfile=args.outfile)
