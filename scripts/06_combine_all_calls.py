from conifertools import CallTable
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--call_files", nargs="+", action="store", required=True)
    parser.add_argument("--outfile", action="store", required=True)
    parser.add_argument("--cols", nargs="+", action="store", default=[], required=False)
    args = parser.parse_args()

    calls = CallTable(args.call_files)

    if len(args.cols) > 0:
        calls.calls[args.cols]\
             .sort(["chromosome", "start"])\
             .to_csv(args.outfile, sep="\t")
    else:
        calls.calls\
             .sort(["chromosome", "start"])\
             .to_csv(args.outfile, sep="\t")
