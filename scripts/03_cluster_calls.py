from conifertools import CallTable
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", "-i", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    parser.add_argument("--gamma", type=float, action="store", required=False, default=0.9)
    parser.add_argument("--cophenetic_cutoff", type=float, action="store", required=False, default=0.85)
    args = parser.parse_args()

    assert args.gamma <= 1, "Gamma must be <= 1.00"
    assert max(args.cophenetic_cutoff) <= 1, "All cophenetic cutoffs must be <= 1.00"

    calls = CallTable(args.infile)

    calls = calls.clusterCalls(gamma=args.gamma,
                               cophenetic_cutoff=args.cophenetic_cutoff)

    calls.save(args.outfile)
