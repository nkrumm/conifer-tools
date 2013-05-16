from conifertools import CallTable
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", "-i", action="store", required=True)
    parser.add_argument("--cohort", action="store", required=True)
    parser.add_argument("--esp_infile", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    parser.add_argument("--gamma", type=float, action="store", required=False, default=0.9)
    parser.add_argument("--cophenetic_cutoff", type=float, action="store", required=False, default=0.85)

    args = parser.parse_args()
    assert args.gamma <= 1, "Gamma must be <= 1.00"
    assert args.cophenetic_cutoff <= 1, "Cophenetic cutoffs must be <= 1.00"

    calls = CallTable(args.infile)
    esp_calls = CallTable(args.esp_infile)

    calls.calls["cohort"] = args.cohort
    esp_calls.calls["cohort"] = "ESP"

    calls.appendCalls(esp_calls)

    calls = calls.clusterCallsByCohort(gamma=args.gamma, cohort_field="cohort", cophenetic_cutoff=args.cophenetic_cutoff)

    #clean up calls table
    del calls.calls["cnvrID_ESP"]
    calls = CallTable(calls.calls.rename(columns={'cnvr_frequency_HSCR': 'cnvr_frequency', 'cnvrID_HSCR': 'cnvrID'}))
    # filter for original cohort and save
    calls.filter(lambda x: x["cohort"] == args.cohort).save(args.outfile)
