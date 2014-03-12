from conifertools import ConiferPipeline
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", "-i", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    parser.add_argument("--chr", "-c", action="store", nargs="+", required=False, default="All")
    parser.add_argument("--samples", "-s", action="store", nargs="+", required=False, default="All")
    parser.add_argument("--families", "-f", action="store", nargs="+", required=False, default="All")
    parser.add_argument("--ncpu", "--ncpus", action="store", type=int, required=False, default=10)
    parser.add_argument("--nretry", "--nretries", action="store", type=int, required=False, default=3)
    parser.add_argument("--verbose", "-v", action="store_true", required=False, default=False)
    args = parser.parse_args()

    p = ConiferPipeline(args.infile)
    if args.chr == "All":
        chromosomes = None  # process all
        print "All Chromosomes"
    else:
        chromosomes = args.chr
        print "Chromosomes: %s" % ", ".join(chromosomes)
    
    if args.samples == "All":
        samples = None  # process all
        print "All samples"
    else:
        samples = args.samples
        print "Samples: ", ", ".join(samples)
    
    if args.families == "All":
        families=None
        print "Not calling by family"
    else:
        families = args.families
        print "Families:", ", ".join(families)

    calls = p.makeCallsMPI(chromosomes=chromosomes,
                           samples=samples,
                           families=families,
                           n_cpus=args.ncpu,
                           n_retry=args.nretry,
                           verbose=args.verbose)

    #print "Done Calling for chromosome(s): %s" % ", ".join(chromosomes)
    #print calls.calls.head(100)
    calls.save(args.outfile)
