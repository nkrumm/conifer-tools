from conifertools import ConiferPipeline, CallTable, CallFilterTemplate
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--conifer_file", action="store", required=True)
    parser.add_argument("--call_file", action="store", required=True)
    parser.add_argument("--outfile", "-o", action="store", required=True)
    args = parser.parse_args()
    p = ConiferPipeline(args.conifer_file)

    SDFilter = CallFilterTemplate(p,
                "/net/eichler/vol8/home/nkrumm/REFERENCE_INFO/hg19genomicSuperDups.bed",
                name="SegDupProbeOverlap",
                filter_type="overlap",
                func=lambda x: x < 0.5)

    SDCount = CallFilterTemplate(p,
                "/net/eichler/vol8/home/nkrumm/REFERENCE_INFO/hg19genomicSuperDups.bed",
                name="SegDupProbeCount",
                filter_type="count")

    PPGFilter = CallFilterTemplate(p,
                "/net/eichler/vol8/home/nkrumm/REFERENCE_INFO/pp_genes.hg19.spans.bed",
                name="PPG_overlap",
                filter_type="overlap",
                func=lambda x: x < 0.1)

    PPG_probe_count = CallFilterTemplate(p,
                      "/net/eichler/vol8/home/nkrumm/REFERENCE_INFO/pp_genes.hg19.spans.bed",
                      name="PPG_exoncount",
                      filter_type="count")

    OtherDupFilter = CallFilterTemplate(p,
                     "/net/eichler/vol8/home/nkrumm/REFERENCE_INFO/3copiesin27of34.bed",
                     name="Dup_overlap",
                     filter_type="overlap",
                     func=lambda x: x < 0.5)

    def signalFilter(x):
        if x["num_probes"] <= 2:
            return np.abs(x["median_svdzrpkm"]) >= 1.5
        elif x["num_probes"] <= 5:
            return np.abs(x["median_svdzrpkm"]) >= 1
        else:
            return np.abs(x["median_svdzrpkm"]) >= 0.5

    calls = CallTable(args.call_file)
    calls = calls.filter(signalFilter)\
                 .filter(lambda x: x["probability"] > 0.99)\
                 .filter(SDFilter)\
                 .filter(PPGFilter)\
                 .filter(OtherDupFilter)\
                 .annotate(SDCount)\
                 .annotate(PPG_probe_count)

    calls.save(args.outfile)
