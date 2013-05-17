from conifertools import *
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--conifer_file", action="store", required=True)
    parser.add_argument("--call_file", action="store", required=True)
    parser.add_argument("--out_dir", action="store", required=True)
    args = parser.parse_args()

    calls = CallTable(args.call_file)

    plotter = ConiferPlotter(args.conifer_file)

    SDtrack = ConiferPlotTrack(plotter,
                               data_in="../../reference_data/hg19genomicSuperDups.bed",
                               name="SD",
                               color="y",
                               position=2.5)

    other_call_track = ConiferPlotTrack(plotter,
                                        data_in=calls,
                                        name="Calls",
                                        collapse=False,
                                        color="g",
                                        linewidth=4,
                                        position=2.3)

    plotter.add_track(SDtrack)
    plotter.add_track(other_call_track)

    for call in calls:
        plotter.basicPlot(call, outdir=args.out_dir)
