#!/bin/bash
##############################################################################
##############################################################################
# conifer-tools: Tools for using CoNIFER
# Developed by Niklas Krumm (C) 2013
# nkrumm@gmail.com
# 
# homepage: http://conifer.sf.net
# This program is described in:
# Krumm et al. 2012. Genome Research. doi:10.1101/gr.138115.112
#
# This file is part of conifer-tools: 
# https://github.com/nkrumm/conifer-tools
##############################################################################

##############################################################################
## REQUIRED OPTIONS

CONIFER_ANALYSIS_FILE="/path/to/conifer_file.hdf5"
N_SAMPLES=
PREFIX="" # prefix for file names
OUTDIR="/base/path/to/out_directory"
NCPUS=100 # for MPI use
CHRS=`seq 1 24` # adjust to process specific chromosomes

QCDIR="$OUTDIR/QC"
PLOTDIR="$OUTDIR/PLOTS"

##############################################################################
##############################################################################
## START PROCESSING
## Trap interrupts and exit instead of continuing the loop
trap "echo Exited!; exit;" SIGINT SIGTERM
mkdir -p $QCDIR
mkdir -p $PLOTDIR

for CHR in $CHRS; 
do
    echo "Calling for chromosome $CHR"
    OUTFILE=$OUTDIR/${PREFIX}.chr${CHR}.calls.csv
    python 01_create_calls.py --infile=$CONIFER_ANALYSIS_FILE --chr $CHR --ncpus $NCPUS --outfile=$OUTFILE
    
    echo "Filtering calls for chromosome $CHR"
    INFILE=$OUTFILE
    OUTFILE="${OUTDIR}/${PREFIX}.chr${CHR}.calls.filtered.csv"   
    python 02_filter_calls.py --conifer_file=$CONIFER_ANALYSIS_FILE --call_file=$INFILE --outfile=$OUTFILE  

    echo "Clustering calls for chromosome $CHR"
    INFILE=$OUTFILE
    OUTFILE="${OUTDIR}/${PREFIX}.chr${CHR}.calls.clustered.csv"
    python 03_cluster_calls.py --infile=$INFILE --outfile=$OUTFILE

    echo "QC plot for chromosome $CHR"
    python plotting/QC01_QCplot.py --infile="$OUTFILE" --outfile="${QCDIR}/${PREFIX}.chr${CHR}.filtered_calls.png" --title="QC Plot: Chr${CHR} (filtered calls)" 

    echo "Plotting all calls for chromosome $CHR"
    mkdir -p $PLOTDIR/$CHR
    python plotting/QC02_plot_cnvrs.py --conifer_file=$CONIFER_ANALYSIS_FILE --min_freq 2 --max_freq 30 --call_file=$OUTFILE --out_dir=$PLOTDIR/$CHR/
done

echo "QC Plot of # of calls per sample"
python plotting/QC03_sample_plot.py --infile ${OUTDIR}/${PREFIX}.chr*.calls.clustered.csv --outfile ${QCDIR}/QC.SamplePlot.png --outfile_list ${QCDIR}/call_counts_per_sample.csv --n_samples=$N_SAMPLES

echo "QC Plot of all calls"
python plotting/QC01_QCplot.py --infile ${OUTDIR}/${PREFIX}.chr*.calls.clustered.csv --outfile="${QCDIR}/${PREFIX}.all_chrs.filtered_calls.png" --title="QC Plot: All Chromosomes (filtered calls)"   

echo "Final combined call file for all chromosomes"
COLS='sampleID chromosome start stop state start_exon stop_exon num_probes size_bp median_svdzrpkm cnvrID cnvr_frequency'
python 06_combine_all_calls.py --cols $COLS --call_files ${OUTDIR}/${PREFIX}.chr*.calls.clustered.csv --outfile ${OUTDIR}/${PREFIX}.all_chrs.calls.clustered.final.csv



#echo "Filtering samples based on max # of calls cutoff determined in QC03..."
#python 04_filter_samples.py --samplefile_in ../${PREFIX}.samples.txt --count_file ${QCDIR}/call_counts_per_sample.csv --samplefile_out ${OUTDIR}/${PREFIX}.samples.passQC.txt
