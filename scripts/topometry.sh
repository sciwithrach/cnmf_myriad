#!/bin/bash -l

# qsub scripts/topometry.sh -o results -c anndatas/adata_ctype_HC.h5ad -x HC -m ctype
# wrapper script to run topometry on subset adata for cNMF on UCL Myriad
# outputs: anndatas/adata_topometry_{metadata}_{subset}.h5ad, {outdir}/hvgs_{metadata}_{subset}.csv

# wallclock time, 12h
#$ -l h_rt=12:00:0

# RAM, 80G (20 * 4)
#$ -l mem=4G

# cores, 20
#$ -pe smp 20

# TMPDIR space
# 10GB is default
#$ -l tmpfs=20G

# job name
#$ -N cnmf-topometry

# working directory
#$ -wd /home/sjjgrww/Scratch/cnmf

# run the script
# $JOB_NAME = job name

##############
# run script #
##############

# load helpers
source scripts/include.sh

# echo arguments
echo "Path to counts:                   $COUNTS"
echo "Obs column:                       $METADATA"
echo "Subset:                           $SUBSET"
echo "Path to output directory:         $OUTDIR"

# work in $TMPDIR
START_DIR=$(pwd)
cd $TMPDIR

# run script
/usr/bin/time --verbose apptainer run $HOME/Scratch/cnmf/envs/utricle-qc.sif python $HOME/Scratch/cnmf/scripts/topometry.py \
	--adata "$HOME/Scratch/cnmf/$COUNTS" \
	--descriptor "$SUBSET" \
	--metadata "$METADATA"

# copy updated adata and HVGs to known paths for prep
TAG="${METADATA}_${SUBSET}"
mkdir -p "$START_DIR/anndatas" "$OUTDIR"
cp "adata_topometry_${SUBSET}.h5ad" "$START_DIR/anndatas/adata_topometry_${TAG}.h5ad"
cp "hvgs_${SUBSET}.csv" "$OUTDIR/hvgs_${TAG}.csv"

# copy plots and topometry object, $TMPDIR is deleted when the job ends
mkdir -p "$OUTDIR/topometry_${TAG}"
cp -r *.png *.pkl figures "$OUTDIR/topometry_${TAG}/"
