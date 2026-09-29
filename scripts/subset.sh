#!/bin/bash -l

# qsub scripts/subset.sh -c anndatas/adata_raw.h5ad -m ctype -x HC -o anndatas
# wrapper script to subset adata for cNMF on UCL Myriad

# wallclock time, 10 mins
#$ -l h_rt=0:10:0

# RAM, 8G
#$ -l mem=8G

# cores, 4
#$ -pe smp 4

# TMPDIR space
# 10GB is default

# job name
#$ -N cnmf-subset

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

# run script
/usr/bin/time --verbose apptainer run envs/cnmf_env.sif python scripts/subset.py \
	--counts $COUNTS \
	--metadata $METADATA \
	--subset $SUBSET \
	--output-dir $OUTDIR \

