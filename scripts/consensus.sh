#!/bin/bash -l

# qsub scripts/consensus.sh -o results -n 20260827 -f 40 -l 0.5
# wrapper script to run cNMF prep on UCL Myriad

# wallclock time, 10 mins
#$ -l h_rt=0:10:0

# RAM, 8G
#$ -l mem=8G

# cores, 4
#$ -pe smp 4

# TMPDIR space
# 10GB is default

# job name
#$ -N cnmf-consensus

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
echo "Run name:                         $RUN_NAME"
echo "Path to output directory:         $OUTDIR"

# consensus
/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf consensus \
	--output-dir $OUTDIR \
	--name $RUN_NAME \
	--components $SELECTED_K \
	--local-density-threshold $THRESHOLD \
	--show-clustering

