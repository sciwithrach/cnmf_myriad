#!/bin/bash -l

# qsub scripts/combine_and_plot.sh -o results -n 20260827
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
#$ -N cnmf-combine

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

# combine
/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf combine \
--output-dir $OUTDIR \
--name $RUN_NAME \

# k selection plot
# run script
/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf k_selection_plot \
--output-dir $OUTDIR \
--name $RUN_NAME \

