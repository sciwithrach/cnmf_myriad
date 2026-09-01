#!/bin/bash -l

# qsub scripts/factorize.sh -o results -g hvgs.csv -i 200 -k 20 30 40
# wrapper script to run cNMF factorisation on UCL Myriad

# email, start and end
#$ -m be
#$ -M sjjgrww@ucl.ac.uk

# wallclock time, each iteration takes ~2m
#$ -l h_rt=5:00:0

# RAM, 2G
#$ -l mem=2G

# TMPDIR space
# 10GB is default

# job name
#$ -N cnmf-factorize

# array
#$ -t 1-600:50

# working directory
#$ -wd /home/sjjgrww/Scratch/cnmf

# run the script
# $SGE_TASK_ID = array number
# $JOB_NAME = job name

##############
# run script #
##############

# load helpers
source scripts/include.sh

# echo arguments
echo "Run name:                 $RUN_NAME"
echo "Path to output directory: $OUTDIR"
echo "Worker index:             $SGE_TASK_ID"
echo "Number of iterations:     $N_ITERS"
echo "Number of jobs for array: $N_JOBS"

# print the full call
set -x

# factorize, run in strides of 100
for (( i=$SGE_TASK_ID; i<$SGE_TASK_ID+50; i++ ))
do
	/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf factorize \
	    --output-dir $OUTDIR \
	    --name $RUN_NAME \
	    --worker-index $i \
	    --total-workers $N_JOBS \
	    --skip-completed-runs
done
