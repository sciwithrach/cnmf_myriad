#!/bin/bash -l

# qsub scripts/factorize.sh -o results -g hvgs.csv -i 200 -k 20 30 40
# wrapper script to run cNMF factorisation on UCL Myriad

# wallclock time, each iteration can take between 1m and 30m depending on k value
# (higher k is slower; longest observed was 26m). Size for the slowest case, as the job is
# killed at the limit: 10 iterations per task * 30m = 5h, plus margin = 6h.
# Killed tasks can be resubmitted (--skip-completed-runs), but if k or the data grow,
# raise this or lower STRIDE
#$ -l h_rt=6:00:0

# RAM, 2G
#$ -l mem=2G

# TMPDIR space
# 10GB is default

# job name
#$ -N cnmf-factorize

# array: default only (3 k x 200 iterations, 10 workers per task)
# directives can't use variables, so for other k/iterations override it:
# qsub -t 1-<number of k * iterations>:10 scripts/factorize.sh ...
# (pipeline.sh does this automatically). The step must equal STRIDE below.
#$ -t 1-600:10

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

# factorize, each task runs STRIDE workers (must match the step in -t)
STRIDE=10
for (( i=$SGE_TASK_ID; i<$SGE_TASK_ID+STRIDE && i<=N_JOBS; i++ ))
do
	/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf factorize \
	    --output-dir $OUTDIR \
	    --name $RUN_NAME \
	    --worker-index $i \
	    --total-workers $N_JOBS \
	    --skip-completed-runs
done
