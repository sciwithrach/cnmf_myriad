#!/bin/bash -l

# qsub scripts/prep.sh -c anndatas/adata_raw.h5ad -k 20 30 40 50 -o results -g hvgs.csv
# wrapper script to run cNMF prep on UCL Myriad

# email, start and end
#$ -m be
#$ -M sjjgrww@ucl.ac.uk

# wallclock time, 10 mins
#$ -l h_rt=1:0:0

# RAM, 8G
#$ -l mem=8G

# cores, 4
#$ -pe smp 4

# TMPDIR space
# 10GB is default

# job name
#$ -N cnmf-prep

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
echo "Seed:                             $SEED"
echo "Components (k) for factorization: $K_VALS"
echo "Path to counts:                   $COUNTS"
echo "Path to output directory:         $OUTDIR"
echo "Path to HVGs:                     $HVG_PATH"

# run script
/usr/bin/time --verbose apptainer run envs/cnmf_env.sif cnmf prepare \
--output-dir $OUTDIR \
--name $RUN_NAME \
-c $COUNTS \
-k $K_VALS \
--n-iter $N_ITERS \
--genes-file $HVG_PATH \
--seed $SEED \
--total-workers 4
