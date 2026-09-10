#!/bin/bash -l

# qsub scripts/topometry.sh -c anndatas/adata_ctype_HC.h5ad -x HC -m ctype
# wrapper script to run topometry on subset adata for cNMF on UCL Myriad

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
echo "Run name:                         $RUN_NAME"
echo "Path to counts:                   $COUNTS"

# work in $TMPDIR
cd $TMPDIR

# run script
/usr/bin/time --verbose apptainer run $HOME/Scratch/cnmf/envs/utricle-qc.sif python $HOME/Scratch/cnmf/scripts/topometry.py \
	--adata $HOME/Scratch/cnmf/$COUNTS \
	--descriptor $SUBSET \
	--metadata $METADATA

# copy files
tar -zcvf $HOME/Scratch/cnmf/results/topometry_$SUBSET_$JOB_ID.tar.gz $TMPDIR

