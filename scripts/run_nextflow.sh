#!/bin/bash -l

# bash scripts/run_nextflow.sh --samplesheet samplesheet.csv --n_iters 200
# starts the Nextflow head process, which submits every other job to SGE itself.
# Run it from the project root, in tmux/screen on a login node (or as a long, small qsub job if
# compute nodes may submit jobs), e.g.
#   mkdir -p logs && qsub -o results/logs -e results/logs scripts/run_nextflow.sh ...
# Extra arguments go to main.nf; use -entry consensus etc. as arguments.

#$ -l h_rt=48:0:0
#$ -l mem=2G
#$ -N cnmf-nextflow
#$ -wd /home/sjjgrww/Scratch/cnmf

# not a login shell under qsub
if ! type module &> /dev/null; then
	source /etc/profile.d/modules.sh
fi
module load java/temurin-17/17.0.2_8
export PATH=$PATH:$HOME/bin

mkdir -p logs

nextflow -log logs/nextflow.log run main.nf -profile ucl_myriad -resume "$@"
