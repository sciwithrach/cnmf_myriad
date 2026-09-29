#!/bin/bash -l

# qsub scripts/pipeline.sh -c anndatas/adata_raw.h5ad -m ctype -x HC -k 20 30 40 -i 200 -o results
# wrapper script to submit the cNMF chain for one subset on UCL Myriad:
# subset -> topometry -> prep -> factorize -> combine (+ k selection plot)
# To skip subset and topometry, leave out -m and -x and give -g (HVGs):
# qsub scripts/pipeline.sh -c anndatas/adata_raw.h5ad -g hvgs.csv -k 20 30 40 -i 200 -o results
# This job only submits the other jobs (linked with -hold_jid), then exits.
# Add -d to print the qsub commands instead of submitting them.
#
# Safe to run for several subsets at once: every subset gets its own run name
# and output folder. Consensus is NOT run here: inspect the k selection plot,
# then run consensus separately (consensus.sh locks the run folder because
# consensus writes a cache there).

# wallclock time, 10 mins
#$ -l h_rt=0:10:0

# RAM, 1G
#$ -l mem=1G

# job name
#$ -N cnmf-pipeline

# working directory
#$ -wd /home/sjjgrww/Scratch/cnmf

# load helpers
source scripts/include.sh

# with -m and -x: subset -> topometry -> prep onwards
# without either: prep onwards, using -c (counts) and -g (HVGs) as given
if [[ -z $COUNTS ]]; then
	echo "Error: -c (counts) is required."
	exit 1
fi
if [[ -n $METADATA && -n $SUBSET ]]; then
	FROM_PREP=
	TAG="${METADATA}_${SUBSET}"
	SUBSET_H5AD="anndatas/adata_${TAG}.h5ad"
	PREP_COUNTS="anndatas/adata_topometry_${TAG}.h5ad"
	HVGS="$OUTDIR/hvgs_${TAG}.csv"
elif [[ -z $METADATA && -z $SUBSET ]]; then
	if [[ -z $HVG_PATH ]]; then
		echo "Error: -g (HVGs) is required when not subsetting."
		exit 1
	fi
	FROM_PREP=1
	TAG=$(basename "$COUNTS" .h5ad)
	PREP_COUNTS=$COUNTS
	HVGS=$HVG_PATH
else
	echo "Error: -m (obs column) and -x (subset) must be given together."
	exit 1
fi

# default run name is per-tag, unless -n was given
if [[ $RUN_NAME == "cnmf_$(date '+%Y%m%d')" ]]; then
	RUN_NAME="cnmf_${TAG}_$(date '+%Y%m%d')"
fi

# submit a job, print its id (fake id on dry run)
submit()
{
	if [[ -n $DRY_RUN ]]; then
		echo "qsub $*" >&2
		echo "DRY"
	else
		qsub -terse "$@" | cut -d. -f1
	fi
}

# hold flag, none for the first job
hold()
{
	[[ -n $1 ]] && echo "-hold_jid $1"
}

COMMON="-o $OUTDIR -n $RUN_NAME"

if [[ -z $FROM_PREP ]]; then
	SUBSET_ID=$(submit -N cnmf-subset-$TAG scripts/subset.sh -c $COUNTS -m $METADATA -x $SUBSET -o anndatas)

	TOPO_ID=$(submit -N cnmf-topometry-$TAG $(hold $SUBSET_ID) scripts/topometry.sh \
		-c $SUBSET_H5AD -m $METADATA -x $SUBSET -o $OUTDIR)
fi

# TOPO_ID is empty when starting from prep, so prep is not held
PREP_ID=$(submit -N cnmf-prep-$TAG $(hold $TOPO_ID) scripts/prep.sh \
	$COMMON -c $PREP_COUNTS -g $HVGS -s $SEED -i $N_ITERS -k $K_VALS)

FACT_ID=$(submit -N cnmf-factorize-$TAG $(hold $PREP_ID) -t 1-${N_JOBS}:50 scripts/factorize.sh \
	$COMMON -i $N_ITERS -k $K_VALS)

COMBINE_ID=$(submit -N cnmf-combine-$TAG $(hold $FACT_ID) scripts/combine_and_plot.sh \
	$COMMON)

echo "Submitted $RUN_NAME: subset ${SUBSET_ID:-skipped}, topometry ${TOPO_ID:-skipped}, prep $PREP_ID, factorize $FACT_ID, combine $COMBINE_ID"
