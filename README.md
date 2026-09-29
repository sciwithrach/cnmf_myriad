# cnmf - lipovsek lab

## Setup

Setup containers as follows:

```bash
# make envs subdirectory in cnmf directory
mkdir -p envs

# load apptainer
module load apptainer

# pull container and move to envs subdirectory
apptainer pull docker://quay.io/biocontainers/cnmf_1.7.1:pyhdfd78af_0 envs/cnmf_env.sif

# for running topometry on subset of data
apptainer build envs/utricle-qc.sif envs/utricle-qc.def
```

## Order of jobs

### Subsetting & clustering data (optional)

1. `subset.sh`

2. `topometry.sh`

### Core cNMF functions

1. `prep.sh`

2. `factorize.sh`
   - Standalone, override the array size: `qsub -t 1-N:50 scripts/factorize.sh ...`, where N = number of k values x iterations (`pipeline.sh` does this for you).

4. `combine_and_plot.sh`

### After k plot inspection

1. `consensus.sh`

- Note: Only run one consensus step is run at a time to avoid re-writing the cache.

## Running the pipeline

`pipeline.sh` submits the whole chain as linked SGE jobs (`-hold_jid`). **Run from the `cnmf` directory.**

Consensus is not included: inspect the k selection plot, then run `consensus.sh` separately.

### Dry run (prints the qsub commands, submits nothing):

```bash
# with subsetting
bash scripts/pipeline.sh -d -c anndatas/adata_raw.h5ad -m ctype -x HC -k 20 30 40 -i 200 -o results

# without subsetting
bash scripts/pipeline.sh -d -c anndatas/adata_raw.h5ad -g hvgs.csv -k 20 30 40 -i 200 -o results
```

### Submit:

```bash
# with subsetting: subset -> topometry -> prep -> factorize -> combine
qsub scripts/pipeline.sh -c anndatas/adata_raw.h5ad -m ctype -x HC -k 20 30 40 -i 200 -o results

# without subsetting: prep -> factorize -> combine (-g is required)
qsub scripts/pipeline.sh -c anndatas/adata_raw.h5ad -g hvgs.csv -k 20 30 40 -i 200 -o results
```

Default run name is `cnmf_{tag}_{YYYYMMDD}` (tag is `{ctype}_{subset}` or the counts file name), unless `-n` is given. Outputs go in `results/{run name}/` unless `-o` is given.

## License

MIT
