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

## Core script order

1. prep.sh

2. factorize.sh

3. consensus.sh

- Note: Only run one consensus step at a time to avoid re-writing the cache.

4. combine_and_plot.sh

## Subsetting & clustering data

Run these *before* running the core scripts above:

1. subset.sh

2. topometry.sh

## License

MIT
