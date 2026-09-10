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

## License

MIT
