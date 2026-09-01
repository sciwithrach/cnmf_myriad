Setup container for scripts:

# make envs subdirectory in cnmf directory
mkdir -p envs

# load apptainer
module load apptainer

# pull container and move to envs subdirectory
apptainer pull docker://quay.io/biocontainers/cnmf_1.7.1:pyhdfd78af_0 envs/cnmf_env.sif

