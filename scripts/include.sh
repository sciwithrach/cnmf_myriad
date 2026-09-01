#!/bin/bash

# source include.sh
# helper script to run cNMF on UCL Myriad

########
# help #
########
Help()
{
    # display arg definitions
    echo "Wrapper script to run cNMF (Kotliar et al., 2025) on UCL Myriad."
    echo
    echo "Syntax:"
    echo "      bash prepare.sh    [-c -h|n|o|g|s|k|i]"
    echo "      bash factorize.sh  [-h|n|o]"
    echo
    echo "Options:"
    echo "h     Print help information      [all]"
    echo "n     Run name                    [all]       (default = {jobName}_{jobID})"
    echo "o     Path to output directory    [all]       (default = working directory)"
    echo "c     Path to count data          [prep]      (nb: required for cnmf prepare)"
    echo "g     Path to HVGs                [prep]      (default = none)"
    echo "s     Seed                        [prep]      (default = 42)"
    echo "k     List of components (k)      [factorize] (default = 20 30 40)"
    echo "i     Number of iterations to run [factorize] (default = 100)"
}

########
# args #
########

# initlialise
RUN_NAME="cnmf_$(date '+%Y%m%d')"
OUTDIR="$(pwd)/results"
SEED=42
K_VALS="20 30 40"
N_ITERS=100

# options
while getopts "hn:o:c:g:s:k:i:" option
do 
    case "${option}" in
        h)  
            Help
            exit;;
        n)
            RUN_NAME=${OPTARG};;
        o)
            OUTDIR=${OPTARG}
            mkdir -p $OUTDIR;;
        c)
            COUNTS=${OPTARG};;
        g)
            HVG_PATH=${OPTARG};;
        s)
            SEED=${OPTARG};;
        k)
            K_VALS=${OPTARG}
            # take very value until the next flag
            while [[ ${!OPTIND} != -* && -n ${!OPTIND} ]]; do
                K_VALS+=" ${!OPTIND}"
                OPTIND=$((OPTIND + 1))
            done;;
        i)
            N_ITERS=${OPTARG};;
        \?) 
            echo "Error: Invalid option -${OPTARG}."
            exit;;
    esac
done

###########
# modules #
###########

# load apptainer
module load apptainer
# add variable for anndata compatibility
export APPTAINERENV_ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1

####################################
# number of required jobs in array #
####################################

# number of values in k list
K_LENGTH=$(wc -w <<< "{$K_VALS}")

# number of k * number of iterations
N_JOBS=$(( $K_LENGTH*$N_ITERS ))

###########
# testing #
###########

#echo "Run name:                         $RUN_NAME"
#echo "Seed:                             $SEED"

#echo "Path to counts:                   $COUNTS"
#echo "Path to output directory:         $OUTDIR"
#echo "Path to HVGs:                     $HVG_PATH"

#echo "Components (k) for factorization: $K_VALS"
#echo "Length of k list:                 $K_LENGTH"
#echo "Number of iterations:             $N_ITERS"
#echo "Number of jobs for array:         $N_JOBS"
