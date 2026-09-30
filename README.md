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

Install Nextflow (needs Java 11+; see the [nf-core Myriad guide](https://github.com/nf-core/configs/blob/master/docs/ucl_myriad.md)):

```bash
# add to ~/.bashrc
module load java/temurin-17/17.0.2_8

# download Nextflow into ~/bin
curl -s https://get.nextflow.io | bash
chmod a+x nextflow
mkdir -p ~/bin && mv nextflow ~/bin/

# add to ~/.bash_profile
export PATH=$PATH:$HOME/bin
```

## Running the pipeline (Nextflow)

Nextflow runs one chain per samplesheet row and submits each step to SGE itself: `subset -> topometry -> prep -> factorize -> combine (+ k selection plot)`. A failed step stops that row's later steps. **Run from the `cnmf` directory** and always give `-profile ucl_myriad` (the [nf-core Myriad config](https://github.com/nf-core/configs/blob/master/docs/ucl_myriad.md), in `conf/ucl_myriad.config`).

### Samplesheet

One row per run (`samplesheet.csv`). The columns are checked against `assets/schema_input.json` before any job is submitted.

| Column | Required | Description |
|---|---|---|
| `counts` | yes | AnnData (`.h5ad`) with raw counts |
| `metadata` | with `subset` | obs column to subset by, e.g. `ctype` |
| `subset` | with `metadata` | value to keep, e.g. `MEL` |
| `k_vals` | no | components (k), space or `;` separated, e.g. `"20 30 40"`. Default: `--k_vals` |
| `hvgs` | if no `metadata`/`subset` | gene list, one per line |
| `run_name` | no | unique run name. Default: `cnmf_{tag}_{YYYYMMDD}` |

```csv
counts,metadata,subset,k_vals,hvgs,run_name
anndatas/adata_ctype.h5ad,ctype,MEL,"20 30 40",,
anndatas/adata_ctype.h5ad,ctype,HC,"20 30 40 50",,
anndatas/adata_raw.h5ad,,,"20 30",hvgs.csv,
```

- Rows with `metadata` and `subset`: the counts are subset (failing if the subset has 0 cells), topometry selects the HVGs, then prep onwards. `tag` is `{metadata}_{subset}`.
- Rows without: prep onwards, using the given `counts` and `hvgs`. `tag` is the counts file name.

### Run

```bash
# from the cnmf directory, e.g. in tmux
bash scripts/run_nextflow.sh --samplesheet samplesheet.csv --n_iters 200

# or directly
nextflow run main.nf -profile ucl_myriad --samplesheet samplesheet.csv --n_iters 200

# check the setup first with a small run (assets/samplesheet_test.csv, results_test/)
nextflow run main.nf -profile ucl_myriad,test
```

`run_nextflow.sh` loads Java, writes the Nextflow log to `logs/nextflow.log` and adds `-resume`, so rerunning the same command after a failure only reruns the missing tasks.

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | none | samplesheet CSV (required) |
| `--outdir` | `results_{run_date}` | one folder per sample (cNMF run, topometry, HVGs, logs) |
| `--run_date` | today (`YYYYMMDD`) | in the default outdir and run names. Pass the same value (or `--outdir`) when resuming on another day |
| `--k_vals` | `"20 30 40"` | k for rows with no `k_vals` |
| `--n_iters` | 100 | iterations per k |
| `--seed` | 42 | seed for `cnmf prepare` |
| `--stride` | 10 | factorize workers per SGE task |

### After k plot inspection: consensus

Consensus is a separate entry point. Only one consensus step runs per run folder at a time (it writes a cache; the process holds a lock).

```bash
nextflow run main.nf -profile ucl_myriad -entry consensus --outdir results_20260929 --sample ctype_MEL --selected_k 40 --threshold 0.5
```

The run name is read from `{outdir}/{sample}/pipeline_info/run_name.txt`; give `--run_name` to override it.

### Outputs

- `results_{date}/{tag}/`: one folder per samplesheet row, holding:
  - `{run_name}/`: cNMF run folder
  - `topometry/`: plots and topometry object (subset rows)
  - `hvgs_{tag}.csv`: HVGs (subset rows)
  - `logs/{PROCESS}.out|err`: this sample's task logs (successful tasks only; for a failed task see the `work/xx/yyyyyy/` folder Nextflow prints, `.command.out` and `.command.err`)
  - `pipeline_info/run_name.txt`: the cNMF run name, read by the consensus step
  - `pipeline_info/trace.txt`: this sample's trace rows, appended after every run
- `results_{date}/pipeline_info/trace.txt`: trace rows for all runs into this folder, appended after each; `trace_{timestamp}.txt` is Nextflow's raw trace for a single run
- `results_{date}/{tag}/anndatas/`: `adata_{tag}.h5ad` and `adata_topometry_{tag}.h5ad`, the subset and topometry AnnData (subset rows)
- `work/`: Nextflow working directory, safe to delete after a run finishes

### Layout

```
main.nf                 workflow and consensus entry point
nextflow.config         parameters, plugins, profiles
nextflow_schema.json    parameter schema (nf-schema)
assets/                 samplesheet schema, test samplesheet
bin/                    subset.py, topometry.py (on PATH inside tasks)
conf/                   base.config (resources), ucl_myriad.config, test.config
modules/local/          one process per step
envs/                   container images and definitions
scripts/                run_nextflow.sh and the older qsub scripts
```

Resources (cpus, total memory, wallclock) are set per label in `conf/base.config`. The Myriad profile divides memory by cpus for the per-core request. Change them there, not in the modules.

## Running with qsub scripts (older route)

The scripts in `scripts/` still work. Each one is a job on its own, or `pipeline.sh` chains them.

## Order of jobs

### Subsetting & clustering data (optional)

1. `subset.sh`

2. `topometry.sh`

### Core cNMF functions

1. `prep.sh`

2. `factorize.sh`
   - Standalone, override the array size: `qsub -t 1-N:10 scripts/factorize.sh ...`, where N = number of k values x iterations (`pipeline.sh` does this for you). If tasks are killed at the wallclock limit, resubmit with the same command: `--skip-completed-runs` skips finished runs. Check with `qacct -j <jobid>` (`ru_wallclock`, `exit_status`).

4. `combine_and_plot.sh`

### After k plot inspection

1. `consensus.sh`

- Note: Only run one consensus step is run at a time to avoid re-writing the cache.

### pipeline.sh

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
