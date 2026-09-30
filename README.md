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

# for the downstream GEP analysis (scanpy, decoupler, seaborn, openpyxl)
apptainer build envs/cnmf-analysis.sif envs/cnmf-analysis.def
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
| `n_iters` | no | iterations per k, e.g. `200`. Default: `--n_iters` |
| `hvgs` | if no `metadata`/`subset` | gene list, one per line |
| `run_name` | no | unique run name. Default: `cnmf_{tag}_{YYYYMMDD}` |

```csv
counts,metadata,subset,k_vals,n_iters,hvgs,run_name
anndatas/adata_ctype.h5ad,ctype,MEL,"20 30 40",200,,
anndatas/adata_ctype.h5ad,ctype,HC,"20 30 40 50",200,,
anndatas/adata_raw.h5ad,,,"20 30",,hvgs.csv,
```

- Rows with `metadata` and `subset`: the counts are subset (failing if the subset has 0 cells), topometry selects the HVGs, then prep onwards. `tag` is `{metadata}_{subset}`.
- Rows without: prep onwards, using the given `counts` and `hvgs`. `tag` is the counts file name.

### Run

Run from the `cnmf` directory, in `tmux` or `screen` (the Nextflow process has to stay alive while it submits jobs), with Java loaded and `nextflow` on your `PATH` (see Setup).

```bash
# real run
nextflow run main.nf -profile ucl_myriad --samplesheet samplesheet.csv

# after a failure or an edit: reuse cached tasks
nextflow run main.nf -profile ucl_myriad --samplesheet samplesheet.csv -resume

# check the setup first with a small run (assets/samplesheet_test.csv, results_test/)
nextflow run main.nf -profile ucl_myriad,test
```

- `-resume` continues the **latest** run in this directory. Without it, a run starts from scratch. Every `nextflow` command that starts a run counts as the latest run, including `-preview` and test-profile runs, so a bare `-resume` after one of those finds nothing cached and reruns everything.
- To continue a specific run, give its **session ID**, not its run name: `-resume 1f50da21-e19b-452b-ba34-87c375069c28`. `-resume` only takes `last` or a full session ID (this is how the launcher parses it in Nextflow 26.04): a run name is silently ignored, and the bare `-resume` resumes the latest run instead. Find the ID in the second-to-last column of `nextflow log` (for example `nextflow log | grep <run name>`); `nextflow log <run name> -f session` works too, but not while that session is running. A double dash (`--resume`) is a pipeline parameter and does nothing.
- Do not resume a run while one of its SGE jobs is still running (`qstat`), or a second copy of that task starts in the same work directory.
- `conf/base.config` (resources) is always loaded; `-profile ucl_myriad` adds the SGE executor and singularity settings.
- `scripts/run_nextflow.sh` is an optional wrapper that loads Java, writes the Nextflow log to `logs/nextflow.log`, and always adds `-profile ucl_myriad` and `-resume`. It can also be submitted as a job (`qsub scripts/run_nextflow.sh --samplesheet samplesheet.csv`) if compute nodes are allowed to submit jobs.

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | none | samplesheet CSV (required) |
| `--outdir` | `results_{run_date}` | one folder per sample (cNMF run, topometry, HVGs, logs) |
| `--run_date` | today (`YYYYMMDD`) | in the default outdir and run names. Pass the same value (or `--outdir`) when resuming on another day |
| `--k_vals` | `"20 30 40"` | k for rows with no `k_vals` |
| `--n_iters` | 100 | iterations per k, for rows with no `n_iters` |
| `--seed` | 42 | seed for `cnmf prepare` |
| `--stride` | 10 | factorize workers per SGE task |

### Cluster stability (clustree)

For rows that go through topometry, the pipeline also draws two plots per sample, into `{outdir}/{tag}/clustering/`. They only need the topometry AnnData, so they start as soon as topometry has finished (listed last in `main.nf`, so resuming an earlier run only adds them):

- `clustree_<series>.png`: cluster stability across the six Leiden resolutions (0.2 to 1.2), from `bin/clustree.R` (container `envs/clustree.sif`: `apptainer build envs/clustree.sif envs/clustree.def`).
- `cluster_comparison_on_embedding_<series>.png`: the projection coloured by `--cluster_color` (default `ctype`; skipped if the column is missing) and by every resolution, from `bin/cluster_compare.py`.

Both are made for every clustering series in the AnnData, so you can compare the latent spaces before choosing `--clusters` for the analysis step: `pca` (`pca_leiden_res*`), `topo` (`topo_clusters_res*`, spectral) and `topo_ms` (`topo_clusters_ms_res*`, multiscale). A series missing from the AnnData is skipped. Within a subset every cell has the same `ctype`, so that panel is one colour; use `--cluster_color clusters` to compare against the earlier cluster labels instead.

### After k plot inspection: consensus

Consensus is a separate step (`--step consensus`). Only one consensus step runs per run folder at a time (it writes a cache; the process holds a lock).

```bash
nextflow run main.nf -profile ucl_myriad --step consensus --outdir results_20260929 --sample ctype_MEL --selected_k 40 --threshold 0.5
```

The run name is read from `{outdir}/{sample}/pipeline_info/run_name.txt`; give `--run_name` to override it.

### After consensus: GEP analysis

`--step analysis` analyses one sample once consensus has been run: it merges the GEP usages into the topometry AnnData, saves the top genes per GEP, runs GO biological process and CollecTRI enrichment for every GEP (decoupler ULM, once for all GEPs), tests differential expression by cluster, and makes the figures.

One-off setup, on a machine with internet access (a login node). Compute nodes may not have any, so the analysis reads these files instead of downloading them:

```bash
singularity exec envs/cnmf-analysis.sif python bin/fetch_resources.py    # writes assets/resources/ (GMT + CollecTRI network)
```

```bash
nextflow run main.nf -profile ucl_myriad --step analysis --outdir results_20260929 --sample ctype_MEL --selected_k 40 --threshold 0.5 --clusters topo_clusters_ms_res0.6
```

`--clusters` is the obs column with the cluster labels to use. Choose it from the clustree plot (`clustering/clustree.png`), e.g. `topo_clusters_ms_res0.6`; the default is the `clusters` column. It is used for the usage-by-cluster plots, the projection summary and the differential expression, and the run stops at the start with the list of clustering columns if the name is not in the AnnData.

Use the same `--selected_k` and `--threshold` as for consensus. The run name is read from `run_name.txt` as for consensus. Options (defaults in `nextflow.config`): `--projection` (obsm key to plot on), `--age` (obs column), `--summary_cols` (columns for the projection summary; default age, clusters and `ctype_detailed`, and missing ones are skipped), `--usage_cutoff` (0.1), `--n_top_genes` (100), `--geps_per_page` (6; 0 for one page), `--formats` (`png,pdf`), `--gmt`, `--collectri`.

Outputs go to `{outdir}/{sample}/analysis/`:

```
csvs/       top100_genes.csv, gsea_gobp.csv, gsea_collectri.csv, dge_by_cluster.csv
excel/      top100_genes.xlsx, gsea_gobp.xlsx, gsea_collectri.xlsx   (gsea: one sheet per GEP)
figures/    gep_usage_pct_by_cluster, gep_usage_pct_by_age, gep_usage_heatmap_cluster_age, projection_summary,
            gep_megaplot_GEP01-06, gep_usage_topgenes_GEP01-06, gep_gobp_GEP01-06, gep_collectri_GEP01-06, ...
            elements/gep_01_usage, gep_01_gene1..3, gep_01_gobp, gep_01_collectri, ...   (every panel on its own)
anndatas/   adata_cnmf_{sample}.h5ad                                 (usages in obs, gene scores in varm)
analysis_call.txt, analysis_call.json                                  (how it was run: see below)
```

`analysis_call.txt` (to read) and `analysis_call.json` (to parse) record the call, so the settings behind any figure or CSV sit next to it: the `nextflow run` command, run name and session ID when it was run through `--step analysis`; the `gep_analysis.py` command with every option and default written out (and as typed); the input files with their paths, sizes and dates; and the versions of Python, scanpy, anndata, decoupler and the other packages. It is written at the start, so it is there even if the run fails, and a rerun for the same sample replaces it (as it replaces the other outputs). The notebook template writes the same record.

Figures are saved as PNG and PDF (points are rasterised, text stays editable). Gene panels and the differential expression use log-normalised expression from the raw counts of all genes: genes found in at least 3 cells are kept first, then `normalize_total` (to 1e4) and `log1p`, with no scaling. `dge_by_cluster.csv` has each cluster against all other cells (`rank_genes_groups`, Wilcoxon with `tie_correct=True`): group, gene, score, log fold change, p-value, adjusted p-value and the fraction of cells expressing it in the cluster and in the rest. Clusters with fewer than 10 cells are left out (`--dge_min_cells` on the command line script), and the step is skipped if fewer than two clusters remain.

The same steps are laid out one by one in `notebooks/gep_analysis_template.ipynb`, using the functions in `bin/gep_lib`, for running or changing any part in a notebook. The analysis step asks for 1 cpu, 8 GB and 1 h, scaled by attempt on a retry (`conf/base.config`). That is sized from runs on subsets of 7-13k cells (3-4 GB peak, 6-12 min, mostly the differential expression), so a much larger sample may need more; the run prints the time each stage takes.

### Outputs

- `results_{date}/{tag}/`: one folder per samplesheet row, holding:
  - `{run_name}/`: cNMF run folder
  - `topometry/`: plots and topometry object (subset rows)
  - `hvgs_{tag}.csv`: HVGs (subset rows)
  - `clustering/`: `clustree_{pca,topo,topo_ms}.png` and `cluster_comparison_on_embedding_{pca,topo,topo_ms}.png` (subset rows)
  - `logs/{PROCESS}.out|err`: this sample's task logs (successful tasks only; for a failed task see the `work/xx/yyyyyy/` folder Nextflow prints, `.command.out` and `.command.err`)
  - `pipeline_info/run_name.txt`: the cNMF run name, read by the consensus step
  - `pipeline_info/trace.txt`: this sample's trace rows, appended after every run (including consensus)
- `results_{date}/pipeline_info/trace.txt`: trace rows for all runs into this folder, appended after each; `trace_{timestamp}.txt` is Nextflow's raw trace for a single run
- `results_{date}/{tag}/anndatas/`: `adata_{tag}.h5ad` and `adata_topometry_{tag}.h5ad`, the subset and topometry AnnData (subset rows)
- `work/`: Nextflow working directory, safe to delete after a run finishes

### Layout

```
main.nf                 pipeline, consensus and analysis workflows (--step)
nextflow.config         parameters, plugins, profiles
nextflow_schema.json    parameter schema (nf-schema)
assets/                 samplesheet schema, test samplesheet, resources/ (GMT, CollecTRI; not in git)
notebooks/              gep_analysis_template.ipynb
bin/                    subset.py, topometry.py, cluster_compare.py, clustree.R, gep_analysis.py + gep_lib/, fetch_resources.py (on PATH inside tasks)
conf/                   base.config (resources), ucl_myriad.config, test.config
modules/local/          one process per step
envs/                   container images and definitions
scripts/                run_nextflow.sh (optional wrapper) and the older qsub scripts
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
