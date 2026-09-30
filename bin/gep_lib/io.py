"""Read cNMF consensus results and merge them into the AnnData."""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import anndata as ad
import pandas as pd


def _dt(threshold: float) -> str:
    """cNMF writes the density threshold into file names as e.g. 0.2 -> '0_2', 2 -> '2_0'."""
    return str(float(threshold)).replace('.', '_')


def load_gep_results(run_dir: str | Path, run_name: str, k: int, threshold: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the consensus usages and gene scores written by `cnmf consensus`.

    Parameters
    ----------
    run_dir : folder holding the cNMF results, e.g. results_test/ctype_MEL/cnmf_test_MEL
    run_name : cNMF run name (the prefix of the files)
    k, threshold : the number of components and local density threshold given to `cnmf consensus`

    Returns
    -------
    usage : cells x GEPs, columns 'GEP_1'..'GEP_k', each row scaled to sum to 1 (as `cNMF.load_results` does)
    scores : genes x GEPs, columns 'GEP_1'..'GEP_k'; high values mean better marker genes
    """
    run_dir, dt = Path(run_dir), _dt(threshold)
    usage_file = run_dir / f'{run_name}.usages.k_{k}.dt_{dt}.consensus.txt'
    score_file = run_dir / f'{run_name}.gene_spectra_score.k_{k}.dt_{dt}.txt'
    for f in (usage_file, score_file):
        if not f.exists():
            raise FileNotFoundError(f'{f} not found. Has `cnmf consensus` been run for k={k}, threshold={threshold}?')

    usage = pd.read_csv(usage_file, sep='\t', index_col=0)
    usage = usage.div(usage.sum(axis=1), axis=0)
    usage.columns = [f'GEP_{int(c)}' for c in usage.columns]

    scores = pd.read_csv(score_file, sep='\t', index_col=0).T
    scores.columns = [f'GEP_{int(c)}' for c in scores.columns]
    return usage, scores


def merge_gep_results(adata: ad.AnnData, usage: pd.DataFrame, scores: pd.DataFrame, run_name: str, k: int,
                      threshold: float) -> ad.AnnData:
    """Add GEP usages to `adata.obs` and GEP gene scores to `adata.varm['gep_scores']` (in place).

    Usages are joined on the cell barcode; cells without a usage get NaN. Gene scores are aligned to `adata.var_names`
    (NaN where the gene was not in the cNMF run). Any earlier GEP_ columns are replaced.
    """
    old = [c for c in adata.obs.columns if c.startswith('GEP_')]
    adata.obs = adata.obs.drop(columns=old).join(usage, how='left')

    missing = usage.index.difference(adata.obs_names)
    if len(missing):
        print(f'WARNING: {len(missing)} cells in the cNMF usages are not in the adata')
    print(f'{adata.obs[usage.columns[0]].notna().sum():,} of {adata.n_obs:,} cells have GEP usages')

    adata.varm['gep_scores'] = scores.reindex(adata.var_names).to_numpy()
    adata.uns['gep_names'] = list(scores.columns)
    adata.uns['cnmf'] = {'run_name': run_name, 'k': int(k), 'threshold': float(threshold)}
    return adata


def top_genes(scores: pd.DataFrame, n: int = 100) -> pd.DataFrame:
    """Top `n` genes per GEP by score: an (n x GEPs) table, one column per GEP, best gene first."""
    return pd.DataFrame(
        {gep: scores[gep].sort_values(ascending=False).index[:n] for gep in scores.columns}
    )


def gep_number(gep: str) -> int:
    """'GEP_12' -> 12"""
    return int(str(gep).split('_')[-1])


def gep_columns(adata: ad.AnnData) -> list:
    """GEP usage columns in adata.obs, in numeric order."""
    return sorted([c for c in adata.obs.columns if c.startswith('GEP_')], key=gep_number)
