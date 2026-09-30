"""Gene set enrichment of each GEP with decoupler's univariate linear model (ULM): GO biological process and CollecTRI."""
from __future__ import annotations

from itertools import chain, repeat
from pathlib import Path

import pandas as pd


def gmt_to_df(path: str | Path) -> pd.DataFrame:
    """Read a GMT file into the `source`/`target` (gene set, gene) table that decoupler expects."""
    pathways = {}
    with Path(path).open('r') as f:
        for line in f:
            name, _, *genes = line.strip().split('\t')
            pathways[name] = genes
    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=['source', 'target'],
    )


def tidy_go_names(net: pd.DataFrame) -> pd.DataFrame:
    """'REGULATION OF X%GOBP%GO:0000001' -> 'Regulation Of X (0000001)'."""
    net = net.copy()
    net['source'] = [x.title().replace('%Gobp%Go:', ' (') + ')' for x in net['source']]
    return net


def filter_gene_set_sizes(net: pd.DataFrame, min_size: int = 15, max_size: int = 500) -> pd.DataFrame:
    """Drop gene sets with <= `min_size` genes (unstable scores) or >= `max_size` genes (too general).

    The defaults match fgsea.
    """
    size = net.groupby('source').size()
    keep = size.index[(size > min_size) & (size < max_size)]
    print(f'{(size <= min_size).sum():,} small and {(size >= max_size).sum():,} large gene sets removed; '
          f'{len(keep):,} remain')
    return net[net['source'].isin(keep)]


def load_gobp(gmt: str | Path) -> pd.DataFrame:
    """GO biological process network (source, target) from a GMT file, tidied and size-filtered."""
    return filter_gene_set_sizes(tidy_go_names(gmt_to_df(gmt)))


def load_collectri(path: str | Path) -> pd.DataFrame:
    """CollecTRI transcription factor network (source, target, weight) saved by bin/fetch_resources.py."""
    return pd.read_parquet(path)


def run_ulm(scores: pd.DataFrame, net: pd.DataFrame, pval: float = 0.05) -> pd.DataFrame:
    """Score every gene set for every GEP at once, and keep the significant ones.

    Parameters
    ----------
    scores : genes x GEPs gene scores (columns 'GEP_1'...)
    net : gene sets as a source/target table (`load_gobp` or `load_collectri`)
    pval : keep gene sets with p < pval

    Returns a long table with columns gep, geneset, score, pval, ordered by GEP number then p-value.
    """
    import decoupler as dc

    mat = scores.dropna(axis='index').T          # GEPs x genes, as decoupler expects
    est, pvals = dc.mt.ulm(mat, net=net)

    long = pd.DataFrame({
        'gep': est.index.repeat(est.shape[1]),
        'geneset': list(est.columns) * est.shape[0],
        'score': est.to_numpy().ravel(),
        'pval': pvals.to_numpy().ravel(),
    })
    long = long[long['pval'] < pval].copy()
    long['gep'] = long['gep'].map(lambda g: int(str(g).split('_')[-1]))
    return long.sort_values(['gep', 'pval']).reset_index(drop=True)


def top_terms(results: pd.DataFrame, gep: int, n: int = 10) -> pd.DataFrame:
    """The `n` strongest gene sets (by absolute score) for one GEP, as a one-row frame ready for `dc.pl.barplot`."""
    sub = results[results['gep'] == gep]
    sub = sub.reindex(sub['score'].abs().sort_values(ascending=False).index).head(n)
    return sub.set_index('geneset')['score'].to_frame(name=str(gep)).T
