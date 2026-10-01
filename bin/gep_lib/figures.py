"""Whole figures: usage summaries, the projection summary, and the per-GEP grids and single-panel elements."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Sequence

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sb

from .io import gep_columns, gep_number
from .panels import plot_enrichment, plot_gene, plot_usage
from .style import (DOUBLE_COLUMN, MM, SINGLE_COLUMN, age_palette, categorical_palette, point_size, usage_cmap)

COLUMN_KINDS = ('usage', 'gene1', 'gene2', 'gene3', 'gobp', 'collectri')
# Figure layouts: which panels go in which figure
MEGAPLOT = ('usage', 'gene1', 'gene2', 'gene3', 'collectri')    # GOBP is too wide to read here; it has its own figure
USAGE_GENES = ('usage', 'gene1', 'gene2', 'gene3')


# --- helpers ----------------------------------------------------------------------------------------------------
def save_figure(fig, path_base: Path, formats: Sequence[str]) -> List[Path]:
    """Save `fig` as path_base.<format> for each format, then close it."""
    written = []
    for fmt in formats:
        path = Path(f'{path_base}.{fmt}')
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path)
        written.append(path)
    plt.close(fig)
    return written


_AGE = re.compile(r'^([EP])(\d+)$')   # embryonic (E14) before postnatal (P0, P4, P119)


def _sorted_categories(values: pd.Series) -> list:
    """Categories in natural order: numbers as numbers (2 before 10), ages E14 < P0 < P4 < P119 chronologically,
    anything else as stored (or alphabetical)."""
    cats = list(values.cat.categories) if hasattr(values, 'cat') else sorted(values.unique())
    present = set(values.astype(str))          # drop categories that no cell has (subsets keep the parent's list)
    cats = [c for c in cats if str(c) in present] or cats
    try:
        return sorted(cats, key=float)
    except (TypeError, ValueError):
        pass
    matches = [_AGE.match(str(c)) for c in cats]
    if all(matches):
        return [c for _, c in sorted(zip([(m.group(1) == 'P', int(m.group(2))) for m in matches], cats))]
    return cats


def group_palettes(adata: ad.AnnData, age_col: str, cluster_col: str, other_cols: Sequence[str] = ()) -> Dict[str, Dict]:
    """Colour dictionaries for age (the notebook's husl colours, in age order), clusters and any other categorical column
    (all colours different). Each palette is made once and also stored in `adata.uns['<column>_colors']` (scanpy's
    convention), so every figure, and any later notebook, shows the same colour for the same group."""
    palettes = {}
    for col in [age_col, cluster_col, *other_cols]:
        if col in palettes or col not in adata.obs:
            continue
        cats = [str(c) for c in _sorted_categories(adata.obs[col])]
        palettes[col] = age_palette(cats) if col == age_col else categorical_palette(cats)
        adata.obs[col] = pd.Categorical(adata.obs[col].astype(str), categories=cats)      # category order = colour order
        adata.uns[f'{col}_colors'] = [palettes[col][c] for c in cats]
    return palettes


def percent_above(adata: ad.AnnData, groupby, cutoff: float = 0.1) -> pd.DataFrame:
    """% of cells in each group whose usage of each GEP is above `cutoff` (groups x GEPs).

    `groupby` is an obs column name or a list of them.
    """
    cols = gep_columns(adata)
    keys = [groupby] if isinstance(groupby, str) else list(groupby)
    frame = sc.get.obs_df(adata, keys=keys + cols)
    for key in keys:
        frame[key] = frame[key].astype(str)
    pct = (frame[cols] > cutoff).groupby([frame[k] for k in keys], observed=True).mean() * 100
    pct.columns = [gep_number(c) for c in pct.columns]
    return pct


# --- usage summaries --------------------------------------------------------------------------------------------
def plot_usage_by_cluster(adata, cluster_col, palettes, cutoff, out_base, formats):
    """Clustered heatmap of % cells with GEP usage above `cutoff`, by cluster."""
    pct = percent_above(adata, cluster_col, cutoff)
    pct = pct.loc[[str(c) for c in _sorted_categories(adata.obs[cluster_col]) if str(c) in pct.index]]
    colours = pd.Series({c: palettes[cluster_col][c] for c in pct.index}, name='Cluster')
    g = sb.clustermap(pct, cmap='viridis', row_cluster=len(pct) > 1, col_cluster=pct.shape[1] > 1, row_colors=colours,
                      dendrogram_ratio=(0.06, 0.06), figsize=(DOUBLE_COLUMN, max(60 * MM, 4.5 * MM * len(pct) + 30 * MM)),
                      cbar_pos=(0.945, 0.3, 0.012, 0.4), linewidths=0)
    g.gs.update(right=0.9)      # leave room on the right for the colour bar, clear of the dendrogram and row colours
    g.cax.set_position([0.935, 0.3, 0.012, 0.4])
    g.cax.set_title(f'% cells with\nusage > {cutoff}', fontsize=6, loc='left')
    g.cax.tick_params(labelsize=5)
    g.ax_row_colors.tick_params(axis='x', length=0)
    g.ax_row_colors.grid(False)
    g.ax_heatmap.grid(False)
    g.ax_heatmap.set_xlabel('GEP')
    g.ax_heatmap.set_ylabel('Cluster')
    return save_figure(g.figure, out_base, formats), pct


def plot_usage_by_age(adata, age_col, palettes, cutoff, out_base, formats):
    """Heatmap of % cells with GEP usage above `cutoff`, by age."""
    pct = percent_above(adata, age_col, cutoff)
    pct = pct.loc[[str(c) for c in _sorted_categories(adata.obs[age_col]) if str(c) in pct.index]]
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN, max(45 * MM, 4.5 * MM * len(pct) + 25 * MM)), layout='constrained')
    sb.heatmap(pct, ax=ax, cmap='viridis', cbar_kws={'label': f'% cells with usage > {cutoff}', 'shrink': 0.6})
    ax.set_xlabel('GEP')
    ax.set_ylabel('Age')
    ax.grid(False)
    ax.tick_params(axis='y', rotation=0, length=0)
    ax.tick_params(axis='x', length=0)
    return save_figure(fig, out_base, formats), pct


def plot_usage_by_cluster_and_age(adata, cluster_col, age_col, palettes, cutoff, out_base, formats, min_cells=10):
    """Heatmap with one row per cluster x age group (sorted by cluster, then age) and side colour bars for both.

    Groups with fewer than `min_cells` cells are left out, since percentages from a handful of cells are unreliable.
    """
    obs = adata.obs
    counts = obs.groupby([obs[cluster_col].astype(str), obs[age_col].astype(str)], observed=True).size()
    pct = percent_above(adata, [cluster_col, age_col], cutoff)
    pct = pct.loc[counts[counts >= min_cells].index.intersection(pct.index)]

    cluster_order = {str(c): i for i, c in enumerate(_sorted_categories(obs[cluster_col]))}
    age_order = {str(c): i for i, c in enumerate(_sorted_categories(obs[age_col]))}
    pct = pct.loc[sorted(pct.index, key=lambda ix: (cluster_order[ix[0]], age_order[ix[1]]))]
    labels = [f'{c} | {a}' for c, a in pct.index]

    height = max(60 * MM, 2.4 * MM * len(pct) + 25 * MM)
    fig, axs = plt.subplots(1, 4, figsize=(DOUBLE_COLUMN, height), layout='constrained',
                            gridspec_kw={'width_ratios': [0.015, 0.015, 1, 0.025], 'wspace': 0.01})
    for ax, (name, col, pal) in zip(axs[:2], [('Cluster', 0, palettes[cluster_col]), ('Age', 1, palettes[age_col])]):
        colours = [mpl.colors.to_rgb(pal[ix[col]]) for ix in pct.index]
        ax.imshow(np.array(colours)[:, None, :], aspect='auto', interpolation='nearest')
        ax.set_xticks([0])
        ax.set_xticklabels([name], rotation=90)
        ax.tick_params(axis='x', length=0)
        ax.grid(False)
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    sb.heatmap(pct.reset_index(drop=True), ax=axs[2], cmap='viridis', cbar_ax=axs[3], yticklabels=labels,
               cbar_kws={'label': f'% cells with usage > {cutoff}'}, linewidths=0)
    axs[2].set_xlabel('GEP')
    axs[2].yaxis.tick_right()
    axs[2].tick_params(axis='y', labelsize=4, length=0, rotation=0)
    axs[2].tick_params(axis='x', length=0, rotation=0)
    axs[2].set_ylabel('')
    axs[2].grid(False)
    axs[3].grid(False)

    # a white line between clusters, across the colour bars and the heatmap
    starts = [i for i in range(1, len(pct)) if pct.index[i][0] != pct.index[i - 1][0]]
    for i in starts:
        for ax in axs[:2]:
            ax.axhline(i - 0.5, color='white', linewidth=1.5)    # imshow rows are centred on whole numbers
        axs[2].axhline(i, color='white', linewidth=1.5)          # heatmap cells span i to i + 1
    return save_figure(fig, out_base, formats), pct


# --- projection summary -----------------------------------------------------------------------------------------
def plot_projection_summary(adata, basis, columns, palettes, cluster_col, out_base, formats):
    """The projection coloured by each column in `columns` (side by side). Legends sit to the right of each panel,
    except clusters, which are labelled on the data. Columns missing from the adata are skipped with a warning."""
    present = []
    for col in columns:
        if col in adata.obs:
            present.append(col)
        else:
            print(f'WARNING: {col} not in adata.obs; skipped in the projection summary')
    if not present:
        return []

    fig, axs = plt.subplots(1, len(present), figsize=(len(present) * 75 * MM, 65 * MM), squeeze=False, layout='constrained')
    for ax, col in zip(axs[0], present):
        order = [str(c) for c in _sorted_categories(adata.obs[col])]
        adata.obs[col] = pd.Categorical(adata.obs[col].astype(str), categories=order)
        cats = order
        palette = palettes.get(col) or categorical_palette(cats)
        on_data = col == cluster_col
        sc.pl.embedding(adata, basis=basis, color=col, palette={c: palette[c] for c in cats if c in palette}, ax=ax,
                        frameon=False, title=col.replace('_', ' '), show=False, size=point_size(adata.n_obs, 75 * MM),
                        legend_loc='on data' if on_data else 'right margin', legend_fontsize=6 if on_data else 6,
                        legend_fontoutline=1.5 if on_data else None)
    return save_figure(fig, out_base, formats)


# --- per-GEP grids and elements ---------------------------------------------------------------------------------
def _draw_panel(kind, adata, expr, gep, topgenes, results, basis, cmap, ax, top, size, title=None):
    """Draw panel `kind` for one GEP on `ax`."""
    gep_col = f'GEP_{gep}'
    if kind == 'usage':
        plot_usage(adata, gep_col, basis, ax, cmap, size)
    elif kind.startswith('gene'):
        genes = list(topgenes[gep_col].iloc[:3])
        plot_gene(expr, genes[int(kind[-1]) - 1], basis, ax, cmap, size)
    elif kind in ('gobp', 'collectri'):
        plot_enrichment(results[kind], gep, ax, cmap, kind, top, title)
    else:
        raise ValueError(f'unknown panel {kind}')


def plot_gep_grid(adata, expr, geps, columns, topgenes, results, basis, out_base, formats, top=10):
    """One row per GEP, one column per panel kind in `columns`, on a single page.

    Embedding panels are square and bar-plot panels 1.8x wider. Wide layouts fit a 180 mm page, narrow ones 85 mm.
    """
    cmap = usage_cmap()
    ratios = [1.8 if c in ('gobp', 'collectri') else 1.0 for c in columns]
    width = DOUBLE_COLUMN if len(columns) > 2 else SINGLE_COLUMN
    unit = width / sum(ratios)
    has_embedding = any(c == 'usage' or c.startswith('gene') for c in columns)
    row_height = unit if has_embedding else max(unit, 36 * MM)   # square cells for projections; 10 GO bar labels need ~30 mm
    title = None if has_embedding else 'GEP {}'                  # the usage panel carries the GEP name; bar-only figures need it
    size = point_size(adata.n_obs, unit)
    fig, axs = plt.subplots(len(geps), len(columns), figsize=(width, row_height * len(geps)), squeeze=False,
                            gridspec_kw={'width_ratios': ratios}, layout='constrained')
    for i, gep in enumerate(geps):
        for j, kind in enumerate(columns):
            _draw_panel(kind, adata, expr, gep, topgenes, results, basis, cmap, axs[i, j], top, size,
                        title.format(gep) if title else None)
    return save_figure(fig, out_base, formats)


def plot_gep_element(adata, expr, gep, kind, topgenes, results, basis, out_base, formats, top=10):
    """A single panel of one GEP as its own figure (e.g. so it can be rearranged in a different layout later)."""
    fig, ax = plt.subplots(figsize=((SINGLE_COLUMN, 50 * MM) if kind in ('gobp', 'collectri') else (50 * MM, 45 * MM)),
                           layout='constrained')
    _draw_panel(kind, adata, expr, gep, topgenes, results, basis, usage_cmap(), ax, top, point_size(adata.n_obs, 50 * MM),
                f'GEP {gep}' if kind in ('gobp', 'collectri') else None)
    return save_figure(fig, out_base, formats)


def page_ranges(geps: Sequence[int], per_page: int) -> List[List[int]]:
    """Split GEP numbers into pages of `per_page` (0 or less means a single page)."""
    geps = list(geps)
    if per_page <= 0 or per_page >= len(geps):
        return [geps]
    return [geps[i:i + per_page] for i in range(0, len(geps), per_page)]


def page_label(page: Sequence[int]) -> str:
    """[1, 2, ..., 10] -> 'GEP01-10'"""
    return f'GEP{page[0]:02d}-{page[-1]:02d}'
