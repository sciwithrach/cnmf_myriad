"""Single plot panels. Each function draws on the `Axes` it is given, so the same code makes the small element figures
and the rows of the big multi-panel figures."""
from __future__ import annotations

import re
import textwrap

import anndata as ad
import matplotlib.axes
import numpy as np
import scanpy as sc

from .enrich import top_terms
from .style import DARKPUR


def _p99_or_max(values) -> float:
    """Colour limit: the 99th percentile, or the maximum when that is 0 (genes expressed in under 1% of the cells)."""
    top = np.nanpercentile(values, 99)
    return top if top > 0 else np.nanmax(values)


def _embedding(adata: ad.AnnData, basis: str, color: str, ax: matplotlib.axes.Axes, cmap, size, title: str, vmin=None) -> None:
    """One projection panel. Every panel is the same square, whatever the shape of the embedding: the axes are given a
    square box, and the colour bar is a small inset to the right instead of taking space from the axes."""
    sc.pl.embedding(adata, basis=basis, color=color, ax=ax, cmap=cmap, vmin=vmin, vmax=_p99_or_max, sort_order=True, size=size,
                    frameon=False, title=title, colorbar_loc=None, show=False)
    ax.set_box_aspect(1)
    cbar = ax.figure.colorbar(ax.collections[0], cax=ax.inset_axes([1.04, 0.15, 0.05, 0.7]))
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(labelsize=5, width=0.4, length=2)


def plot_usage(adata: ad.AnnData, gep_col: str, basis: str, ax: matplotlib.axes.Axes, cmap, size: float | None = None) -> None:
    """GEP usage on the projection (cells with the highest usage drawn on top, colour capped at the 99th percentile)."""
    _embedding(adata, basis, gep_col, ax, cmap, size, gep_col.replace('_', ' '))


def plot_gene(expr: ad.AnnData, gene: str, basis: str, ax: matplotlib.axes.Axes, cmap, size: float | None = None) -> bool:
    """Log-normalised expression of one gene on the projection. `expr` is a small AnnData of log-normalised counts (see bin/gep_analysis.py).

    Returns False, and leaves a note on the axes, if the gene is not in `expr`.
    """
    if gene not in expr.var_names:
        ax.axis('off')
        ax.text(0.5, 0.5, f'{gene}\nnot in data', ha='center', va='center', transform=ax.transAxes)
        return False
    _embedding(expr, basis, gene, ax, cmap, size, gene, vmin=0)
    return True


def id_first(term: str, width: int = 42) -> str:
    """'Long-Term Synaptic Depression (0060292)' -> '0060292: Long-Term Synaptic Depression'; only the name is shortened."""
    match = re.match(r'^(.*) \((\d+)\)$', term)
    if not match:
        return textwrap.shorten(term, width=width, placeholder='...')
    name, go_id = match.groups()
    return f'{go_id}: ' + textwrap.shorten(name, width=width - len(go_id) - 2, placeholder='...')


def plot_enrichment(results, gep: int, ax: matplotlib.axes.Axes, cmap, kind: str, top: int = 10,
                    title: str | None = None) -> None:
    """Bar plot of the strongest gene sets for one GEP.

    `results` is the long table from `run_ulm`; `kind` is 'gobp' (long names, drawn inside the bars like the notebook)
    or 'collectri' (short transcription factor names, drawn outside). `title` is written above the bars.
    """
    import decoupler as dc

    wide = top_terms(results, gep, top)
    if wide.empty:
        ax.axis('off')
        ax.text(0.5, 0.5, 'no significant\ngene sets', ha='center', va='center', transform=ax.transAxes)
        return

    if kind == 'gobp':   # 'ID: Name', with very long names shortened so they stay inside the panel
        wide.columns = [id_first(c) for c in wide.columns]
    dc.pl.barplot(wide, name=str(gep), top=top, cmap=cmap, ax=ax)

    ax.set_title(title or '', loc='left', fontsize=8)
    ax.grid(False)
    ax.set_xlabel('ULM score', fontsize=6)
    ax.tick_params(axis='x', labelsize=6, width=0.5)
    if kind == 'gobp':
        for label in ax.yaxis.get_majorticklabels():
            label.set_horizontalalignment('left')
            label.set_bbox(dict(facecolor='white', alpha=0.35, edgecolor='none', pad=0.4))   # readable over the bars
        ax.tick_params(axis='y', length=0, direction='in', labelsize=5, labelcolor=DARKPUR, pad=-2)
    else:
        ax.tick_params(axis='y', labelsize=5, width=0.5)     # 10 names in a ~22 mm tall panel
