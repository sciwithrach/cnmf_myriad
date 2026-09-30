"""Single plot panels. Each function draws on the `Axes` it is given, so the same code makes the small element figures
and the rows of the big multi-panel figures."""
from __future__ import annotations

import textwrap

import anndata as ad
import matplotlib.axes
import scanpy as sc

from .enrich import top_terms
from .style import DARKPUR


def plot_usage(adata: ad.AnnData, gep_col: str, basis: str, ax: matplotlib.axes.Axes, cmap, size: float | None = None) -> None:
    """GEP usage on the projection (cells with the highest usage drawn on top, colour capped at the 99th percentile)."""
    sc.pl.embedding(adata, basis=basis, color=gep_col, ax=ax, cmap=cmap, vmax='p99', sort_order=True, size=size,
                    frameon=False, title=gep_col.replace('_', ' '), show=False)


def plot_gene(expr: ad.AnnData, gene: str, basis: str, ax: matplotlib.axes.Axes, cmap, size: float | None = None) -> bool:
    """Log-normalised expression of one gene on the projection. `expr` is a small AnnData of log-normalised counts (see bin/gep_analysis.py).

    Returns False, and leaves a note on the axes, if the gene is not in `expr`.
    """
    if gene not in expr.var_names:
        ax.axis('off')
        ax.text(0.5, 0.5, f'{gene}\nnot in data', ha='center', va='center', transform=ax.transAxes)
        return False
    sc.pl.embedding(expr, basis=basis, color=gene, ax=ax, cmap=cmap, vmin=0, vmax='p99', sort_order=True, size=size,
                    frameon=False, title=gene, show=False)
    return True


def plot_enrichment(results, gep: int, ax: matplotlib.axes.Axes, cmap, kind: str, top: int = 10) -> None:
    """Bar plot of the strongest gene sets for one GEP.

    `results` is the long table from `run_ulm`; `kind` is 'gobp' (long names, drawn inside the bars like the notebook)
    or 'collectri' (short transcription factor names, drawn outside).
    """
    import decoupler as dc

    wide = top_terms(results, gep, top)
    if wide.empty:
        ax.axis('off')
        ax.text(0.5, 0.5, 'no significant\ngene sets', ha='center', va='center', transform=ax.transAxes)
        return

    if kind == 'gobp':   # shorten very long GO names so they stay inside the panel
        wide.columns = [textwrap.shorten(c, width=42, placeholder='...') for c in wide.columns]
    dc.pl.barplot(wide, name=str(gep), top=top, cmap=cmap, ax=ax)

    ax.set_title('')
    ax.grid(False)
    ax.set_xlabel('ULM score', fontsize=6)
    ax.tick_params(axis='x', labelsize=6, width=0.5)
    if kind == 'gobp':
        for label in ax.yaxis.get_majorticklabels():
            label.set_horizontalalignment('left')
            label.set_bbox(dict(facecolor='white', alpha=0.75, edgecolor='none', pad=0.4))   # readable over the bars
        ax.tick_params(axis='y', length=0, direction='in', labelsize=5, labelcolor=DARKPUR, pad=-2)
    else:
        ax.tick_params(axis='y', labelsize=6, width=0.5)
