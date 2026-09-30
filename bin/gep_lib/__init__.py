"""Helpers for downstream analysis of cNMF gene expression programmes (GEPs).

Used by bin/gep_analysis.py (the command line / Nextflow step) and by notebooks/gep_analysis_template.ipynb.
In a notebook: `sys.path.insert(0, 'bin')` then `from gep_lib import ...`.
"""
from .enrich import load_collectri, load_gobp, run_ulm
from .excel import write_gep_workbook, write_table_workbook
from .figures import (MEGAPLOT, USAGE_GENES, group_palettes, page_label, page_ranges, plot_gep_element, plot_gep_grid,
                      plot_projection_summary, plot_usage_by_age, plot_usage_by_cluster, plot_usage_by_cluster_and_age)
from .io import gep_columns, gep_number, load_gep_results, merge_gep_results, top_genes
from .record import write_call_record
from .style import set_publication_style, ucl_palette, usage_cmap

__all__ = [
    'load_gep_results', 'merge_gep_results', 'top_genes', 'gep_columns', 'gep_number',
    'load_gobp', 'load_collectri', 'run_ulm',
    'write_gep_workbook', 'write_table_workbook', 'write_call_record',
    'set_publication_style', 'ucl_palette', 'usage_cmap',
    'group_palettes', 'plot_usage_by_cluster', 'plot_usage_by_age', 'plot_usage_by_cluster_and_age',
    'plot_projection_summary', 'plot_gep_grid', 'plot_gep_element', 'page_ranges', 'page_label',
    'MEGAPLOT', 'USAGE_GENES',
]
