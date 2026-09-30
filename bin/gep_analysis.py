#!/usr/bin/env python3
"""Downstream analysis of one cNMF run: merge GEP usages into the AnnData, top genes, GOBP + CollecTRI enrichment,
and the summary figures. Outputs go to <out>/{csvs,excel,figures,anndatas}/.

Example (test sample):
    gep_analysis.py --adata results_test/ctype_MEL/anndatas/adata_topometry_ctype_MEL.h5ad \\
        --run_dir results_test/ctype_MEL/cnmf_test_MEL --run_name cnmf_test_MEL --k 3 --threshold 2 \\
        --sample ctype_MEL --out results_test/ctype_MEL/analysis
"""
import argparse
import sys
import warnings
import time
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # so `gep_lib` imports work wherever this is run from

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
matplotlib.use('Agg')
warnings.filterwarnings('ignore', category=FutureWarning, module='scanpy')   # scanpy's own deprecation notices

import gep_lib as gl
from gep_lib import figures


@contextmanager
def stage(name):
    """Print how long each stage takes (used to size the cluster request)."""
    print(f'--- {name}', flush=True)
    start = time.time()
    yield
    print(f'    {name}: {time.time() - start:.1f} s', flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group('inputs')
    g.add_argument('--adata', required=True, help='AnnData with the projection and metadata (topometry output)')
    g.add_argument('--run_dir', required=True, help='cNMF run folder, e.g. results/ctype_MEL/cnmf_run')
    g.add_argument('--run_name', required=True, help='cNMF run name (prefix of the files in --run_dir)')
    g.add_argument('--k', type=int, required=True, help='number of components used for consensus')
    g.add_argument('--threshold', type=float, required=True, help='local density threshold used for consensus')
    g.add_argument('--gmt', required=True, help='GO biological process GMT (see bin/fetch_resources.py)')
    g.add_argument('--collectri', required=True, help='CollecTRI network parquet (see bin/fetch_resources.py)')
    g.add_argument('--sample', default='sample', help='sample name, used in the output AnnData file name')
    g.add_argument('--out', required=True, help='output folder (csvs/, excel/, figures/, anndatas/ are made inside)')

    o = p.add_argument_group('options')
    o.add_argument('--projection', default='projection', help='obsm key of the embedding to plot on')
    o.add_argument('--clusters', default='clusters', help='obs column with cluster labels')
    o.add_argument('--age', default='age_pretty', help='obs column with age labels')
    o.add_argument('--summary_cols', default='age_pretty,clusters,ctype_detailed',
                   help='comma-separated obs columns for the projection summary; missing ones are skipped')
    o.add_argument('--usage_cutoff', type=float, default=0.1, help='usage above which a cell counts as using a GEP')
    o.add_argument('--min_cells', type=int, default=10, help='smallest cluster x age group shown in the heatmap')
    o.add_argument('--n_top_genes', type=int, default=100, help='genes saved per GEP')
    o.add_argument('--target_sum', type=float, default=1e4, help='normalize_total target for the gene panels')
    o.add_argument('--top_terms', type=int, default=10, help='gene sets shown per GEP in the bar plots')
    o.add_argument('--pval', type=float, default=0.05, help='keep gene sets with p below this')
    o.add_argument('--geps_per_page', type=int, default=6, help='GEPs per page in the multi-GEP figures (0 = one page)')
    o.add_argument('--formats', default='png,pdf', help='comma-separated figure formats')
    o.add_argument('--dpi', type=int, default=300)
    o.add_argument('--no_elements', action='store_true', help='skip the single-panel element figures')
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.out)
    csvs, excel, figs, anndatas = (out / d for d in ('csvs', 'excel', 'figures', 'anndatas'))
    for d in (csvs, excel, figs, anndatas):
        d.mkdir(parents=True, exist_ok=True)
    formats = [f.strip() for f in args.formats.split(',') if f.strip()]
    gl.set_publication_style(args.dpi)
    t0 = time.time()

    with stage('load and merge'):
        adata = ad.read_h5ad(args.adata)
        usage, scores = gl.load_gep_results(args.run_dir, args.run_name, args.k, args.threshold)
        gl.merge_gep_results(adata, usage, scores, args.run_name, args.k, args.threshold)
        geps = [gl.gep_number(c) for c in gl.gep_columns(adata)]

    with stage('top genes'):
        top = gl.top_genes(scores, args.n_top_genes)
        top.to_csv(csvs / f'top{args.n_top_genes}_genes.csv', index=False)
        gl.write_table_workbook(top, excel / f'top{args.n_top_genes}_genes.xlsx')

    results = {}
    with stage('enrichment'):
        info = {'method': 'decoupler ULM', 'p-value cut-off': args.pval, 'GEP scores from': args.run_dir,
                'k': args.k, 'density threshold': args.threshold}
        for name, net, note in (('gobp', gl.load_gobp(args.gmt), args.gmt), ('collectri', gl.load_collectri(args.collectri), args.collectri)):
            results[name] = gl.run_ulm(scores, net, args.pval)
            results[name].to_csv(csvs / f'gsea_{name}.csv', index=False)
            gl.write_gep_workbook(results[name], excel / f'gsea_{name}.xlsx', {**info, 'gene sets': note, 'name': name}, len(geps))
            print(f'    {name}: {len(results[name]):,} significant gene sets across {results[name]["gep"].nunique()} GEPs')

    with stage('usage summaries'):
        other = [c for c in args.summary_cols.split(',') if c not in (args.age, args.clusters)]
        palettes = gl.group_palettes(adata, args.age, args.clusters, other)
        gl.plot_usage_by_cluster(adata, args.clusters, palettes, args.usage_cutoff, figs / 'gep_usage_pct_by_cluster', formats)
        gl.plot_usage_by_age(adata, args.age, palettes, args.usage_cutoff, figs / 'gep_usage_pct_by_age', formats)
        gl.plot_usage_by_cluster_and_age(adata, args.clusters, args.age, palettes, args.usage_cutoff,
                                         figs / 'gep_usage_heatmap_cluster_age', formats, args.min_cells)

    with stage('projection summary'):
        gl.plot_projection_summary(adata, args.projection, [c for c in args.summary_cols.split(',') if c],
                                   palettes, args.clusters, figs / 'projection_summary', formats)

    with stage('gene expression for the top genes'):
        # scanpy's standard normalize_total + log1p (no scaling), from the raw counts of ALL genes, for the plotted genes only
        genes = list(dict.fromkeys(g for gep in top.columns for g in top[gep].iloc[:3]))
        total = np.asarray(adata.raw.X.sum(axis=1)).ravel()                  # each cell's counts over all genes
        counts = adata.raw[:, genes].X
        counts = counts.toarray() if sp.issparse(counts) else counts
        expr = ad.AnnData(np.log1p(counts / total[:, None] * args.target_sum), obs=adata.obs,
                          var=pd.DataFrame(index=genes), obsm={args.projection: adata.obsm[args.projection]})

    pages = gl.page_ranges(geps, args.geps_per_page)
    layouts = {'gep_megaplot': gl.MEGAPLOT, 'gep_usage_topgenes': gl.USAGE_GENES,
               'gep_gobp': ('gobp',), 'gep_collectri': ('collectri',)}
    for stem, columns in layouts.items():
        with stage(stem):
            for page in pages:
                gl.plot_gep_grid(adata, expr, page, columns, top, results, args.projection,
                                 figs / f'{stem}_{gl.page_label(page)}', formats, args.top_terms)

    if not args.no_elements:
        with stage('single-panel elements'):
            for gep in geps:
                for kind in figures.COLUMN_KINDS:
                    gl.plot_gep_element(adata, expr, gep, kind, top, results, args.projection,
                                        figs / 'elements' / f'gep_{gep:02d}_{kind}', formats, args.top_terms)

    with stage('save AnnData'):
        adata.write_h5ad(anndatas / f'adata_cnmf_{args.sample}.h5ad')

    print(f'Done in {time.time() - t0:.1f} s. Outputs in {out}')


if __name__ == '__main__':
    main()
