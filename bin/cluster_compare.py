#!/usr/bin/env python3
"""Cluster comparison plots and clustree input for one topometry AnnData, for every clustering series it has:

    pca      pca_leiden_res*        PCA-based latent space
    topo     topo_clusters_res*     Topometry - spectral latent space
    topo_ms  topo_clusters_ms_res*  Topometry - multiscale latent space

Writes, in the current directory:
    cluster_comparison_on_embedding_<series>.png   the projection coloured by --color and by every resolution in the series
    cluster_table.csv                              cells x resolutions of all series, the input for bin/clustree.R
    clustree_series.txt                            one tab-separated line per series: name, column prefix, plot title
"""
import argparse

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc

# name, column prefix (followed by the resolution), title of the clustree plot
SERIES = [
    ('pca', 'pca_leiden_res', 'PCA-based latent space'),
    ('topo', 'topo_clusters_res', 'Topometry - spectral latent space'),
    ('topo_ms', 'topo_clusters_ms_res', 'Topometry - multiscale latent space'),
]


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--adata', required=True, help='topometry AnnData (adata_topometry_*.h5ad)')
    p.add_argument('--color', default='ctype', help='obs column shown next to the clusterings (skipped if missing)')
    p.add_argument('--projection', default='projection', help='obsm key of the embedding to plot on')
    return p


def main():
    args = build_parser().parse_args()

    # only obs, obsm and uns are needed, so leave the expression matrix on disk
    adata = ad.read_h5ad(args.adata, backed='r')
    print(f'latent space chosen by topometry: {adata.uns["basis"]}')

    colours = []
    if args.color in adata.obs.columns:
        colours.append(args.color)
    else:
        print(f'WARNING: {args.color} is not in adata.obs; plotting the clusterings only')

    found, table_cols = [], []
    for name, prefix, title in SERIES:
        # the series in order of resolution (topo_clusters_res0.2, ..., topo_clusters_res1.2)
        cols = sorted((c for c in adata.obs.columns if c.startswith(prefix)), key=lambda c: float(c[len(prefix):]))
        if not cols:
            print(f'{name}: no {prefix}* columns, skipped')
            continue
        print(f'{name}: {len(cols)} resolutions ({prefix}*)')

        adata.obs[cols] = adata.obs[cols].astype(str).astype('category')
        sc.pl.embedding(adata, basis=args.projection, color=colours + cols, legend_loc='on data',
                        legend_fontsize='xx-small', frameon=False, show=False)
        plt.savefig(f'cluster_comparison_on_embedding_{name}.png')
        plt.close()

        found.append((name, prefix, title))
        table_cols += cols

    if not found:
        raise SystemExit(f'None of the clustering series {[s[1] + "*" for s in SERIES]} are in {args.adata}')

    adata.obs[table_cols].to_csv('cluster_table.csv')
    with open('clustree_series.txt', 'w') as f:
        f.writelines(f'{name}\t{prefix}\t{title}\n' for name, prefix, title in found)


if __name__ == '__main__':
    main()
