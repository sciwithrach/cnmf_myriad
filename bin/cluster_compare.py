#!/usr/bin/env python3
"""Cluster comparison plot and clustree input for one topometry AnnData.

Picks the clustering series that matches the latent space topometry chose (`adata.uns['basis']`):
    ms_spectral*  -> topo_clusters_ms_res*  (Topometry - multiscale latent space)
    spectral*     -> topo_clusters_res*     (Topometry - spectral latent space)
    otherwise     -> pca_leiden_res*        (PCA-based latent space)

Writes, in the current directory:
    cluster_comparison_on_embedding.png   the projection coloured by --color and by every resolution in the series
    cluster_table.csv                     cells x resolutions, the input for bin/clustree.R
    clustree_series.txt                   the column prefix and plot title for bin/clustree.R
"""
import argparse

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc


def choose_series(basis: str):
    """Column prefix and clustree title for the latent space named `basis`."""
    if 'ms_spectral' in basis:
        return 'topo_clusters_ms_res', 'Topometry - multiscale latent space'
    if 'spectral' in basis:
        return 'topo_clusters_res', 'Topometry - spectral latent space'
    return 'pca_leiden_res', 'PCA-based latent space'


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
    basis = str(adata.uns['basis'])
    prefix, title = choose_series(basis)

    # the series, in order of resolution (topo_clusters_res0.2, ..., topo_clusters_res1.2)
    cluster_cols = sorted((c for c in adata.obs.columns if c.startswith(prefix)), key=lambda c: float(c[len(prefix):]))
    if not cluster_cols:
        raise SystemExit(f'No columns starting with {prefix!r} in {args.adata} (basis: {basis})')
    print(f'basis: {basis} -> {prefix}* ({len(cluster_cols)} resolutions)')

    colours = []
    if args.color in adata.obs.columns:
        colours.append(args.color)
    else:
        print(f'WARNING: {args.color} is not in adata.obs; plotting the clusterings only')

    adata.obs[cluster_cols] = adata.obs[cluster_cols].astype(str).astype('category')
    sc.pl.embedding(
        adata,
        basis=args.projection,
        color=colours + cluster_cols,
        legend_loc='on data',
        legend_fontsize='xx-small',
        frameon=False,
        show=False,
    )
    plt.savefig('cluster_comparison_on_embedding.png')
    plt.close()

    adata.obs[cluster_cols].to_csv('cluster_table.csv')
    with open('clustree_series.txt', 'w') as f:
        f.write(f'{prefix}\n{title}\n')


if __name__ == '__main__':
    main()
