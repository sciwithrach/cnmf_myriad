#!/usr/bin/env python3
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb
import topo as tp

sc.settings.figdir = '.'


def pca_umap_baseline(adata, npcs=40, savename=None):
    sc.pp.pca(adata, layer="scaled", n_comps=100)
    sc.pl.pca_variance_ratio(adata, n_pcs=100, show=False, save=f'_{savename}.png')
    plt.close()

    sc.pp.neighbors(adata, n_pcs=npcs, use_rep='X_pca', metric='cosine', key_added='PCA')

    for res in [0.2, 0.4, 0.6, 0.8, 1, 1.2]:
        sc.tl.leiden(
            adata,
            resolution=res,
            neighbors_key='PCA',
            flavor='igraph',
            key_added=f'pca_leiden_res{res}',
            random_state=42
        )

    sc.tl.umap(adata, neighbors_key='PCA')
    adata.obsm['X_UMAP_on_PCA'] = adata.obsm['X_umap'].copy()
    del adata.obsm['X_umap']

    sc.pp.neighbors(adata, n_pcs=0, use_rep='X', metric='cosine', key_added='X')
    sc.tl.umap(adata, neighbors_key='X')
    adata.obsm['X_UMAP_on_adataX'] = adata.obsm['X_umap'].copy()
    del adata.obsm['X_umap']


def run_topometry(adata_path, descriptor, npcs=30):
    # read adata
    adata = ad.read_h5ad(adata_path)

    # select raw if available (all genes), else fall back to the counts layer
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    elif 'counts' in adata.layers:
        adata.X = adata.layers['counts']

    # topometry computes its own HVGs, so it needs unnormalised counts -- check
    # the first few nonzero values are whole numbers (all cells/genes go
    # through the same upstream processing, so a sample is enough), and report
    # the gene count for debugging
    sample = adata.X.data[:100] if sp.issparse(adata.X) else np.asarray(adata.X)[np.asarray(adata.X) != 0][:100]
    is_whole_numbers = np.allclose(sample, np.round(sample), atol=1e-6)
    print(f'{adata_path}: n_genes = {adata.n_vars:,}')
    assert is_whole_numbers, f'{adata_path}: adata.X does not look like unnormalised counts (topometry needs raw counts for its own HVG selection)'

    # tp.sc.preprocess would otherwise overwrite .raw with a log-normalised
    # snapshot (taken after normalize_total/log1p, before HVG subsetting) --
    # set the true unnormalised, all-gene raw ourselves beforehand, and tell
    # preprocess not to overwrite it
    adata.raw = adata.copy()

    # run standard pre-processing
    adata = tp.sc.preprocess(adata, save_to_raw=False)

    # generate the PCAxUMAP baseline
    pca_umap_baseline(adata, npcs=npcs, savename=descriptor)

    # set seed and resolutions for clustering
    seed = 42
    leiden_resolutions = (0.2, 0.4, 0.6, 0.8, 1, 1.2)

    # fit topograph object
    tg = tp.sc.fit_adata(
        adata,
        projections=("MAP", "PaCMAP"),
        do_leiden=True,
        leiden_resolutions=leiden_resolutions,
        n_jobs=-1,
        verbosity=0,
        random_state=seed
    )

    tp.sc.intrinsic_dim(
        adata,
        tg=tg,
        n_jobs=-1,
        id_methods=['fsa', 'mle'],
        id_k_values=None
    )

    tp.sc.evaluate_representations(
        adata,
        tg,
        return_df=False,
        print_results=True,
        plot_results=True,
        plot_path=f'geometry_preservation_{descriptor}.png',
        n_neighbors=npcs,
        n_jobs=-1,
        times=(1, 2, 4),
        r=32,
        k_for_pf1=None,
    )

    return adata.copy(), tg


def find_best_embedding(adata, descriptor):
    df = adata.uns['topometry_representation_eval'].copy()
    tag = 'qc' if descriptor == 'postqc' else f'before_{descriptor}'

    # rename any existing basis/projection from a prior run so they're included in comparison
    for rep in ['basis', 'projection']:
        renamed = f'{rep}_{tag}'
        if rep in df['representation'].values:
            df.loc[df['representation'] == rep, 'representation'] = renamed
        if rep in adata.uns:
            adata.uns[renamed] = adata.uns.pop(rep)
        if rep in adata.obsm:
            adata.obsm[f'X_{renamed}'] = adata.obsm[rep]
            del adata.obsm[rep]

    base_list = ['spectral_scaffold', 'ms_spectral_scaffold', 'pca', f'basis_{tag}']
    bases = df[df.representation.isin(base_list)]
    projections = df[~df.representation.isin(base_list)].copy()
    # projection_{tag} naturally falls into projections

    best_base = bases.sort_values(['SP', 'PF1', 'PJS'], ascending=False)['representation'].iloc[0]
    print(f'Selected basis : {best_base}')

    projections['mean_score'] = projections[['PF1', 'PJS', 'SP']].mean(axis=1)
    best_projection = projections.sort_values('mean_score', ascending=False)['representation'].iloc[0]
    print(f'Selected projection : {best_projection}')

    adata.uns["basis"] = best_base
    adata.obsm["basis"] = adata.obsm[f'X_{best_base}']
    adata.uns["projection"] = best_projection
    adata.obsm["projection"] = adata.obsm[f'X_{best_projection}']

    return adata.copy()


def plot_comparisons(adata, descriptor, metadata=None):
    representations = [
        'pca', 'UMAP_on_PCA', 'TopoMAP', 'TopoPaCMAP', 'msTopoMAP', 'msTopoPaCMAP'
    ]

    fig, axs = plt.subplots(3, 2, figsize=(12, 15))
    axes = [axs[0,0], axs[0,1], axs[1,0], axs[1,1], axs[2,0], axs[2,1]]

    for i, rep in enumerate(representations):
        sc.pl.embedding(adata, basis=rep, ax=axes[i], show=False, title=rep,
                        color=metadata if metadata else None, frameon=False)
    plt.savefig(f'embeddings_comparison_{descriptor}.png')
    plt.close()


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--descriptor', required=True)
    p.add_argument('--metadata', default=None)
    return p


def main():

    # parse arguments
    args = build_parser().parse_args()
    
    # run topo
    adata, tg = run_topometry(args.adata, args.descriptor)

    # find best embedding
    adata = find_best_embedding(adata, args.descriptor)

    # plot comparisons
    plot_comparisons(adata, args.descriptor, args.metadata)

    # save adata
    adata.write(f'adata_topometry_{args.descriptor}.h5ad')

    # save topograph object
    tp.save_topograph(tg, f'topometry_object_{args.descriptor}.pkl')

    # save list of highly variable genes
    hvgs = list(adata.var_names)    
    pd.Dataframe(hvgs).to_csv('hvgs_{args.descriptor}.csv', header=False, index=False)

if __name__ == '__main__':
    main()
