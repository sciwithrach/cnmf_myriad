#!/usr/bin/env nextflow

/*
 * dimensionality reduction
 */

 process run_topometry {

    container "$projectDir/containers/utricle-qc.sif"

    input :
    path adata
    val metadata
    val descriptor

    output:
    path "adata_topometry_${descriptor}.h5ad", emit : adata
    path "topometry_object_${descriptor}.pkl", emit : tg
    path "pca_variance_ratio_${descriptor}.png", emit : plt_pca
    path "geometry_preservation_${descriptor}.png", emit : plt_geometry
    path "embeddings_comparison_${descriptor}.png", emit : plt_comparison

    script:
    """
    topometry.py --adata ${adata} --descriptor ${descriptor} --metadata ${metadata}
    """
 }
