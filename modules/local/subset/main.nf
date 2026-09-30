process SUBSET {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/anndatas" }, pattern: 'adata_*.h5ad', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("adata_${meta.tag}.h5ad"), emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    """
    subset.py \\
        --counts ${counts} \\
        --metadata ${meta.metadata} \\
        --subset ${meta.subset} \\
        --output-dir .

    # fail here rather than passing an empty AnnData to topometry
    python - <<PY
    import anndata as ad, sys
    a = ad.read_h5ad('adata_${meta.tag}.h5ad', backed='r')
    print('subset ${meta.tag}:', a.shape)
    if a.n_obs == 0:
        sys.exit("Subset ${meta.metadata} == ${meta.subset} has 0 cells in ${counts}")
    PY
    """
}
