process TOPOMETRY {
    label 'process_high'
    container "${projectDir}/envs/utricle-qc.sif"
    containerOptions '--bind /usr/bin/time:/usr/bin/time'   // GNU time is not in the container
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/anndatas" }, pattern: 'adata_topometry_*.h5ad', mode: 'copy',
        saveAs: { fn -> "adata_topometry_${task.tag}.h5ad" }
    publishDir { "${params.outdir}/${task.tag}" }, pattern: 'hvgs_*.csv', mode: 'copy',
        saveAs: { fn -> "hvgs_${task.tag}.csv" }
    publishDir { "${params.outdir}/${task.tag}/topometry" }, pattern: '{*.png,*.pkl,figures}', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("counts_${meta.subset}.h5ad"), path("hvgs_${meta.subset}.csv"), emit: res
    tuple val(meta), path("adata_topometry_${meta.subset}.h5ad"), emit: adata
    path '*.png', optional: true
    path '*.pkl', optional: true
    path 'figures', optional: true
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    // bin/topometry.py includes the sc.pp.filter_genes(min_cells=3) step
    """
    /usr/bin/time --verbose topometry.py \\
        --adata ${adata} \\
        --descriptor ${meta.subset} \\
        --metadata ${meta.metadata}
    """
}
