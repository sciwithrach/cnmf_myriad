// Cluster comparison plot, and the table of clusterings that CLUSTREE plots. Needs only the topometry AnnData.
process CLUSTER_COMPARE {
    container "${projectDir}/envs/cnmf-analysis.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/clustering" }, pattern: 'cluster_comparison_on_embedding.png', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path('cluster_table.csv'), path('clustree_series.txt'), emit: table
    path 'cluster_comparison_on_embedding.png'
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    """
    # matplotlib / numba want writable cache folders; the container's home may not be
    export MPLCONFIGDIR=\$PWD/.matplotlib NUMBA_CACHE_DIR=\$PWD/.numba

    cluster_compare.py \\
        --adata ${adata} \\
        --color ${params.cluster_color} \\
        --projection ${params.projection}
    """
}
