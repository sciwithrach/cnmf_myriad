// Cluster stability across resolutions (clustree), from the table written by CLUSTER_COMPARE.
process CLUSTREE {
    container "${projectDir}/envs/clustree.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/clustering" }, pattern: 'clustree.png', mode: 'copy'
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(table), path(series)

    output:
    path 'clustree.png'
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    """
    clustree.R ${table} ${series}
    """
}
