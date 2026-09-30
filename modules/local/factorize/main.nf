// each task runs one chunk of worker indices (1..n_jobs, params.stride at a time),
// as the old SGE array did. cNMF takes the index modulo total workers, so 1..N covers every run.
process FACTORIZE {
    label 'process_single'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/logs/${task.process}" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.tag}.${task.index}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), val(workers)

    output:
    val meta, emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file(params.outdir).toAbsolutePath()
    """
    for i in ${workers.join(' ')}; do
        cnmf factorize \\
            --output-dir ${outdir} \\
            --name ${meta.run_name} \\
            --worker-index \$i \\
            --total-workers ${meta.n_jobs} \\
            --skip-completed-runs
    done
    """
}
