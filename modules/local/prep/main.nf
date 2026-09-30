process PREP {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/logs/${task.process}" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.tag}${fn.replace('.command', '')}" }

    input:
    tuple val(meta), path(counts), path(hvgs)

    output:
    val meta, emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file(params.outdir).toAbsolutePath()
    """
    cnmf prepare \\
        --output-dir ${outdir} \\
        --name ${meta.run_name} \\
        -c ${counts} \\
        -k ${meta.k_vals.join(' ')} \\
        --n-iter ${params.n_iters} \\
        --genes-file ${hvgs} \\
        --seed ${params.seed} \\
        --total-workers ${meta.n_jobs}
    """
}
