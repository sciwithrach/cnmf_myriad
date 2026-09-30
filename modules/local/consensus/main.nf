// consensus writes a cache in the run folder, so combine + consensus run under an exclusive lock
process CONSENSUS {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${run_name}"
    publishDir { "${params.outdir}/logs/${task.process}" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.tag}${fn.replace('.command', '')}" }

    input:
    val run_name

    output:
    val run_name, emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file(params.outdir).toAbsolutePath()
    """
    (
        flock -x 9
        cnmf combine --output-dir ${outdir} --name ${run_name}
        cnmf consensus \\
            --output-dir ${outdir} \\
            --name ${run_name} \\
            --components ${params.selected_k} \\
            --local-density-threshold ${params.threshold} \\
            --show-clustering
    ) 9>"${outdir}/${run_name}/.consensus.lock"
    """
}
