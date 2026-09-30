process COMBINE {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${meta.tag}"
    publishDir { "${params.outdir}/${task.tag}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    val meta

    output:
    val meta, emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file("${params.outdir}/${meta.tag}").toAbsolutePath()   // cnmf adds /${meta.run_name}
    """
    cnmf combine --output-dir ${outdir} --name ${meta.run_name}
    cnmf k_selection_plot --output-dir ${outdir} --name ${meta.run_name}
    """
}
