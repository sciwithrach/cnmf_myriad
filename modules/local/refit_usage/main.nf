// cNMF 1.7.1 writes wrong usages (its final refit pairs the spectra and the data on different genes), so the final
// step is repeated with the genes in matching order. Reads the consensus outputs and writes
// <run_name>.usages.k_<k>.dt_<threshold>.consensus.refit.txt next to them; the cNMF files are not changed.
process REFIT_USAGE {
    label 'process_medium'
    container "${projectDir}/envs/cnmf_env.sif"
    tag "${params.sample}"
    publishDir { "${params.outdir}/${params.sample}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    val run_name

    output:
    val run_name, emit: res
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    def outdir = file("${params.outdir}/${params.sample}").toAbsolutePath()   // the run folder is ${outdir}/${run_name}
    """
    (
        flock -x 9
        refit_usage.py \\
            --output_dir ${outdir} \\
            --name ${run_name} \\
            --k ${params.selected_k} \\
            --threshold ${params.threshold}
    ) 9>"${outdir}/${run_name}/.consensus.lock"
    """
}
