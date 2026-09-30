// Downstream GEP analysis of one sample after consensus: merges the GEP usages into the AnnData, top genes,
// GOBP + CollecTRI enrichment, and the figures. Writes <outdir>/<sample>/analysis/{csvs,excel,figures,anndatas}.
process ANALYSIS {
    container "${projectDir}/envs/cnmf-analysis.sif"
    tag "${sample}"
    publishDir { "${params.outdir}/${sample}" }, pattern: 'analysis', mode: 'copy'
    publishDir { "${params.outdir}/${sample}/logs" }, pattern: '.command.{out,err}', mode: 'copy',
        saveAs: { fn -> "${task.process}${fn.replace('.command', '')}" }

    input:
    tuple val(sample), val(run_name), path(adata), val(run_dir)
    path gmt
    path collectri

    output:
    path 'analysis'
    path '.command.out', hidden: true, optional: true
    path '.command.err', hidden: true, optional: true

    script:
    """
    # matplotlib / numba want writable cache folders; the container's home may not be
    export MPLCONFIGDIR=\$PWD/.matplotlib NUMBA_CACHE_DIR=\$PWD/.numba

    gep_analysis.py \\
        --adata ${adata} \\
        --run_dir ${run_dir} \\
        --run_name ${run_name} \\
        --k ${params.selected_k} \\
        --threshold ${params.threshold} \\
        --gmt ${gmt} \\
        --collectri ${collectri} \\
        --sample ${sample} \\
        --out analysis \\
        --projection ${params.projection} \\
        --clusters ${params.clusters} \\
        --age ${params.age} \\
        --summary_cols ${params.summary_cols} \\
        --usage_cutoff ${params.usage_cutoff} \\
        --n_top_genes ${params.n_top_genes} \\
        --geps_per_page ${params.geps_per_page} \\
        --formats ${params.formats}
    """
}
