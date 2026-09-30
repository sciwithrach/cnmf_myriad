#!/usr/bin/env nextflow
// cNMF pipeline: subset -> topometry -> prep -> factorize -> combine (+ k selection plot), one chain per samplesheet row.
//   nextflow run main.nf -profile ucl_myriad --samplesheet samplesheet.csv --n_iters 200
// Consensus is a separate step, run after inspecting the k selection plot; then --step analysis for the GEP analysis:
//   nextflow run main.nf -profile ucl_myriad --step consensus --outdir results_20260929 --sample ctype_MEL --selected_k 40 --threshold 0.5

nextflow.enable.dsl = 2

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

include { SUBSET    } from './modules/local/subset/main'
include { TOPOMETRY } from './modules/local/topometry/main'
include { PREP      } from './modules/local/prep/main'
include { FACTORIZE } from './modules/local/factorize/main'
include { COMBINE   } from './modules/local/combine/main'
include { CONSENSUS } from './modules/local/consensus/main'
include { ANALYSIS  } from './modules/local/analysis/main'
include { CLUSTER_COMPARE } from './modules/local/cluster_compare/main'
include { CLUSTREE  } from './modules/local/clustree/main'

// nf-schema has already validated the samplesheet (assets/schema_input.json); this turns a row into a meta map.
// samplesheetToList returns each row as a list in schema property order; empty cells may be null or [].
def makeMeta(row) {
    def (counts, metadata, subset, k_str, n_iters_row, hvgs, run_name) = row.collect { (it instanceof Collection && it.isEmpty()) || it == '' ? null : it }

    def counts_f = file(counts)
    def n_iters  = (n_iters_row ?: params.n_iters) as int      // from the row, else --n_iters (CLI values are strings)
    def tag      = metadata ? "${metadata}_${subset}" : counts_f.baseName
    // k values from the row (space or ; separated), else the --k_vals default
    def k_vals   = (k_str ?: params.k_vals).toString().trim().split(/[\s;]+/).collect { it as int }

    return [
        tag      : tag,
        counts   : counts_f,
        metadata : metadata,
        subset   : subset,
        hvgs     : hvgs ? file(hvgs) : null,
        k_vals   : k_vals,
        run_name : run_name ?: "cnmf_${tag}_${params.run_date}",
        n_iters  : n_iters,
        n_jobs   : k_vals.size() * n_iters,
        n_chunks : Math.ceil(k_vals.size() * n_iters / (params.stride as int)) as int
    ]
}

// Nextflow cannot append to a trace, so each run writes its own trace_<stamp>.txt and its rows are appended
// here to the run-wide pipeline_info/trace.txt and, matched on the tag column, to each sample's own trace.txt
def appendRows(f, header, rows) {
    if( !rows ) return
    f.parent.mkdirs()
    if( !f.exists() ) f.append(header + '\n')
    f.append(rows.join('\n') + '\n')
}

def appendTraces() {
    def raw = file("${params.outdir}/pipeline_info/trace_${params.trace_stamp}.txt")
    if( !raw.exists() ) return
    def lines  = raw.readLines()
    def tagCol = lines[0].split('\t').toList().indexOf('tag')
    if( lines.size() < 2 || tagCol < 0 ) return
    def rows = lines.tail()
    appendRows(file("${params.outdir}/pipeline_info/trace.txt"), lines[0], rows)
    def tags = params.step == 'consensus'
        ? [params.sample]
        : samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json").collect { row -> makeMeta(row).tag }
    tags.each { tag ->
        appendRows(file("${params.outdir}/${tag}/pipeline_info/trace.txt"), lines[0],
                   rows.findAll { it.split('\t', -1)[tagCol] == tag })
    }
}

workflow pipeline {
    if( !params.samplesheet )
        error "Give a samplesheet with --samplesheet"

    // rows are checked against assets/schema_input.json before any job is submitted
    metas = Channel.fromList( samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json") )
        .map { row -> makeMeta(row) }
        .toList()
        .map { list ->
            def dup = list.countBy { it.run_name }.findAll { it.value > 1 }.keySet()
            if( dup ) error "Samplesheet rows give the same run name (set run_name to tell them apart): ${dup}"
            list
        }
        .flatMap { it }

    rows = metas.branch { m ->
        subsetting : m.metadata
        direct     : true
    }

    // rows with metadata + subset: subset -> topometry. These two only get the fields they use, not the whole
    // meta: it carries n_jobs / k_vals / run_name, and Nextflow hashes val inputs, so changing --n_iters or
    // --stride would otherwise invalidate their cache
    SUBSET( rows.subsetting.map { m -> tuple([tag: m.tag, metadata: m.metadata, subset: m.subset], m.counts) } )
    TOPOMETRY( SUBSET.out.res )

    // the full meta is joined back on the tag for prep
    topo_out = TOPOMETRY.out.res
        .map { sm, counts, hvgs -> tuple(sm.tag, counts, hvgs) }
        .join( metas.map { m -> tuple(m.tag, m) } )
        .map { tag, counts, hvgs, m -> tuple(m, counts, hvgs) }

    // rows without: prep straight from the given counts and HVGs
    prep_in = topo_out
        .mix( rows.direct.map { m -> tuple(m, m.counts, m.hvgs) } )
    PREP( prep_in )

    // chunks of params.stride worker indices per task
    FACTORIZE( PREP.out.res.flatMap { m ->
        (1..m.n_jobs).collate(params.stride as int).collect { chunk -> tuple(m, chunk) }
    } )

    // release each row's combine as soon as ITS chunks are done
    COMBINE( FACTORIZE.out.res
        .map { m -> tuple(groupKey(m.tag, m.n_chunks), m) }
        .groupTuple()
        .map { key, ms -> ms[0] } )

    // cluster comparison plots and clustree: they only need the topometry AnnData, so they start as soon as it exists
    // and are listed last so that a resume of an earlier run only adds them
    CLUSTER_COMPARE( TOPOMETRY.out.adata )
    CLUSTREE( CLUSTER_COMPARE.out.table )
}

// the run name defaults to the one PREP recorded for this sample
def runNameFor(sample) {
    def name_file = file("${params.outdir}/${sample}/pipeline_info/run_name.txt")
    def run_name  = params.run_name ?: (name_file.exists() ? name_file.text.trim() : null)
    if( !run_name )
        error "No run name: give --run_name, or check ${name_file} exists (is --outdir the run's results folder?)"
    return run_name
}

workflow consensus {
    if( !params.sample || !params.selected_k )
        error "Consensus needs --sample and --selected_k (and optionally --run_name, --threshold, --outdir)"
    CONSENSUS( Channel.of(runNameFor(params.sample)) )
}

// after consensus: GEP analysis for one sample (needs the same --selected_k / --threshold as the consensus step)
workflow analysis {
    if( !params.sample || !params.selected_k )
        error "The analysis step needs --sample and --selected_k (and the --threshold used for consensus)"
    def run_name = runNameFor(params.sample)
    def sample_dir = "${params.outdir}/${params.sample}"
    ANALYSIS(
        Channel.of( tuple(params.sample, run_name, file("${sample_dir}/anndatas/adata_topometry_${params.sample}.h5ad"),
                          "${file(sample_dir).toAbsolutePath()}/${run_name}") ),
        file(params.gmt),
        file(params.collectri)
    )
}

workflow {
    // Nextflow ignores a stray word without complaint, so a typo such as `threshold 0.18` (no --) would run with the
    // default value. The same happens to `-resume <run name>`, which only accepts `last` or a full session id.
    if( args )
        error "Unexpected argument(s): ${args.join(' ')}\n" +
              "Pipeline options need two dashes (--threshold 0.18). Nextflow options take one, and -resume only accepts a full session id (-resume 1f50da21-...), not a run name."

    validateParameters()
    if( params.step == 'consensus' )
        consensus()
    else if( params.step == 'analysis' )
        analysis()
    else if( params.step == 'pipeline' )
        pipeline()
    else
        error "Unknown --step '${params.step}' (use pipeline, consensus or analysis)"

    workflow.onComplete = { appendTraces() }
}
