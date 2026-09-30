#!/usr/bin/env nextflow
// cNMF pipeline: subset -> topometry -> prep -> factorize -> combine (+ k selection plot), one chain per samplesheet row.
//   nextflow run main.nf -profile ucl_myriad --samplesheet samplesheet.csv --n_iters 200
// Consensus is a separate entry point, run after inspecting the k selection plot:
//   nextflow run main.nf -profile ucl_myriad -entry consensus --run_name cnmf_ctype_MEL_20260929 --selected_k 40 --threshold 0.5

nextflow.enable.dsl = 2

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

include { SUBSET    } from './modules/local/subset/main'
include { TOPOMETRY } from './modules/local/topometry/main'
include { PREP      } from './modules/local/prep/main'
include { FACTORIZE } from './modules/local/factorize/main'
include { COMBINE   } from './modules/local/combine/main'
include { CONSENSUS } from './modules/local/consensus/main'

// nf-schema has already validated the samplesheet (assets/schema_input.json); this turns a row into a meta map.
// samplesheetToList returns each row as a list in schema property order; empty cells may be null or [].
def makeMeta(row) {
    def (counts, metadata, subset, k_str, hvgs, run_name) = row.collect { (it instanceof Collection && it.isEmpty()) || it == '' ? null : it }

    def counts_f = file(counts)
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
        n_jobs   : k_vals.size() * params.n_iters,
        n_chunks : Math.ceil(k_vals.size() * params.n_iters / params.stride) as int
    ]
}

workflow {
    validateParameters()
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

    // rows with metadata + subset: subset -> topometry
    SUBSET( rows.subsetting.map { m -> tuple(m, m.counts) } )
    TOPOMETRY( SUBSET.out.res )

    // rows without: prep straight from the given counts and HVGs
    prep_in = TOPOMETRY.out.res
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
}

workflow consensus {
    validateParameters()
    if( !params.run_name || !params.selected_k )
        error "Consensus needs --run_name and --selected_k (and optionally --threshold)"
    CONSENSUS( Channel.of(params.run_name) )
}
