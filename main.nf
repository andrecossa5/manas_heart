#!/usr/bin/env nextflow

include { prep_bams     } from './subworkflows/prep_bams/main'
include { call_variants } from './subworkflows/call_variants/main'

workflow {

    // ── Placenta: Trophoblasts only, all samples merged into one chunk ──────
    placenta_ch = Channel.fromPath(params.placenta_csv)
        .splitCsv(header: true)
        .filter  { row -> row.Histo == 'Trophoblasts' }
        .map     { row -> tuple("placenta", row.Histo, row.Sample) }
        .groupTuple(by: [0, 1])
    // → ["placenta", "Trophoblasts", [s1, s2, ...]]

    // ── Heart: group samples by Chunk field ─────────────────────────────────
    heart_ch = Channel.fromPath(params.heart_csv)
        .splitCsv(header: true)
        .map     { row -> tuple("heart", row.Chunk, row.Sample) }
        .groupTuple(by: [0, 1])
    // → ["heart", "Myocardium.0", [s1, s2, s3, s4]], ["heart", "Myocardium.1", [...]], ...

    // ── Merge BAMs for all chunks ────────────────────────────────────────────
    chunks_ch = placenta_ch.mix(heart_ch)
    prep_bams(chunks_ch)

    // ── Split merged output by tissue ────────────────────────────────────────
    merged_placenta = prep_bams.out.merged_bams.filter { it[0] == "placenta" }
    merged_heart    = prep_bams.out.merged_bams.filter { it[0] == "heart"    }

    // ── Somatic variant calling: each heart chunk vs merged placenta ─────────
    call_variants(merged_heart, merged_placenta)

}
