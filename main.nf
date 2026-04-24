#!/usr/bin/env nextflow

include { prep_bams     } from './subworkflows/prep_bams/main'
include { call_variants } from './subworkflows/call_variants/main'
include { prep_sites    } from './subworkflows/prep_sites/main'
include { force_call    } from './subworkflows/force_call/main'


// ── Default workflow: somatic calling on chunk-merged heart vs placenta ─────
workflow {

    // ── Placenta: all samples grouped by chunk ───────────────────────────────
    placenta_ch = Channel.fromPath(params.placenta_csv)
        .splitCsv(header: true)
        .map     { row -> tuple("placenta", row.chunk, row.Sample_ID) }
        .groupTuple(by: [0, 1])
    // → ["placenta", "Placenta", [s1, s2, ...]]

    // ── Heart: group samples by Chunk field ─────────────────────────────────
    heart_ch = Channel.fromPath(params.heart_csv)
        .splitCsv(header: true)
        .map     { row -> tuple("heart", row.chunk, row.Sample_ID) }
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


// ── Named workflow: force-call a known mutation list across individual BAMs ─
workflow forcecall {

    if ( !params.sites ) {
        error "forcecall entrypoint requires --sites <TSV with columns CHROM POS REF ALT>"
    }

    // Per-sample channels (no groupTuple — one element per individual BAM)
    heart_ch = Channel.fromPath(params.heart_csv)
        .splitCsv(header: true)
        .map { row -> tuple("heart", row.chunk, row.Sample_ID) }

    placenta_ch = Channel.fromPath(params.placenta_csv)
        .splitCsv(header: true)
        .map { row -> tuple("placenta", row.chunk, row.Sample_ID) }

    samples_ch = heart_ch.mix(placenta_ch)
        .map { tissue, chunk, sample_id ->
            def bam = file("${params.bam_folder}/${tissue}")
                .listFiles()
                .findAll { it.name.startsWith(sample_id) && it.name.endsWith('.bam') }
                .first()
            tuple(tissue, chunk, sample_id, bam)
        }

    // Reference value channel (same shape as call_variants)
    ref_ch = Channel.value(
        tuple(
            file(params.ref),
            file("${params.ref}.fai"),
            file("${params.ref}".replaceAll(/\.fa(sta)?$/, '.dict'))
        )
    )

    // Prepare sorted sites + BED + bgzipped alleles targets
    prep_sites(Channel.fromPath(params.sites), ref_ch)

    // Force-call per BAM × site-shard → gather long-format allelic table
    force_call(samples_ch, prep_sites.out.sites, ref_ch)

}
