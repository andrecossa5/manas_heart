include { CREATE_INTERVALS  } from './modules/create_intervals'
include { MUTECT2           } from './modules/mutect2'
include { FILTER_MUTS       } from './modules/filter_muts'
include { GATHER_TABLE      } from './modules/gather_table'
include { GATHER_STATS      } from './modules/gather_stats'
include { COMPUTE_3NT_CTX   } from './modules/compute_3nt_ctx'

workflow call_variants {

    take:
    merged_heart    // [tissue="heart",    chunk_name, [bams], [bais]]
    merged_placenta // [tissue="placenta", chunk_name, [bams], [bais]]

    main:

    // Drop tissue field; keep [chunk_name, [bams], [bais]] for each channel
    heart_ch    = merged_heart.map    { tissue, chunk, bams, bais -> tuple(chunk, bams, bais) }
    placenta_ch = merged_placenta.map { tissue, chunk, bams, bais -> tuple(chunk, bams, bais) }

    // Reference as value channel
    ref_ch = Channel.value(
        tuple(
            file(params.ref),
            file("${params.ref}.fai"),
            file("${params.ref}".replaceAll(/\.fa(sta)?$/, '.dict'))
        )
    )

    // Create genomic intervals
    CREATE_INTERVALS(ref_ch)

    interval_ch = CREATE_INTERVALS.out.intervals_dir
        .flatMap { dir ->
            dir.listFiles()
               .findAll { it.name.endsWith('-scattered.interval_list') }
               .sort     { a, b -> a.name <=> b.name }
               .collect  { f ->
                   def name = f.name.replace('-scattered.interval_list', '')
                   tuple(name, f)
               }
        }

    // Scatter: (heart × placenta) × interval → one Mutect2 job per combination
    scattered_ch = heart_ch
        .combine(placenta_ch)
        .combine(interval_ch)
    // → [heart_chunk, [heart_bams], [heart_bais], placenta_chunk, [placenta_bams], [placenta_bais], interval_name, interval_file]

    MUTECT2(scattered_ch, ref_ch)

    // Filter each per-interval VCF directly — no VCF merge step
    filter_ch = MUTECT2.out.vcf
        .map { heart_chunk, interval_name, vcf -> tuple(heart_chunk, interval_name, vcf) }

    FILTER_MUTS(filter_ch)

    // Collect all TSVs and stats across all intervals and chunks
    all_tsvs  = FILTER_MUTS.out.tsv.map   { chunk, interval, tsv   -> tsv   }.collect()
    all_stats = FILTER_MUTS.out.stats.map  { chunk, interval, stats -> stats }.collect()

    GATHER_TABLE(all_tsvs)
    GATHER_STATS(all_stats)

    COMPUTE_3NT_CTX(GATHER_TABLE.out, ref_ch)

    emit:
    filtered_tsv    = COMPUTE_3NT_CTX.out
    filtered_stats  = GATHER_STATS.out

}
