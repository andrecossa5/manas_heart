include { CREATE_INTERVALS } from './modules/create_intervals'
include { MUTECT2          } from './modules/mutect2'
include { GATHER_VCFS      } from './modules/gather_vcfs'
include { FILTER_MUTS      } from './modules/filter_muts'
include { GATHER_TABLE     } from './modules/gather_table'
include { GATHER_STATS     } from './modules/gather_stats'

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

    // Gather: group all interval VCFs back per heart_chunk, then merge
    gathered_ch = MUTECT2.out.vcf
        .map { heart_chunk, interval_name, vcf -> tuple(heart_chunk, vcf) }
        .groupTuple(by: 0)

    GATHER_VCFS(gathered_ch)

    // Filter each gathered VCF
    FILTER_MUTS(GATHER_VCFS.out.unfiltered_vcf)

    // Collect all TSVs and stats across chunks, then gather into single files
    all_tsvs  = FILTER_MUTS.out.tsv.map   { chunk, tsv   -> tsv   }.collect()
    all_stats = FILTER_MUTS.out.stats.map  { chunk, stats -> stats }.collect()

    GATHER_TABLE(all_tsvs)
    GATHER_STATS(all_stats)

    emit:
    unfiltered_vcf  = GATHER_VCFS.out.unfiltered_vcf
    filtered_tsv    = GATHER_TABLE.out
    filtered_stats  = GATHER_STATS.out

}
