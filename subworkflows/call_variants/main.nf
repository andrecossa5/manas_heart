include { MUTECT2 } from './modules/mutect2'

workflow CALL_VARIANTS {

    take:
    merged_heart    // [tissue="heart",    chunk_name, bam, bai]
    merged_placenta // [tissue="placenta", chunk_name, bam, bai]

    main:

    // Drop tissue field; keep [chunk_name, bam, bai] for each channel
    heart_ch    = merged_heart.map    { tissue, chunk, bam, bai -> tuple(chunk, bam, bai) }
    placenta_ch = merged_placenta.map { tissue, chunk, bam, bai -> tuple(chunk, bam, bai) }

    // Cartesian product: each heart chunk paired with the single placenta
    combined_ch = heart_ch.combine(placenta_ch)
    // → [heart_chunk, heart_bam, heart_bai, placenta_chunk, placenta_bam, placenta_bai]

    // Reference as value channel
    ref_ch = Channel.value(
        tuple(
            file(params.ref),
            file("${params.ref}.fai"),
            file("${params.ref}".replaceAll(/\.fa(sta)?$/, '.dict'))
        )
    )

    MUTECT2(combined_ch, ref_ch)

    emit:
    unfiltered_vcf = MUTECT2.out.unfiltered_vcf
    // → [heart_chunk, vcf, tbi, stats]

}
