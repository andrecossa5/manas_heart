include { INDEX_BAM } from './modules/index_bam'
include { MERGE_BAM } from './modules/merge_bam'

workflow prep_bams {

    take:
    chunks_ch   // [tissue, chunk_name, [sample_list]]

    main:

    // Flatten each chunk into one element per sample
    samples_ch = chunks_ch
        .flatMap { tissue, chunk_name, samples ->
            samples.collect { sample -> tuple(tissue, chunk_name, sample) }
        }

    // Resolve BAM path from bam_folder glob
    bams_ch = samples_ch
        .map { tissue, chunk_name, sample ->
            def bam = file("${params.bam_folder}/${tissue}/${sample}*.bam")
                .findAll { it.name.endsWith('.bam') }
                .first()
            tuple(tissue, chunk_name, bam)
        }

    // Index each BAM
    INDEX_BAM(bams_ch)

    // Group indexed BAMs back by tissue + chunk
    grouped_ch = INDEX_BAM.out.groupTuple(by: [0, 1])
    // → [tissue, chunk_name, [bams], [bais]]

    // Merge grouped BAMs and reheader SM tags
    MERGE_BAM(grouped_ch)

    emit:
    merged_bams = MERGE_BAM.out
    // → [tissue, chunk_name, merged_bam, merged_bai]

}
