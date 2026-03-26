include { INDEX_BAM    } from './modules/index_bam'
include { REHEADER_BAM } from './modules/reheader_bam'

workflow prep_bams {

    take:
    chunks_ch   // [tissue, chunk_name, [sample_list]]

    main:

    // Flatten each chunk into one element per sample
    samples_ch = chunks_ch
        .flatMap { tissue, chunk_name, samples ->
            samples.collect { sample -> tuple(tissue, chunk_name, sample) }
        }

    // Resolve BAM path by listing the directory (handles symlinks on Lustre)
    bams_ch = samples_ch
        .map { tissue, chunk_name, sample ->
            def bam = file("${params.bam_folder}/${tissue}")
                .listFiles()
                .findAll { it.name.startsWith(sample) && it.name.endsWith('.bam') }
                .first()
            tuple(tissue, chunk_name, bam)
        }

    // Index each BAM
    INDEX_BAM(bams_ch)

    // Reheader each BAM: set SM tag to tissue name
    REHEADER_BAM(INDEX_BAM.out)

    // Group reheadered BAMs back by tissue + chunk
    grouped_ch = REHEADER_BAM.out.groupTuple(by: [0, 1])
    // → [tissue, chunk_name, [rh_bam1, rh_bam2, ...], [rh_bai1, rh_bai2, ...]]

    emit:
    merged_bams = grouped_ch

}
