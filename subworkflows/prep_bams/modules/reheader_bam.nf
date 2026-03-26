process REHEADER_BAM {

    tag "${tissue}.${chunk_name}.${bam.baseName}"

    input:
    tuple val(tissue), val(chunk_name), path(bam), path(bai)

    output:
    tuple val(tissue), val(chunk_name), path("${bam.baseName}.rh.bam"), path("${bam.baseName}.rh.bam.bai")

    stub:
    """
    touch "${bam.baseName}.rh.bam"
    touch "${bam.baseName}.rh.bam.bai"
    """

    script:
    """
    samtools addreplacerg \
        -r "ID:${bam.baseName}\tSM:${tissue}" \
        -o ${bam.baseName}.rh.bam \
        ${bam}
    samtools index ${bam.baseName}.rh.bam
    """

}
