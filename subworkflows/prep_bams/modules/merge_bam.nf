process MERGE_BAM {

    tag "${tissue}.${chunk_name}"

    input:
    tuple val(tissue), val(chunk_name), path(bams), path(bais)

    output:
    tuple val(tissue), val(chunk_name), path("${tissue}.${chunk_name}.bam"), path("${tissue}.${chunk_name}.bam.bai")

    stub:
    """
    touch ${tissue}.${chunk_name}.bam
    touch ${tissue}.${chunk_name}.bam.bai
    """

    script:
    """
    samtools merge -@ ${task.cpus} -f merged_raw.bam ${bams}
    samtools reheader -c "sed 's/SM:[^\\t]*/SM:${tissue}/g'" merged_raw.bam > ${tissue}.${chunk_name}.bam
    samtools index ${tissue}.${chunk_name}.bam
    """

}
