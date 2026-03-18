process INDEX_BAM {

    tag "${tissue}.${chunk_name}.${bam.baseName}"

    input:
    tuple val(tissue), val(chunk_name), path(bam)

    output:
    tuple val(tissue), val(chunk_name), path(bam), path("${bam}.bai")

    stub:
    """
    touch ${bam}.bai
    """

    script:
    """
    samtools index -@ ${task.cpus} ${bam} ${bam}.bai
    """

}
