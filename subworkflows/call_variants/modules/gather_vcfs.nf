process GATHER_VCFS {

    tag "${heart_chunk}"

    publishDir "${params.output_folder}/vcfs", mode: 'copy'

    input:
    tuple val(heart_chunk), path(vcfs)

    output:
    tuple val(heart_chunk), path("${heart_chunk}_unfiltered.vcf.gz"), emit: unfiltered_vcf

    stub:
    """
    touch "${heart_chunk}_unfiltered.vcf.gz"
    """

    script:
    """
    ls *.vcf.gz | sort > vcf.list

    gatk MergeVcfs \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" \\
        -I vcf.list \\
        -O "${heart_chunk}_unfiltered.vcf.gz"
    """
}
