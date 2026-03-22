process MUTECT2 {

    tag "${heart_chunk}"

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    tuple val(heart_chunk),    path(heart_bam),    path(heart_bai),
          val(placenta_chunk), path(placenta_bam), path(placenta_bai)
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    tuple val(heart_chunk), path("${heart_chunk}_unfiltered.vcf.gz"), emit: unfiltered_vcf

    stub:
    """
    touch "${heart_chunk}_unfiltered.vcf.gz"
    """

    script:
    def germline_arg = params.germline_resource ? "--germline-resource ${params.germline_resource}" : ""
    
    """
    gatk Mutect2 \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" \\
        -R "${ref}" \\
        -I "${heart_bam}"    --tumor-sample  heart \\
        -I "${placenta_bam}" --normal-sample placenta \\
        ${germline_arg} \\
        -O "${heart_chunk}_unfiltered.vcf.gz"
    """

}
