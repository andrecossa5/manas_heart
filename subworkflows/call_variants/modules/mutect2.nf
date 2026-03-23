process MUTECT2 {

    tag "${heart_chunk}.${interval_name}"

    input:
    tuple val(heart_chunk),    path(heart_bam),    path(heart_bai),
          val(placenta_chunk), path(placenta_bam), path(placenta_bai),
          val(interval_name),  path(interval_file)
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    tuple val(heart_chunk), val(interval_name),
          path("${heart_chunk}.${interval_name}.unfiltered.vcf.gz"), emit: vcf

    stub:
    """
    touch "${heart_chunk}.${interval_name}.unfiltered.vcf.gz"
    """

    script:
    def germline_arg = params.germline_resource ? "--germline-resource ${params.germline_resource}" : ""

    """
    gatk Mutect2 \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" \\
        -R "${ref}" \\
        -I "${heart_bam}"    --tumor-sample  heart \\
        -I "${placenta_bam}" --normal-sample placenta \\
        -L "${interval_file}" \\
        ${germline_arg} \\
        -O "${heart_chunk}.${interval_name}.unfiltered.vcf.gz"
    """

}
