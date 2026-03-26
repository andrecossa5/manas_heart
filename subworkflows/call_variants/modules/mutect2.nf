process MUTECT2 {

    tag "${heart_chunk}.${interval_name}"

    input:
    tuple val(heart_chunk),    path(heart_bams),    path(heart_bais),
          val(placenta_chunk), path(placenta_bams), path(placenta_bais),
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
    def heart_i    = (heart_bams    instanceof List ? heart_bams    : [heart_bams])
                     .collect { "-I ${it} --tumor-sample  heart"    }.join(" \\\n        ")
    def placenta_i = (placenta_bams instanceof List ? placenta_bams : [placenta_bams])
                     .collect { "-I ${it} --normal-sample placenta" }.join(" \\\n        ")

    """
    gatk Mutect2 \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" \\
        -R "${ref}" \\
        ${heart_i} \\
        ${placenta_i} \\
        -L "${interval_file}" \\
        ${germline_arg} \\
        -O "${heart_chunk}.${interval_name}.unfiltered.vcf.gz"
    """

}
