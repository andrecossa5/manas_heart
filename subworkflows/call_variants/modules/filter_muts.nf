process FILTER_MUTS {

    tag "${heart_chunk}"

    input:
    tuple val(heart_chunk), path(vcf)

    output:
    tuple val(heart_chunk), path("${heart_chunk}_filtered.tsv"),    emit: tsv
    tuple val(heart_chunk), path("${heart_chunk}_filtered.stats"),   emit: stats

    stub:
    """
    touch "${heart_chunk}_filtered.tsv"
    touch "${heart_chunk}_filtered.stats"
    """

    script:
    """
    filter_mutect2.py \\
        ${vcf} \\
        --chunk             ${heart_chunk} \\
        --ad-placenta-max   ${params.ad_placenta_max} \\
        --ad-heart-min      ${params.ad_heart_min} \\
        --dp-min            ${params.dp_min} \\
        --af-heart-max      ${params.af_heart_max} \\
        --mmq-min           ${params.mmq_min} \\
        --mbq-min           ${params.mbq_min} \\
        --sb-pval-min       ${params.sb_pval_min} \\
        --mpos-min          ${params.mpos_min} \\
        --alt2-ad-ratio-max ${params.alt2_ad_ratio_max} \\
        -o ${heart_chunk}_filtered.tsv
    """

}
