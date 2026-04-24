process PARSE_PILEUP {

    tag "${tissue}.${sample_id}.${shard_name}"

    input:
    tuple val(tissue), val(chunk), val(sample_id), val(shard_name),
          path(vcf), path(shard_tsv)

    output:
    path("${sample_id}.${shard_name}.allelic.tsv.gz")

    stub:
    """
    touch "${sample_id}.${shard_name}.allelic.tsv.gz"
    """

    script:
    """
    parse_forcecall_vcf.py \\
        --vcf        ${vcf} \\
        --sites      ${shard_tsv} \\
        --sample-id  ${sample_id} \\
        --tissue     ${tissue} \\
        --chunk      ${chunk} \\
        --output     "${sample_id}.${shard_name}.allelic.tsv.gz"
    """

}
