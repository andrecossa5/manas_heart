process MPILEUP_FORCECALL {

    tag "${tissue}.${sample_id}.${shard_name}"

    input:
    tuple val(tissue), val(chunk), val(sample_id),
          path(bam), path(bai),
          val(shard_name), path(shard_tsv)
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    tuple val(tissue), val(chunk), val(sample_id), val(shard_name),
          path("${sample_id}.${shard_name}.vcf.gz"),
          path(shard_tsv),                                 emit: vcf

    stub:
    """
    touch "${sample_id}.${shard_name}.vcf.gz"
    """

    script:
    """
    # Build per-shard BED and alleles file
    awk 'BEGIN{FS=OFS="\t"} {print \$1, \$2-1, \$2}' ${shard_tsv} > shard.bed
    awk 'BEGIN{FS=OFS="\t"} {print \$1, \$2, \$3","\$4}' ${shard_tsv} \\
        | bgzip -c > shard_alleles.tsv.gz
    tabix -s1 -b2 -e2 -f shard_alleles.tsv.gz

    bcftools mpileup \\
        -f ${ref} \\
        -R shard.bed \\
        -a FORMAT/AD,FORMAT/DP \\
        -q ${params.min_mq_forcecall} \\
        -Q ${params.min_bq_forcecall} \\
        -Ou \\
        ${bam} \\
      | bcftools call \\
          -Aim \\
          -C alleles \\
          -T shard_alleles.tsv.gz \\
          -Oz \\
          -o "${sample_id}.${shard_name}.vcf.gz"

    tabix -p vcf "${sample_id}.${shard_name}.vcf.gz"
    """

}
