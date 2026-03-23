process CREATE_INTERVALS {

    input:
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    path("genomic_intervals/"), emit: intervals_dir

    stub:
    """
    mkdir -p genomic_intervals/
    touch genomic_intervals/0000-scattered.interval_list
    """

    script:
    """
    if grep -qP '\\tSN:chr1\\t' ${ref_dict}; then
        prefix="chr"
    else
        prefix=""
    fi

    if [ "\${prefix}" = "chr" ]; then mt="chrM"; else mt="MT"; fi
    chrom_args=\$(for c in \$(seq 1 22) X Y; do printf -- "-L \${prefix}\${c} "; done)
    chrom_args="\${chrom_args} -L \${mt}"

    mkdir -p genomic_intervals/
    gatk SplitIntervals \\
        -R ${ref} \\
        --scatter-count ${params.scatter_count} \\
        \$chrom_args \\
        --dont-mix-contigs \\
        -O genomic_intervals/
    """
}
