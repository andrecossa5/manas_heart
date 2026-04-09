process COMPUTE_3NT_CTX {

    publishDir "${params.output_folder}/muts", mode: 'copy'

    input:
    path(tsv)
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    path("ALL_FILTERED_ctx.tsv.gz")

    stub:
    """
    touch ALL_FILTERED_ctx.tsv.gz
    """

    script:
    """
    add_context.py \\
        --input     ${tsv} \\
        --ref       ${ref} \\
        --output    ALL_FILTERED_ctx.tsv.gz \\
        --n-workers ${task.cpus}
    """

}
