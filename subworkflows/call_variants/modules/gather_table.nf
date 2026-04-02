process GATHER_TABLE {

    publishDir "${params.output_folder}/muts", mode: 'copy'

    input:
    path(tsvs)

    output:
    path("ALL_FILTERED.tsv.gz")

    stub:
    """
    touch ALL_FILTERED.tsv.gz
    """

    script:
    """
    gather_tables.py --mode muts --output ALL_FILTERED.tsv.gz *.tsv.gz
    """

}
