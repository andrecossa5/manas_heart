process GATHER_ALLELIC_TABLE {

    publishDir "${params.output_folder}/forcecall", mode: 'copy'

    input:
    path(tsvs)

    output:
    path("ALLELIC_TABLE.tsv.gz")

    stub:
    """
    touch ALLELIC_TABLE.tsv.gz
    """

    script:
    """
    gather_tables.py --mode allelic --output ALLELIC_TABLE.tsv.gz *.allelic.tsv.gz
    """

}
