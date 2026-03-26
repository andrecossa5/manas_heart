process GATHER_STATS {

    publishDir "${params.output_folder}/muts", mode: 'copy'

    input:
    path(stats)

    output:
    path("ALL_FILTERED.stats")

    stub:
    """
    touch ALL_FILTERED.stats
    """

    script:
    """
    gather_tables.py --mode stats --output ALL_FILTERED.stats *.stats
    """

}
