process CHUNK_SITES {

    tag "chunk_sites"

    input:
    tuple path(sites_sorted_tsv), path(sites_bed), path(alleles_tsv_gz), path(alleles_tbi)
    val   n_shards

    output:
    path("shard_*.tsv")

    stub:
    """
    touch shard_000.tsv
    """

    script:
    """
    # Drop header, split body into N shards with equal line counts
    tail -n +2 ${sites_sorted_tsv} > body.tsv
    total=\$(wc -l < body.tsv)
    n=${n_shards}
    if [ "\$total" -eq 0 ]; then
        echo "ERROR: sites file is empty" >&2; exit 1
    fi
    if [ "\$n" -gt "\$total" ]; then n=\$total; fi
    # Ceiling division so all rows are covered
    per=\$(( (total + n - 1) / n ))
    split -d -a 3 -l "\$per" --additional-suffix=.tsv body.tsv shard_
    """

}
