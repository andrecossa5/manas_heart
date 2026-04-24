include { PREP_SITES } from './modules/prep_sites'

workflow prep_sites {

    take:
    sites_tsv_ch    // path to sites TSV (CHROM POS REF ALT)
    ref_ch          // tuple(ref, ref_fai, ref_dict)

    main:
    PREP_SITES(sites_tsv_ch, ref_ch)

    emit:
    sites = PREP_SITES.out   // tuple(sites.sorted.tsv, sites.bed, alleles.tsv.gz, alleles.tsv.gz.tbi)

}
