process PREP_SITES {

    tag "prep_sites"

    publishDir "${params.output_folder}/forcecall", mode: 'copy'

    input:
    path(sites_tsv)
    tuple path(ref), path(ref_fai), path(ref_dict)

    output:
    tuple path("sites.sorted.tsv"),
          path("sites.bed"),
          path("alleles.tsv.gz"),
          path("alleles.tsv.gz.tbi")

    stub:
    """
    touch sites.sorted.tsv sites.bed alleles.tsv.gz alleles.tsv.gz.tbi
    """

    script:
    """
    # Validate header
    read -r HDR < <(zcat -f ${sites_tsv} | head -n1)
    echo "[prep_sites] header: \$HDR" >&2
    case "\$HDR" in
        CHROM*POS*REF*ALT*) ;;
        *) echo "ERROR: sites TSV must start with header 'CHROM<TAB>POS<TAB>REF<TAB>ALT'" >&2; exit 1 ;;
    esac

    # Build a contig order from the .fai (stable with reference)
    awk 'BEGIN{OFS="\t"} {print \$1, NR}' ${ref_fai} > contig_order.tsv

    # Sort sites by reference contig order, then POS (numeric). Keep only CHROM,POS,REF,ALT.
    zcat -f ${sites_tsv} \\
      | awk 'BEGIN{FS=OFS="\t"} NR>1 {print \$1,\$2,\$3,\$4}' \\
      | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ord[\$1]=\$2; next} (\$1 in ord){print ord[\$1], \$0}' contig_order.tsv - \\
      | sort -k1,1n -k3,3n \\
      | cut -f2- \\
      > sites.sorted.body.tsv

    { printf "CHROM\tPOS\tREF\tALT\n"; cat sites.sorted.body.tsv; } > sites.sorted.tsv

    # BED (0-based, half-open) for bcftools mpileup -R
    awk 'BEGIN{FS=OFS="\t"} {print \$1, \$2-1, \$2}' sites.sorted.body.tsv > sites.bed

    # Alleles file for bcftools call -C alleles -T : CHROM \\t POS \\t REF,ALT
    awk 'BEGIN{FS=OFS="\t"} {print \$1, \$2, \$3","\$4}' sites.sorted.body.tsv \\
      | bgzip -c > alleles.tsv.gz
    tabix -s1 -b2 -e2 -f alleles.tsv.gz
    """

}
