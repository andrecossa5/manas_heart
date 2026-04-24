include { INDEX_BAM             } from '../prep_bams/modules/index_bam'
include { CHUNK_SITES           } from './modules/chunk_sites'
include { MPILEUP_FORCECALL     } from './modules/mpileup_forcecall'
include { PARSE_PILEUP          } from './modules/parse_pileup'
include { GATHER_ALLELIC_TABLE  } from './modules/gather_allelic_table'

workflow force_call {

    take:
    samples_ch   // [tissue, chunk, sample_id, bam_path]
    sites_ch     // tuple(sites.sorted.tsv, sites.bed, alleles.tsv.gz, alleles.tsv.gz.tbi)
    ref_ch       // tuple(ref, ref_fai, ref_dict)

    main:

    // Index each BAM once — reuse INDEX_BAM by feeding the (tissue, chunk, bam) it expects
    index_in_ch = samples_ch.map { tissue, chunk, sample_id, bam ->
        tuple(tissue, "${chunk}::${sample_id}", bam)
    }
    INDEX_BAM(index_in_ch)

    // Re-attach sample_id; INDEX_BAM emits (tissue, chunk_name, bam, bai) where chunk_name == "${chunk}::${sample_id}"
    indexed_ch = INDEX_BAM.out.map { tissue, tagged, bam, bai ->
        def parts     = tagged.toString().split('::', 2)
        def chunk     = parts[0]
        def sample_id = parts[1]
        tuple(tissue, chunk, sample_id, bam, bai)
    }

    // Shard the sites
    shards_raw_ch = CHUNK_SITES(sites_ch, params.sites_scatter_count)
        .flatten()
        .map { f -> tuple(f.name.replace('.tsv', ''), f) }
    // → [shard_name, shard_tsv]

    // Cartesian product: every BAM × every shard
    scattered_ch = indexed_ch.combine(shards_raw_ch)
    // → [tissue, chunk, sample_id, bam, bai, shard_name, shard_tsv]

    MPILEUP_FORCECALL(scattered_ch, ref_ch)

    PARSE_PILEUP(MPILEUP_FORCECALL.out.vcf)

    all_allelic = PARSE_PILEUP.out.collect()

    GATHER_ALLELIC_TABLE(all_allelic)

    emit:
    allelic_table = GATHER_ALLELIC_TABLE.out

}
