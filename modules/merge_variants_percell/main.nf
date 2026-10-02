// MERGE_VARIANTS_PERCELL: each per-cell pileup parquet merged together with its
// genotype sibling:
//
//   raw_variants.parquet       + variants_count.parquet       -> variants_snv.parquet
//   raw_variants_indel.parquet + variants_indel_count.parquet -> variants_indel.parquet
//
// The genotype parquets are optional inputs (`[]` when a sample produced none),
// because count_variant*.py stops before writing them when there is no genotype
// to report while always writing the pileup table. The merged file then carries
// the pileup rows with null genotype columns.
process MERGE_VARIANTS_PERCELL {
    label 'cpu_low_mem'
    tag "${meta.sample_name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    tuple val(meta), path(snv_raw), path(snv_genotypes), path(indel_raw), path(indel_genotypes)

    output:
    tuple val(meta), path("variants_snv.parquet"), emit: snv_parquet
    tuple val(meta), path("variants_indel.parquet"), emit: indel_parquet
    tuple val(meta), path("status.parquet"), topic: status, optional: true

    script:
    def snv_gt = snv_genotypes ? "--snv-genotypes ${snv_genotypes}" : ''
    def indel_gt = indel_genotypes ? "--indel-genotypes ${indel_genotypes}" : ''
    """
    merge_variants_percell.py ${meta.sample_name} ${snv_raw} ${indel_raw} ${snv_gt} ${indel_gt}
    """

    stub:
    """
    touch variants_snv.parquet
    touch variants_indel.parquet
    """
}
