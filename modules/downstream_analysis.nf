process DOWNSTREAM_ANALYSIS {
    label 'small_job'
    container "community.wave.seqera.io/library/gcc_pip_mudata_polars:966b9cfe2b439554"
    tag "${meta.sample_name}"
    debug false

    input:
    tuple val(meta), path(mudata_h5mu), path(mudata_zarr)

    output:
    tuple val(meta), path("*.png"), emit: qc

    script:
    """
    downstream_analysis.py ${meta.sample_name} ${mudata_zarr}
    """

    stub:
    """
    touch foo.png
    """
}
