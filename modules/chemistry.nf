#!/usr/bin/env nextflow

process UPDATE_DEFINITIONS {
    label 'small_job'
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    path barcode_def_path
    path pipeline_root_dir
    path chemistry_schema

    output:
    path ("chemistry_defs.json"), emit: updated_config

    script:

    """
        update_defs.py ${barcode_def_path} ${pipeline_root_dir} ${chemistry_schema}
        """
}
