#!/usr/bin/env nextflow

include { UPDATE_DEFINITIONS } from './modules/chemistry.nf'


nextflow.preview.output = true

workflow {
    
    main:
        def barcode_config_path_ch = channel.value(file("${projectDir}/assets/chemistry/chemistry_defs.yaml", checkIfExists: true))
        def pipeline_root_dir = channel.value(file("${projectDir}", checkIfExists: true))
        def chemistry_schema_path = channel.value(file("${projectDir}/assets/chemistry/chemistry_schema.json", checkIfExists: true))
        // Update barcode definitions
        UPDATE_DEFINITIONS(barcode_config_path_ch, pipeline_root_dir, chemistry_schema_path)
    
    publish:
        updated_config = UPDATE_DEFINITIONS.out.updated_config

}

output {
    updated_config{
        path { json -> "${projectDir}/assets/chemistry/" }
    }
        
}