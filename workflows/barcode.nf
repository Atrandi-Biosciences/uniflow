include { EXTRACT_BARCODE ; CORRECT_BARCODE ; QC_BARCODE } from '../modules/barcode.nf'


def getLibraryLengths(Map chemistry) {
    chemistry.libraries.collectEntries { libraryName, library ->

        def length = library
            .collect { segment ->
                def part = segment.values().first()
                part.start + part.length
            }
            .max()

        [(libraryName): length]
    }
}
workflow BARCODE {
    take:
    read2_libraries_ch
    inputs_validated
    whitelist_file_ch
    chemistry_json_path

    main:



    def modality_name = "barcode"
    def chemistry_def = new groovy.json.JsonSlurper().parse(file(chemistry_json_path, checkIfExists: true))
    r2_lengths = getLibraryLengths(chemistry_def)


    EXTRACT_BARCODE(read2_libraries_ch, inputs_validated, channel.value(r2_lengths))

    EXTRACT_BARCODE.out.barcode_extracted.combine(whitelist_file_ch)
        | CORRECT_BARCODE

    CORRECT_BARCODE.out.corrected_barcodes
        | combine(whitelist_file_ch)
        | QC_BARCODE

    emit:
    barcode_mapping          = CORRECT_BARCODE.out.barcode_mapping
    metrics_parquet          = CORRECT_BARCODE.out.metrics_parquet
    qc                       = channel.empty().mix(QC_BARCODE.out.barcode_qc).groupTuple(by: 0).map { meta, files ->
        return [meta, modality_name, files]
    }
}
