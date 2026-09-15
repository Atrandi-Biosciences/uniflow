#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
nextflow.preview.output = true

// Include modules
include { FALCO              } from './modules/falco/main.nf'
include { validateParameters ; paramsHelp ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'


// Workflows entry points
include { INPUT_VALIDATION   } from "./workflows/input_validation.nf"
include { BARCODE            } from "./workflows/barcode.nf"
include { DEMULTIPLEXING     } from "./workflows/demultiplexing.nf"
include { DNA                } from "./workflows/dna.nf"
include { RNA                } from "./workflows/rna.nf"
include { SAMPLE             } from "./workflows/sample.nf"
include { LIBRARY            } from "./workflows/library.nf"
include { EXPERIMENT         } from "./workflows/experiment.nf"


workflow {

    main:
    // Create a map to hold runtime parameters for reproducibility and provenance
    def runtime_params = [:]

    def products = new groovy.json.JsonSlurper().parse(
        file("${projectDir}/assets/config/products.json", checkIfExists: true)
    )
    if (!products.containsKey(params.product_id)) {
        error("Unknown product_id: ${params.product_id}")
    }
    def chemistry_id = products["${params.product_id}"].chemistry_id
    def modalities = products["${params.product_id}"].modalities
    def counts_format = products["${params.product_id}"].counts_format
    // Pick and choose modalities to run based on the product configuration
    def run_gene_expression_ch = modalities.contains('gene_expression')
        ? channel.value(true)
        : channel.empty()
    def run_cross_processing_ch = modalities.contains('cross_processing')
        ? channel.value(true)
        : channel.empty()
    def run_amplicon_ch = modalities.contains('amplicon')
        ? channel.value(true)
        : channel.empty()
    def run_variant_ch = modalities.contains('variant')
        ? channel.value(true)
        : channel.empty()

    def barcode_num = products["${params.product_id}"].barcode_num

    runtime_params.chemistry_id = chemistry_id
    runtime_params.modalities = modalities
    runtime_params.input_params = params

    def metadataDir = file("${workflow.outputDir}/metadata")
    metadataDir.mkdirs()
    def paramsFile = file("${workflow.outputDir}/metadata/params.json")
    paramsFile.text = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(runtime_params)
    )
    param_path_ch = channel.value(paramsFile)

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info(paramsSummaryLog(workflow))

    def amplicon_fasta_file_ch = params.amplicon_fasta
        ? channel.value(file(params.amplicon_fasta, checkIfExists: true))
        : channel.empty()

    def star_index_folder_ch = params.star_index
        ? channel.value(file(params.star_index, checkIfExists: true))
        : channel.empty()

    def bc_long = "${projectDir}/assets/barcodes/B4_${barcode_num}_v01/long.csv"
    def whitelist_file_ch = channel.fromPath(file(bc_long, checkIfExists: true))
    def chemistries_path = file("${projectDir}/assets/chemistry/chemistry_defs.json", checkIfExists: true)

    def p = nextflow.file.FileHelper.asPath(params.input_csv)
    def firstLine = p.newReader().withReader { r -> r.readLine() }
    def header = firstLine.replaceAll(/\r$/, '').split(',') as List<String>
    def demultiplexingEnabled = header.contains('demultiplexing_indices')
    def forceCellsEnabled = header.contains('force_cells')
    def metrics_reference_ch = channel.value(file("${projectDir}/assets/metrics/metrics_reference.csv", checkIfExists: true))
    def multiqc_static_config_ch = channel.value(file("${projectDir}/assets/config/multiqc_config.yaml", checkIfExists: true))
    def multiqc_css_ch = channel.value(file("${projectDir}/assets/config/multiqc_report.css", checkIfExists: true))
    def internal_flag = channel.value(params.internal)

    def chemistries = new groovy.json.JsonSlurper().parse(file(chemistries_path, checkIfExists: true))

    if (!chemistries.containsKey(chemistry_id)) {
        error("Unknown chemistry_id: '${chemistry_id}'. Valid options are: ${chemistries.keySet().sort().join(', ')}")
    }
    def chemistry_file = file("${workflow.outputDir}/metadata/chemistry.json")

    chemistry_file.text = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(chemistries[chemistry_id])
    )


    // Generate sample Metas
    channel.fromList(samplesheetToList(params.input_csv, "assets/schema_input.json"))
        | map { meta, read1, read2 ->

            if (demultiplexingEnabled && meta.demultiplexing_indices) {
                meta.demultiplexing_indices = meta.demultiplexing_indices
            }
            else {
                meta.demultiplexing_indices = null
            }
            if (forceCellsEnabled && meta.force_cells) {
                meta.force_cells = meta.force_cells
            }
            else {
                meta.force_cells = null
            }

            [meta, read1, read2]
        }
        | map { meta, read1, read2 ->
            def library_id_r2 = (read2.getName() =~ /^(.*?)(?:[_\.][Rr]2(?:_\d{3})?|[_\.]2)(?:\.f(?:ast)?q\.gz)$/)
            def library_id_r1 = (read1.getName() =~ /^(.*?)(?:[_\.][Rr]1(?:_\d{3})?|[_\.]1)(?:\.f(?:ast)?q\.gz)$/)
            // if library_id_r1 is not the same as library_id_r2, throw an error
            if (library_id_r1[0][1] != library_id_r2[0][1]) {
                throw new IllegalArgumentException("Mismatched library IDs for read pair: ${read1.getName()} and ${read2.getName()}")
            }
            [meta + [library_id: library_id_r1[0][1]], read1, read2]
        }
        | set { libraries }



    if (params.disable_validation) {
        log.warn("Input validation is disabled. This is not recommended unless you are confident in your input CSV and FASTQ files.")
        inputs_validated = channel.value(true)
        ch_amplicon_similarity = channel.empty()
        
    }
    else {
        INPUT_VALIDATION(amplicon_fasta_file_ch, libraries, demultiplexingEnabled)
        inputs_validated = INPUT_VALIDATION.out.validation_done
        ch_amplicon_similarity = INPUT_VALIDATION.out.amplicon_similarity
        
    }

    whitelist_ch = channel.fromList(
            [
                file("${projectDir}/assets/barcodes/B4_${barcode_num}_v01/bcD_${barcode_num}.txt", checkIfExists: true),
                file("${projectDir}/assets/barcodes/B4_${barcode_num}_v01/bcC_${barcode_num}.txt", checkIfExists: true),
                file("${projectDir}/assets/barcodes/B4_${barcode_num}_v01/bcB_${barcode_num}.txt", checkIfExists: true),
                file("${projectDir}/assets/barcodes/B4_${barcode_num}_v01/bcA_${barcode_num}_wRT.txt", checkIfExists: true),
            ]
        )
        .collect()

    ch_multiqc_files = channel.empty()
    ch_sample_metrics_files = channel.empty()
    ch_library_metrics_files = channel.empty()


    // Run basic barcode processing for all libraries
    libraries
        | map { meta, _read1, read2 ->
            [meta.library_id, meta.library_type, read2]
        }
        | unique
        | set { read2_libraries_ch }

    libraries
        | map { meta, read1, read2 ->
            [meta.library_id, meta.library_type, read1, read2]
        }
        | unique
        | flatMap { library_id, library_type, read1, read2 ->
            [
                [library_id, library_type, read1],
                [library_id, library_type, read2],
            ]
        }
        | combine(inputs_validated)
        | FALCO

    BARCODE(
        read2_libraries_ch,
        inputs_validated,
        whitelist_file_ch,
        chemistry_file,
    )

    DEMULTIPLEXING(
        libraries,
        BARCODE.out.barcode_mapping,
        params.input_csv,
        inputs_validated,
        whitelist_file_ch,
        chemistry_file,
    )

    DEMULTIPLEXING.out.per_sample_barcode_mappings
        | transpose
        | set { sample_barcode_mappings }

    RNA(
        DEMULTIPLEXING.out.rna_libraries,
        star_index_folder_ch,
        whitelist_ch,
        sample_barcode_mappings,
        params.trimming_length,
        amplicon_fasta_file_ch,
        run_gene_expression_ch,
        run_cross_processing_ch,
        internal_flag,
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        RNA.out.metrics.map { _meta, _modality, it -> it.flatten() }
    )

    DNA(
        DEMULTIPLEXING.out.dna_libraries,
        amplicon_fasta_file_ch,
        sample_barcode_mappings,
        run_amplicon_ch,
        run_variant_ch,
    )


    ch_multiqc_files = ch_multiqc_files.mix(
        DNA.out.metrics.mix(DNA.out.internal_metrics).map { _meta, _modality, it -> it.flatten() }
    )

    // Every task of a process produces an identical versions.yml, so dedup on content
    // (not on path, which is per-workdir).
    ch_multiqc_files = ch_multiqc_files.mix(
        channel.topic("versions").map { it.text }.unique().collect().map { texts ->
            def parser = new org.yaml.snakeyaml.Yaml()
            def tools = new TreeMap<String, TreeSet<String>>()
            texts.each { text ->
                (parser.load(text) ?: [:]).each { _process, entries ->
                    entries.each { tool, version ->
                        def v = version?.toString()?.trim()
                        if (v) {
                            tools
                                .computeIfAbsent(tool as String) { new TreeSet<String>() }
                                .add(v)
                        }
                    }
                }
            }
            tools
                .collect { tool, versions ->
                    "${tool}:\n" + versions.collect { "  - \"${it}\"" }.join("\n")
                }
                .join("\n")
        }.filter { it }.map { "${it}\n" }.collectFile(name: "uniflow_mqc_versions.yml")
    )

    ch_h5ad_files = channel.empty()
        .mix(
            DNA.out.counts_h5ad,
            RNA.out.counts_h5ad,
        )

    // Gather internally gathered metrics
    ch_sample_metrics_files = ch_sample_metrics_files.mix(
        DNA.out.metrics_parquet,
        RNA.out.metrics_parquet,
    )

    ch_library_metrics_files = ch_library_metrics_files.mix(
        BARCODE.out.metrics_parquet
    )
    // Format and generate per sample outputs
    SAMPLE(
        ch_h5ad_files,
        ch_sample_metrics_files,
        counts_format,
        metrics_reference_ch,
        multiqc_static_config_ch,
        multiqc_css_ch,
        internal_flag
    )
    // Format and generate per library outputs
    LIBRARY(
        ch_library_metrics_files
    )

    // Aggregate metrics across samples and libraries
    ch_merged_metrics_csv = channel.empty()
        .mix(
            SAMPLE.out.per_sample_metrics,
            LIBRARY.out.merged_metrics_csv,
        )
    
    qc_merged_channels = channel.empty().mix(RNA.out.qc).mix(DNA.out.qc).mix(BARCODE.out.qc
    | map {library_id, modality, plots -> def meta = [
        sample_name: library_id
    ]
    return [meta, modality, plots]})
    
    // Generate aggregated reports for the entire experiment
    EXPERIMENT(
        ch_merged_metrics_csv,
        metrics_reference_ch,
        ch_multiqc_files,
        multiqc_static_config_ch,
        multiqc_css_ch,
        internal_flag,
        qc_merged_channels,
        param_path_ch,

    )

    publish:
    status               = SAMPLE.out.per_sample_status
    metadata             = DEMULTIPLEXING.out.metadata
    amplicon_similarity  = ch_amplicon_similarity
    sample_metrics       = SAMPLE.out.per_sample_metrics
    sample_reports       = SAMPLE.out.per_sample_report
    counts               = SAMPLE.out.counts
    fastqc_reports       = FALCO.out.tagged_html
    library_metrics      = LIBRARY.out.merged_metrics_csv
    barcode_qc           = BARCODE.out.qc
    corrected_barcodes   = BARCODE.out.barcode_mapping
    metrics              = channel.empty().mix(RNA.out.metrics).mix(DNA.out.metrics)
    internal_metrics     = channel.empty().mix(RNA.out.internal_metrics).mix(DNA.out.internal_metrics)
    qc                   = channel.empty().mix(RNA.out.qc).mix(DNA.out.qc)
    dev                  = channel.empty().mix(RNA.out.dev).mix(DNA.out.dev)
    internal_dev         = DNA.out.internal_dev
    rna_out              = RNA.out.out
    dna_out              = DNA.out.out
    multiqc_report       = EXPERIMENT.out.report
    experiment_metrics   = EXPERIMENT.out.aggregated_metrics_csv
    qc_tarball           = EXPERIMENT.out.qc_tarball
}


output {
    multiqc_report {
        path { report -> "experiment/" }
    }
    sample_metrics {
        path { sample_id, metrics_csv -> "samples/${sample_id}/metrics/" }
    }
    sample_reports {
        path { sample_id, _report -> "samples/${sample_id}" }
    }
    counts {
        path { meta, _files -> "samples/${meta.sample_name}/counts/" }
    }
    metadata {
        path { csv -> "metadata/" }
    }
    amplicon_similarity {
        enabled params.internal
        path { _csv -> "experiment/" }
    }
    fastqc_reports {
        path { library_id, library_type, html -> "libraries/${library_id}/reports" }
    }
    library_metrics {
        path { library_id, metrics_csv -> "libraries/${library_id}/metrics/" }
    }
    corrected_barcodes {
        enabled params.internal
        path { library_id, library_type, barcode_mappings -> "libraries/${library_id}/dev" }
    }
    barcode_qc {
        path { library_id, modality, plots -> "libraries/${library_id}/qc/${modality}/" }
    }
    rna_out {
        path { meta, files -> "samples/${meta.sample_name}/RNA/" }
    }
    metrics {
        path { meta, modality, files -> "samples/${meta.sample_name}/metrics/${modality}/" }
    }
    internal_metrics {
        enabled params.internal
        path { meta, modality, files -> "samples/${meta.sample_name}/metrics/${modality}/" }
    }
    qc {
        enabled params.internal
        path { meta, modality, plots -> "samples/${meta.sample_name}/qc/${modality}/" }
    }
    dev {
        path { meta, modality, files -> "samples/${meta.sample_name}/dev/${modality}/" }
    }
    internal_dev {
        enabled params.internal
        path { meta, modality, files -> "samples/${meta.sample_name}/dev/${modality}/" }
    }
    dna_out {
        path { meta, bam -> "samples/${meta.sample_name}/DNA/" }
    }
    experiment_metrics {
        path { metrics_csv -> "experiment/" }
    }
    status {
        path { sample_id, status_csv -> "samples/${sample_id}/" }
    }
    qc_tarball {
        path { tarball -> "experiment/" }
    }
}
