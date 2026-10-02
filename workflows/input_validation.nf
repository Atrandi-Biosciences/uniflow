include { FASTA_VALIDATION } from '../subworkflows/fasta_validation/main.nf'


workflow INPUT_VALIDATION {
    take:
    amplicon_reference_fasta_ch
    libraries
    demultiplexingEnabled

    main:

    // The amplicon FASTA is validated against the FASTA specification by
    // bin/validate_amplicon_fasta.py
    FASTA_VALIDATION(amplicon_reference_fasta_ch)

    // Validate that meta.library_id does not contain R1 or R2. Error out if it does.
    libraries
        | map { meta, _fastq1, _fastq2 ->
            if (meta.library_id =~ /R[12]/) {
                error("Invalid library_id: ${meta.library_id}\nPlease provide a valid filename that does not contain R1 or R2.")
            }
        }
    // Check that each library channel and their fastq files is unique. If not, error out.
    libraries
        .map { meta, fastq1, fastq2 ->
            tuple(
                meta.sample_name,
                meta.library_type,
                fastq1,
                fastq2,
                meta.demultiplexing_indices,
            )
        }
        .collect(flat: false)
        .map { rows ->

            def seen = [:]
            def duplicates = []

            rows.eachWithIndex { row, idx ->

                def processingKey = tuple(
                    row[1],
                    row[2],
                    row[3],
                    row[4],
                )

                if (seen.containsKey(processingKey)) {
                    duplicates << [
                        line: idx + 2,
                        sample_name: row[0],
                        original_line: seen[processingKey].line,
                        original_name: seen[processingKey].sample_name,
                        key: processingKey,
                    ]
                }
                else {
                    seen[processingKey] = [
                        line: idx + 2,
                        sample_name: row[0],
                    ]
                }
            }

            if (duplicates) {
                error(
                    "Duplicate samplesheet inputs found:\n" + duplicates.collect { d ->
                        "Line ${d.line} ('${d.sample_name}') duplicates " + "line ${d.original_line} ('${d.original_name}'):\n" + "  library_type=${d.key[0]}\n" + "  fastq_1=${d.key[1]}\n" + "  fastq_2=${d.key[2]}\n" + "  demultiplexing_indices=${d.key[3]}"
                    }.join('\n\n')
                )
            }

            true
        }

    // Compare demultiplexing indices across library types for each sample. If they do not match, error out.
    if (demultiplexingEnabled) {
        libraries
            | map { meta, _fastq1, _fastq2 ->
                [meta.sample_name, meta.library_type, meta.demultiplexing_indices]
            }
            | groupTuple(by: [0])
            | map { sample_name, _library_type, indices ->
                if (indices.unique().size() > 1) {
                    error("Demultiplexing indices do not match across library types for sample: ${sample_name}\n${indices}\nPlease provide valid demultiplexing indices that match across library types.")
                }
            }
        libraries
            .flatMap { meta, fastq1, fastq2 ->
                def indices = meta.demultiplexing_indices ? meta.demultiplexing_indices.tokenize(';') : [null]
                indices.collect { index ->
                    tuple(meta.sample_name, meta.library_type, fastq1, fastq2, index)
                }
            }
            .collect(flat: false)
            .map { rows ->
                def seen = [:]
                def duplicates = []

                rows.each { row ->
                    def key = tuple(row[1], row[2], row[3], row[4])
                    if (seen.containsKey(key)) {
                        duplicates << "  '${row[0]}' collides with '${seen[key]}' on index ${row[4] ?: 'none'}"
                    }
                    else {
                        seen[key] = row[0]
                    }
                }

                if (duplicates) {
                    error("Duplicate demultiplexing inputs:\n" + duplicates.join('\n'))
                }
            }
    }

    emit:
    validation_done      = FASTA_VALIDATION.out.done
    amplicon_similarity  = FASTA_VALIDATION.out.similarity
    fasta_validation_log = FASTA_VALIDATION.out.report
}
