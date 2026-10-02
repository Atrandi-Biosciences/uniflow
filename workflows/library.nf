include { MERGE_METRICS as MERGE_LIBRARY_METRICS } from '../modules/shared.nf'

workflow LIBRARY {
    take:
    ch_library_metrics_files

    main:
    ch_library_metrics_files
        | groupTuple(by: 0)
        | MERGE_LIBRARY_METRICS

    emit:
    merged_metrics_csv = MERGE_LIBRARY_METRICS.out.metrics_csv
}
