The metrics module is used throughout the pipeline to collect metrics into a parquet file to be later aggregated and provided as an output.

To make it simple, each metrics output provides the same output structure:

Example:
```
sample_id: sampleA
library_type: RNA
process_name: barcode_correction
metrics:
  number_of_A_barcodes_corrected: 3049
  number_of_B_barcodes_corrected: 3157
  number_of_C_barcodes_corrected: 3179
  number_of_D_barcodes_corrected: 2728
  number_of_primers_corrected: 3700
  total_reads: 51155042
  corrected_filtered_reads: 46360770
  discarded_reads: 4794272
```

The top level keys are identifier keys.

* sample_id: sequencing library name
* library_type: what kind of modality is in the library
* process_name: Name of the process that populated the metrics. This is automatically uses the name of the file where the MetricsRecord is being instantiated.
* metrics: any metrics produced in the process.

Here is an example on how to use it for a process located at: `bin/PROCESS_NAME.py`

In the code use: 

```
from lib.metrics.metrics import MetricsRecord

metrics = MetricsRecord(sample_id=SAMPLE_ID, library_type=LIBRARY_TYPE)

metrics.add(NAME_OF_METRIC, VALUE_OF_METRIC)

metrics.write_records("metrics.parquet")
```

This is then emitted in the nextflow process as:
` path("metrics.parquet"), emit: metrics_parquet`