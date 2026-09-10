### Pipeline Testing with nf-test

This repository uses [nf-test] https://www.nf-test.com/docs to ensure that modules, subworkflows, and the entire pipeline work correctly through automated unit and integration tests.

All test files are located in the `/tests/` folder.  
The main test configuration (`nf-test.config`) is in the `/uniflow/` folder. This file is required by the nf-test framework.  
You can also load additional configurations. The `nextflow.config` file located in the `/tests/` folder will be loaded automatically.

---

## Prerequisites

Before running the tests, make sure the following tools are installed:

* **Nextflow** (version 25.04.8)
* **nf-test** (version 0.9.5)
* **Docker** (required for container environments)

---

### Running Tests

Tests are executed using the `nf-test` command-line tool. By default, tests use the Nextflow profile defined in `nf-test.config`.

```bash
# Run All Tests
# To run the complete test suite (modules, subworkflows, and end-to-end tests):
nf-test test

# Run a Single Test
# You can run only the specific test file you are currently working on:

# Test a single module
nf-test test tests/modules/some_module.nf.test

# Test a single subworkflow
nf-test test tests/subworkflows/some_subworkflow.nf.test

# Run all process tests
nf-test test --tag process 

# Run all workflow tests
nf-test test --tag workflow

# Run local tests before committing. Important!
nf-test test --only-changed --follow-dependencies
--only-changed
#When enabled, this parameter instructs nf-test to execute tests only for files that have been modified within the current git working tree.
--follow-dependencies
# When this flag is set, nf-test will traverse all dependencies when the related-tests flag is set. This option is particularly useful when you need to ensure that all  dependent tests are executed.

# Switch profiles
# If you want to run the tests a specific environment (for example, docker)
nf-test test --profile docker
```

---

## Folder structure

The sctructure inside the `/tests/` folder mirrors the layout of our pipeline. This makes tests easy and quick to find_
Die Struktur im `/tests/`-Ordner spiegelt den Aufbau unserer Pipeline wider, um Tests schnell auffindbar zu machen:

tests/
├── config/                 # Central test configurations and nf-test profiles
│   └── nextflow-test.config     # Specific nf-test settings for test data
├── modules/                # Unit tests for individual processes
│   └── some_module.nf.test
├── subworkflows/           # Integration tests for combined module workflow
│   └── some_subworkflow.nf.test
├── pipeline/               # End-to-End (E2E) tests for the entire pipeline
|   └── main.nf.test
├───params.json             # Paths for entering test data
|
test_data/
└───modules/                # Input of test data for process tests

### Test data
## Unit Tests
All unit tests (process tests) use their own test data as input. The test data is very small to ensure that the large number of tests runs quickly. This minimal test data is sufficient to verify that the processes, as individual components, run without errors.

## Integration Tests
All integration tests retrieve test data from S3. The paths to the input data are defined in the params.json file. The test data is larger compared to the test data of the unit tests, which makes the tests slightly slower.

---

### Create new test

For every new module and every new subworkflow, a corresponding `.nf.test` file **must** be created

## Automatically generate a test skeleton
Use the nf-test CLI to create a standardized test template for an existing module:
```bash
# Process test
nf-test generate process modules/some_module.nf

# Workflow test
nf-test generate workflow subworkflows/some_workflow.nf
```

### Best practices for assertions (test expectations)
For stable tests, always use snapshots instead of fixed file paths, MD5 checksums (which often change due to tool updates), or timestamps.

Example of a clean test file (`tests/modules/some_module.nf.test`):

```groovy
nextflow_process {

    name "Test for the module SOME_MODULE"
    script "modules/some_module.nf"
    process "SOME_MODULE"

    test("Should create correct output files with default parameters") {

        when {
            process {
                """
                // samples_ch
                input[0] = Channel.fromList([
                    [
                        [sample_name: 'sample1', library_id: 'Lib1_RNA_L001', library_type: 'RNA'],
                        file("${path}/Lib1_RNA_L001_R1_001.fastq.gz")
                    ]
                """
            }
        }

        then {
            // Check if the Nextflow process completed successfully
            assert process.success

            // Create an check a snapshot of the entire output (directories, files, channels)
            assert snapshot(process.out).match()
        }

    }
}
```

---

## Updating snapshots

If you intentionally change the source code of a module (for example, a tool update or additional output parameters), the old snapshots will fail. You can automatically update them like this:

```bash
# Update snapshots for all tests
nf-test test --update-snapshot

# Update the snapshot for a specific test
nf-test test tests/modules/some_module.nf.test --update-snapshot
```

## EC2 instance
For end-to-end testing, the EC2 instance i-056f6547c9b51393b must be started. The service `/ect/systemd/system/start_runner.service` automatically starts the Docker container for Bitbucket after the instance starts.

In Bitbucket, under Repository Settings - Pipelines - Runners, the status should be "online".

After the tests start, the directory 61335a8d-xxx etc. is created under `~/atlassian-bitbucket-pipelines-runner/temp/`. The repository is cloned in this directory, and the tests are executed there.

After the test, the directory is automatically cleaned up. If the instance hangs and fails to clean up, the instance status in Bitbucket is set to "unhalthy".
In this case, the entire directory 61335a8d-xxx etc. must be deleted manually. Connect to the instance via SSH and delete the directory.
After the directory has been deleted, the service must be restarted with `sudo systemctl stop start_runner`
`sudo systemctl start start_runner`.