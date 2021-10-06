# Nextflow Configuration

Nextflow pipelines look for a configuration file `nextflow.conf` at runtime to determine resource allocation and
define any environment variables required for processes within the script. The longitudinal GWAS pipeline 
provides a default configuration file for users to manage resource allocation and define environment
variables specifying the paths of the outputs and cache directory. Users can override the default configuration 
by specifying a custom configuation file with the `-params-file` option or specifying individual parameters via
the commandline `--something value` at runtime.

By default, each task is launched using the `gwas-pipeline` docker image and assigned 2 cpus. The process labels
`small` and `medium` are used to manage resource allocation between 2 types of tasks within the pipeline. `small`
tasks run in parallel each with access to 2 cpus and 20GB of working memory. This label is typically assigned to 
the model fitting jobs. `medium` tasks have access to 48 cpus and 225GB of working memory. These tasks benefit 
from access to multiple cores and are typically run consecutively to take advantage of the capcacity of the 
compute resources.

Below is the default `nextflow.config` file used by the longitudinal GWAS pipeline

```ini
process.container = 'gwas-pipeline'

// environment variables for GWAS pipeline
env {
  GWAS_RESOURCE_DIR = '/files'
  GWAS_OUTPUT_DIR = '/files/longGWAS_pipeline/results'
  GWAS_STORE_DIR = '/files/longGWAS_pipeline/cache'
  ADDI_QC_PIPELINE = '/usr/src/ADDI-GWAS-QC-pipeline/addi_qc_pipeline.py'
}

executor {
  name = 'local'
  cpus = 48
  memory = '480 GB'
}

docker {
  enabled = true
  temp = 'auto'
}

process {
  cpus = 2

  withLabel: small {
    cpus = 2
    memory = '20 GB'
  }

  withLabel: medium {
    cpus = 48
    memory = '225 GB'
  }
}
```
