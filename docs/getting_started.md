# Getting Started

This page will walk you through the process of building the docker image and setting up Nextflow to run the 
GWAS pipeline. If running this pipeline from ADWB, skip the installation section and go to the basic usage
section.

## Installation

### Requirements

The pipeline requires the following software prerequisites and has been tested on a Linux Ubuntu 18.04 environment

- Nextflow 
- Docker
- Java JRE 8

To enable Nextflow to manage memory resources within Docker containers set the following option in
`/etc/default/grub`

```ini
GRUB_CMDLINE_LINUX="cgroup_enable=memory swapaccount=1"
```

### Clone the Repository

The repository for the pipeline can be found at 
[https://github.com/michael-ta/longitudinal-GWAS-pipeline](https://github.com/michael-ta/longitudinal-GWAS-pipeline).
Run the command below to clone the repository

```sh
git clone https://github.com/michael-ta/longitudinal-GWAS-pipeline.git 
```

### Build the Docker Image

This repository comes with the `Dockerfile` and local resources needed to build the docker image used by the 
pipeline. After cloning the repository, you can build the image with

```sh
cd longitudinal-GWAS-pipeline
sudo docker build --build-arg BUILD_VAR=$(date +%Y%m%d-%H%M%S) -t gwas-pipeline .
```

The parameter `BUILD_VAR` sets the environment variable `IMAGE_BUILD_VAR` to the date and time of the current 
build. This can be used to track different versions of the Docker image once built. Using the default 
`nextflow.config` the pipeline will launch containers using a local image with the label `gwas-pipeline`. This 
behavior can be adjusted by changing the `nextflow.config` or setting the option at runtime. For more details, 
see the Nextflow Configuration page.

_Note:_ invoking sudo is not necessary if the Docker user has been previously added to sudoers

## Basic Usage

Once the Docker image is built, the pipeline can be called by running the following command within the local 
cloned irepository

```sh
sudo nextflow gwas-pipeline.nf \
  --input_vcf "example/data/genetic/*.vcf.gz" \
  --covarfile "examples/basic/covar.tsv" \
  --phenofile "examples/basic/pheno.tsv"
```

The outputs from the pipeline will be saved to the directory defined by the environment variable 
`GWAS_OUTPUT_DIR` in the Nextflow configuration. Within the output directory you'll see the following folders
and files after running the pipeline on the basic example files

```text
results/
  cor_timestamp/
    *.linear
    plots/
      *.

cache/
  p1_run_cache/
  p2_qc_pipeline_cache/ 
```
