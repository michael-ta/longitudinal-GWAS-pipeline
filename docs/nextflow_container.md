# Nextflow Container

The GWAS pipeline can be run from a containerized Nextflow installation using a Docker image built from the 
example Dockerfile below or similar. This Dockerfile installs the necessary Java runtime environment and grabs 
the latest Nextflow binary from the official website. Create a Dockerfile by copying the script below.

```Dockerfile
FROM ubuntu:18.04

RUN set -eux; \
        apt-get update; \
        apt-get install -y --no-install-recommends \
          openjdk-8-jre-headless \
          wget \
        && wget -qO- https://get.nextflow.io | bash \
        && chmod +x nextflow \
        && mv nextflow /usr/local/bin
```

You can build this image using the Docker `build` command below if the above Dockerfile is created in the 
current user directory. Otherwise, set the path of the Dockerfile in the build command.

```bash
sudo docker build -t local:nextflow .
```

Once this image has been built, create a nextflow container with access to all the user mounts needed to run
the pipeline along with the docker socket and binary required for Nextflow to launch containers. If running this 
from the ADWB platform, the following configuration is needed to mount the user home directory, file share, and
docker dependencies. Other mount points required to access analysis input or output files can be supplied using 
the `-v` option. 

```sh
sudo docker run -it -d \
  -v /usr/bin/docker:/usr/bin/docker \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v /files:/files \
  -v $HOME:$HOME \
  -v /tmp:/tmp \
  --name nf local:nextflow
```

The GWAS pipeline can then be run with the following command. Please note, the full path to the input files must
be provided for the pipeline to work in this manner. Relative paths may not be evaluated and passed correctly to
the working containers by the Dockerized Nextflow installation.

```sh
sudo docker exec -w $PWD nf \
  nextflow gwas-pipeline.nf -w /files/longGWAS_pipeline/work \
    --input_vcf '/files/ADNI/example/example.chr*.vcf' \
    --phenofile /files/ADNI/example/phenotype.lt.tsv \
    --longitudinal_flag --covarfile /files/ex.covars.tsv \
    --out ex_new-lt \
    --dataset ex_new-dock-nf \
    --covariates 'SEX PC1 PC2 PC3' \
    --assembly hg19
```

