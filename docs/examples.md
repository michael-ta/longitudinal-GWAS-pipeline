# Examples

The following commands can be used to run the example dataset provided in the github repository. 
First untar the example dataset

```sh
tar -xzf example/genotype/example.vcf.tar.gz
```

Second run the pipeline for cross-sectional analysis using the following command

```sh
sudo nextflow gwas-pipeline.nf \
  --input_vcf "example/genotype/*.vcf.gz" \
  --phenofile example/phenotype.cs.tsv \
  --covarfile example/covariates.tsv \
  --assembly "hg19"
```

Run the pipeline for longitudinal analysis

```sh
sudo nextflow gwas-pipeline.nf \
  --input_vcf "example/genotype/*.vcf.gz" \
  --phenofile example/phenotype.lt.tsv \
  --covarfile example/covariates.tsv \
  --assembly "hg19" \
  --longitudinal_flag
```

