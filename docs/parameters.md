# Commandline Parameters

Users can specify the following parameters when running the Nextflow pipeline

## --ancestry

This parameter is used to specify the ancestry group to perform the analysis on. Valid values for this parameter
are `['EUR', 'SAS', ...]`

## --assembly

This parameter is used to specify the genome assembly of the input genotyping files. The outputs of the pipeline
will all be in relation to hg38. If you supply genotypes from a different reference assembly specify one of the 
following options `['hg18', 'hg19']` so the positions can be lifted over to hg38. The default value is `'hg38'`

## --dataset

This parameter is used to specify an identifier for the input genotype files so that subsequent re-runs of the 
pipeline can use cached results from the variant filtering and ancestry + outlier detection steps. The default
value is `''` and it is highly recommended to set this parameter to avoid misusing cached results.

## --input_vcf

This parameter is used to specify the paths of the genotyping files. The pipeline requires genotypes in 
uncompressed (\*.vcf) or compressed (\*.vcf.gz) files. Genotype files should be shredded at the chromosome level
and each file should contain the chromosome number prefix with 'chr' case insensitive.

Acceptable filenames

```text
/path/to/dataset_prefix_chr*.vcf
/path/to/dataset_chr*_suffix.vcf.gz
/path/to/chr*.vcf
```

## --longitudinal_flag

This flag is used to specify the estimation of longitudinal associations via the GALLOP algorithm. If the option 
is not supplied then the pipeline performs a cross-sectional GWAS using plink2.

_Note:_ performing longitudinal analysis requires the inclusion of a `study_days` variable in the input 
phenotype. This variable should correspond to the timepoint (in days) since the start of the study at which the 
measurement for the outcome was taken. Initial measurements taken at baseline will have a `study_days` of 0.

## --out

This parameter is used to specify the output suffix of the files to distinguish results or re-runs of the 
pipeline.

## --pheno_file

## --r2thres

This parameter is used to filter out imputed geneotypes of low quality if the input genotyping files include 
imputed variants.


