# Input File Format

## Input Genotype

The pipeline requires genotype information to be given in the form of uncompressed VCF or compressed VCF.gz 
files.

If missing variants have been imputed, specify the optional parameter `--r2thres` to filter out poorly imputed
variants

## Input Covariates

Covariates for each subject to be passed into the model can be provided via a tab-delimited file (\*.tsv)

For both cross-sectional and longitudinal analysis, the pipeline expects covariates to be defined in the 
following format

_Note:_ the Plink style columns `#FID` and `PHENO` can be populated with 0

```text
#FID  IID SEX PHENO study_arm apoe4 levodopa_usage age_at_baseline
0 sid-1 1 0 control 0 0 35
0 sid-2 1 0 control 0 0 40
0 sid-3 0 0 control 1 0 32
.
.
.
0 sid-98  1 0 PD  1 0 55
0 sid-99  0 0 PD  0 1 66
0 sid-100 1 0 PD  0 0 58
```


## Input Phenotype / Outcomes

Phenotype and measured outcomes can be passed into the pipeline via a tab-delimited file (\*.tsv)

For cross-sectional analysis, the pipeline expects a minimum of 2 columns in the following format

```text
IID y
sid-1 1
sid-2 0
sid-3 1
.
.
.
sid-98 0
sid-99 0
sid-100 1
```

