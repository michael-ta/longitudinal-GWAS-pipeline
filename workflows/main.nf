#!/usr/bin/env nextflow 


/*
 * Enables modules
 */
nextflow.enable.dsl = 2


/*
 * Default pipeline parameters. 
 * To override these, we are going to use the params.yml file 
 */
params.input = "data/*"
params.dataset = ""
params.ancestry = "EUR"
params.assembly = "hg38"
params.out = ""
params.r2thres = -9
params.chunk_size = 30000
params.plink_chunk_size = 20000

params.covarfile = ""
params.study_col = "study_arm"
params.time_col = "study_days"
params.minor_allele_freq = "0.05"
params.minor_allele_ct = "20"
params.kinship = "0.177"

params.model = ""
params.pheno_name = ""
params.pheno_name_file = ""
params.covariates = "SEX PC1 PC2 PC3"
params.phenofile = ""

params.longitudinal_flag = false
params.survival_flag = false
params.linear_flag = false
params.chunk_flag = false
params.mh_plot = false


/*
 * Main workflow log
 */
if ( params.longitudinal_flag) {
    MODEL = "LMM GALLOP"
} 
else if ( params.survival_flag ) {
    MODEL = "Cox model"
}
else {
    MODEL = "Cox model"
}

log.info """\
 LONG-GWAS - GWAS P I P E L I N E
 ======================================
 Chunk size for genetic processing        : ${params.chunk_size}
 Kinship matrix threshold                 : ${params.kinship}
 R2 threshold                             : ${params.r2thres}
 MAF threshold                            : ${params.minor_allele_freq}
 data ancestry                            : ${params.ancestry}
 genetic data assemble                    : ${params.assembly}
 phenotype name                           : ${params.pheno_name}
 covariates                               : ${params.covariates}
 analysis                                 : ${MODEL}
 outdir                                   : ${PWD}/files/longGWAS_pipeline/results/cache/${params.dataset}
 """

 
/*
 * Datetime
 */
datetime = new java.util.Date()
params.datetime = new java.text.SimpleDateFormat("YYYY-MM-dd'T'HHMMSS").format(datetime)


/* 
 * Import modules
 */
include { DOQC }              from '../subworkflows/fullqc.nf'
include { GWASDATA_PREP }     from '../subworkflows/gwasinputs.nf'
include { GWAS_RUN }          from '../subworkflows/rungwas.nf'
include { SAVE_RESULTS }      from '../subworkflows/saveresults.nf'


/* 
 * Get the cache and the input check channels
 */
Channel
  .fromPath( "${STORE_DIR}/${params.dataset}/p1_run_cache/*" )
  .map{ f -> tuple(f.getSimpleName(), f) }
  .set{ cache }

Channel
   .fromPath(params.input)
   .map{ f -> tuple(f.getSimpleName(), f) }
   .set{ input_check_ch }


/* 
 * main script flow
 */
workflow {
  
  DOQC( input_check_ch, cache )

  GWASDATA_PREP(params.covarfile, DOQC)
  
  GWAS_RUN(params.phenofile, 
          GWASDATA_PREP.out.CHUNKS, 
          GWASDATA_PREP.out.PLINK_SLIST)

  SAVE_RESULTS(GWAS_RUN.out)
}