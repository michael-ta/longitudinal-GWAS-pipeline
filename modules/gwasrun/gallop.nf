
process GWASGALLOP {
  scratch true
  label 'medium'

  //publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    tuple val(fSimple), path(samplelist), path(rawfile) //from gwas_rawfile_gallop
    path x, stageAs: 'phenotypes.tsv' //from "${params.phenofile}"
  
  output:
    tuple env(KEY), path("*.gallop") //into gallop_results
  
  script:
    def m = []
    def cohort = rawfile.getName()
    m = cohort =~ /(.*).raw/
    outfile = "${m[0][1]}.${params.out}"


    def model = ""
    def pheno_name = ""
    def pheno_name_file = ""

    if (params.model != '') {
      model = "--model '${params.model}'"
    }
    
    if (params.pheno_name != '') {
      pheno_name = "--pheno-name '${params.pheno_name}'"
    }

    if (params.pheno_name_file != '') {
      pheno_name_file = "--pheno-name-file '${params.pheno_name_file}'"
    }

    """
    set -x
    KEY="${cohort}_${fSimple}"
    gallop --gallop \
           --rawfile ${rawfile} \
           --pheno "phenotypes.tsv" \
           --covar ${samplelist} \
           --covar-name ${params.covariates} \
           ${model} ${pheno_name} ${pheno_name_file} \
           --time-name ${params.time_col} \
           --out "${outfile}"
    """
}