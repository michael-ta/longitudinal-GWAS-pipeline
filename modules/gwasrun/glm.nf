
process GWASGLM {
  scratch true
  label 'medium'
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"

  input:
    tuple val(fSimple), path(plog), path(pgen), path(psam), path(pvar) //from p3_in_files_plink
    each samplelist //from plink_samplelist
  output:
    tuple env(KEY), path("*.linear")


  script:
    def covariates = "${params.covariates}".replaceAll(/ /, ",")
    def m = []
    def cohort = samplelist.getName()
    m = cohort =~ /(.*)_analyzed.tsv/
    cohort = m[0][1]
    
    def outfile = "${cohort}_${fSimple}.${params.out}"      

    def pheno_name = "y"
    if (params.pheno_name != '') {
      pheno_name = "${params.pheno_name}"
    }

    """
    set -x
    KEY="${cohort}_${fSimple}"

    plink2 --pfile ${fSimple} \
           --glm hide-covar omit-ref cols=+beta,+a1freq \
           --pheno ${samplelist} \
           --pheno-name ${pheno_name} \
           --covar ${samplelist} \
           --covar-name ${covariates} \
           --covar-variance-standardize \
           --keep ${samplelist} \
           --output-chr chrM \
           --mac ${params.minor_allele_ct} \
           --hwe 1e-6 \
           --threads ${task.cpus} \
           --memory ${task.memory.toMega()} \
           --out ${outfile}
    """
}
