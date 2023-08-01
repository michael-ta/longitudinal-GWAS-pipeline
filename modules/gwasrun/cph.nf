process GWASCPH {
  scratch true
  label 'small'

  input:
    tuple val(fSimple), path(samplelist), path(rawfile) //from gwas_rawfile_coxph
    path x, stageAs: 'phenotypes.tsv' //from "${params.phenofile}"

  output:
    tuple env(KEY), path("*.coxph") //into coxph_results
    //path("*.coxph")
  script:
    def m = []
    def cohort = rawfile.getName()
    m = cohort =~ /(.*).raw/
    outfile = "${m[0][1]}.${params.out}.coxph"

    def getkey = []
    def uniqpheno = samplelist.getName()
    getkey = uniqpheno =~ /(.*)_filtered.pca.tsv/
    //println getkey

    def pheno_name = ""
    if (params.pheno_name != '') {
      pheno_name = "--pheno-name '${params.pheno_name}'"
    }

    """
    set -x
    #KEY="${cohort}_${fSimple}"
    KEY="${getkey[0][1]}"
    /srv/GWAS-Pipeline/References/Scripts/survival.R \
               --rawfile ${rawfile} \
               --pheno "phenotypes.tsv" \
               --covar ${samplelist} \
               --covar-name "${params.covariates}" \
               ${pheno_name} \
               --out ${outfile}
    """
}