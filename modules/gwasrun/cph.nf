process GWASCPH {
  scratch true
  label 'medium'

  input:
    tuple val(fSimple), path(samplelist), path(rawfile)
    path x, stageAs: 'phenotypes.tsv'
    each phenoname

  output:
    tuple env(KEY), path("*.coxph")
  
  script:
    def m = []
    def outfile = rawfile.getName()
    m = outfile=~ /(.*).raw/
    outfile = "${m[0][1]}.${params.out}.coxph"

    def getkey = []
    def  pop_pheno = samplelist.getName()
    getkey = pop_pheno =~ /(.*)_filtered.pca.tsv/
    pop_pheno = getkey[0][1]

    // def pheno_name = ""
    // if (params.pheno_name != '') {
    //   pheno_name = "--pheno-name '${params.pheno_name}'"
    // }

    """
    set -x
    KEY="${pop_pheno}_${phenoname}"
    
    survival.R --rawfile ${rawfile} \
               --pheno "phenotypes.tsv" \
               --covar ${samplelist} \
               --covar-name "${params.covariates}" \
               --pheno-name "${phenoname}" \
               --out ${outfile}
    """
}