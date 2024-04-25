
process GWASGALLOP {
  scratch true
  label 'medium'

  input:
    tuple val(fSimple), path(samplelist), path(rawfile)
    path x, stageAs: 'phenotypes.tsv'
    each phenoname

  output:
    tuple env(KEY), path("*.gallop")

  script:
    def m = []
    def outfile = rawfile.getName()
    m = outfile =~ /(.*).raw/
    outfile = "${m[0][1]}.${params.out}"

    def getkey = []
    def pop_pheno = samplelist.getName()
    getkey = pop_pheno =~ /(.*)_filtered.pca.tsv/
    pop_pheno = getkey[0][1]

    def model = ""
    def pheno_name = ""
    def pheno_name_file = ""

    if (params.model != '') {
      model = "--model '${params.model}'"
    }

    // if (params.pheno_name != '') {
    //   pheno_name = "--pheno-name '${params.pheno_name}'"
    // }

    // if (params.pheno_name_file != '') {
    //   pheno_name_file = "--pheno-name-file '${params.pheno_name_file}'"
    // }


    """
    set -x
    KEY="${pop_pheno}_${phenoname}"

    gallop --gallop \
           --rawfile ${rawfile} \
           --pheno "phenotypes.tsv" \
           --pheno-name "${phenoname}" \
           --covar ${samplelist} \
           --covar-name ${params.covariates} \
           --time-name ${params.time_col} \
           --out "${outfile}"
    """
}