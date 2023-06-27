
process RAWFILE_EXPORT {
  scratch true
  label 'small'
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"

  input:
    tuple val(fSimple), path(plog), path(pgen), path(psam), path(pvar), path(plink_chunk) //from p3_plink_chunks
    each path(samplelist) //from gwas_samplelist_gallop.flatten()
  
  output:
    tuple val(fSimple), path(samplelist), path('*.raw'), emit: gwas_rawfile //into gwas_rawfile 
    path "*.log", emit: gwas_rawfile_log //into p3_export_rawfile_log
  
  when:
    params.longitudinal_flag || params.survival_flag

  script:
    def cohort = ""
    cohort = samplelist.getName()
    m = cohort =~ /(.*)_filtered.pca.tsv/
    cohort = m[0][1]

    def outfile = "${cohort}_${fSimple}"

    """
    set -x
    from=\$(cat $plink_chunk | cut -f 1)
    to=\$(cat $plink_chunk | cut -f 2)
    
    plink2 --pfile ${fSimple} \
           --keep ${samplelist} \
           --export A \
           --from \${from} \
           --to \${to} \
           --mac ${params.minor_allele_ct} \
           --update-sex ${samplelist} \
           --pheno ${samplelist} \
           --pheno-col-nums 4 \
           --hwe 1e-6 \
           --out "${outfile}"  \
           --threads ${task.cpus} \
           --memory ${task.memory.toMega()}
    """
}



process EXPORT_PLINK {
  debug true
  scratch true
  label 'small'

  input:
    val samplelist //from gwas_samplelist_plink.flatten()
    path x, stageAs: 'phenotypes.tsv' //from "${params.phenofile}"
  output:
    path "*_analyzed.tsv" optional true //into plink_samplelist
  when:
    ! params.longitudinal_flag && ! params.survival_flag

  script:
    //println "Samplelist: ${samplelist}"
    //println "Phenotypes: phenotypes.tsv"

    def m = []
    def cohort = ""
    cohort = samplelist.getName()
    m = cohort =~ /(.*)_filtered.pca.tsv/
    outfile = "${m[0][1]}"
    
    def pheno_name = "y"
    if (params.pheno_name != '') {
      pheno_name = "${params.pheno_name}"
    }

    """
    #!/usr/bin/env python3
    import pandas as pd
    import time
    import sys

    covars = "${params.covariates}"
    covars = covars.split(' ')

    d_pheno = pd.read_csv("phenotypes.tsv", sep="\t", engine='c')
    d_sample = pd.read_csv("${samplelist}", sep="\t", engine='c')

    d_result = pd.merge(d_pheno, d_sample, on='IID', how='inner')
    
    if d_result.shape[0] > 0:    
      d_set = d_result.loc[:, ["#FID", "IID", "${pheno_name}"] + covars].copy()
      d_set.to_csv("${outfile}_analyzed.tsv", sep="\t", index=False)

    time.sleep(10)
    """
}