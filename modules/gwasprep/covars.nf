
process COMPUTE_PCA {
  scratch true
  label 'large_mem'

  storeDir "${STORE_DIR}/${params.dataset}/p3_PCA_QC/"
  publishDir "${OUTPUT_DIR}/${params.dataset}/LOGS/COMPUTE_PCA_${params.datetime}/", mode: 'copy', overwrite: true, pattern: "*.log"

  input:
    each file(samplelist) //from gwas_samplelist.flatten()
    file "*" //from input_p3_pca.collect()
    
  output:
    tuple file(samplelist), file("${cohort_prefix}.pca.eigenvec"), emit: eigenvec //into p3_cohort_pca_processed
    //tuple file("*.log") file("*.eigenval") //into p3_cohort_pca_logs
    //file "*.eigenval", emit: eigenvalues //into p3_cohort_pca_eval
    //file "*.log", emit: logs

  script:
    def m = []
    def cohort = ""
    cohort = samplelist.getName()
    m = cohort =~ /(.*)_filtered.tsv/
    cohort = m[0][1]
    cohort_prefix = "${cohort}"
    
    """
    plink2 \
          --indep-pairwise 50 .2 \
          --maf ${params.minor_allele_freq} \
          --pfile "allchr_${params.dataset}_p2in" \
          --out ${cohort}.ld

    plink2 \
          --keep ${samplelist} \
          --out ${cohort}.pca \
          --extract ${cohort}.ld.prune.in \
          --pca 10 \
          --threads ${task.cpus} \
          --memory ${task.memory.toMega()} \
          --pfile "allchr_${params.dataset}_p2in"
    """
}


process MERGE_PCA {
  scratch true
  label 'small'
  
  input:
    tuple file(samplelist), file(cohort_pca) //from p3_cohort_pca_processed
  
  output:
    file "${params.ancestry}_*_filtered.pca.tsv" //into p3_merge_pca_processed
  
  script:
    """
     #!/usr/bin/env python3
     import pandas as pd
     import time   
     import os
     
     print(os.listdir())
     sample_fn = "${samplelist.getName()}"
     cohort = sample_fn[:-(len('_filtered.tsv'))]
     pc_fn = cohort + '.pca.eigenvec'

     print(pc_fn, cohort)
     pc_df = pd.read_csv(pc_fn, sep="\t")
     samples_df = pd.read_csv(sample_fn, sep="\t")
     pc_df.rename(columns={"#IID": "IID"}, inplace=True)
     samples_df = samples_df.merge(pc_df, on="IID")
     
     samples_df.to_csv(cohort + '_filtered.pca.tsv', sep="\t", index=False)
     time.sleep(5)
     """
}