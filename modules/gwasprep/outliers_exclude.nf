
process GETPHENOS {
  scratch true
  label 'small'
  
  input:
    path covarfile, stageAs: 'covariates.tsv' //from "${params.covarfile}"
  output:
    path "phenos_list.txt", emit: allphenos//into p3_get_cohorts_processed
  
  script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import time
    
    study_id_colname = "${params.study_col}"

    data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')
    cohorts = data_df[study_id_colname].unique().tolist()
    
    with open("phenos_list.txt", 'w') as f:
      f.write('\\n'.join(cohorts))
    
    time.sleep(5)
    """
}


process REMOVEOUTLIERS {
  scratch true
  label 'medium'
  storeDir "${STORE_DIR}/${params.dataset}/p3_COVARIATES_QC/${params.out}"

  input:
    path samplelist //from p2_qc_processed
    path covarfile, stageAs: 'covariates.tsv' //from "${params.covarfile}"
    each cohort //from p3_in_cohort_list 
  output:
    path "${params.ancestry}_${cohort}_filtered.tsv" //into gwas_samplelist

  script:
    //println "Separating IID based on ancestry (${ancestry})"
    """
    #!/usr/bin/env python3
    import pandas as pd
    import time

    ancestry = "${params.ancestry}"
    study_id_colname = "${params.study_col}"

    ancestry_df = pd.read_hdf("${samplelist}", key="ancestry_keep")
    outlier_df = pd.read_hdf("${samplelist}", key="outliers")
    kin_df = pd.read_hdf("${samplelist}", key="kin")
    data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')

    cohorts = data_df[study_id_colname].unique().tolist()
    cohorts = filter(lambda x: x == "${cohort}", cohorts)

    kin_df = kin_df[kin_df.KINSHIP >= ${params.kinship}]
    # TODO: address case when single cohort present
    # TODO: currently does not address longitudinal covariates
    for cohort in cohorts:
      print(f'---- {cohort} ----')

      if not outlier_df.empty:
        samples = data_df[ (data_df.IID.isin(ancestry_df.IID)) &
                           (data_df[study_id_colname] == cohort) &
                           ~(data_df.IID.isin(outlier_df.IID)) ].copy(deep=True)
        print(f'Samples removed (outliers) = {data_df.IID.isin(outlier_df.IID).sum()}')
      else:
        samples = data_df[ (data_df.IID.isin(ancestry_df.IID)) &
                           (data_df[study_id_colname] == cohort)].copy(deep=True)

      r = kin_df[(kin_df['#IID1'].isin(samples.IID)) & (kin_df.IID2.isin(samples.IID))].copy()
      samples = samples[~samples.IID.isin(r.IID2)].copy()
      samples.to_csv(f"{ancestry}_{cohort}_filtered.tsv", sep="\t", index=False)
      print(f'Samples removed (kinship) = {r.shape[0]}')
      print(f'Samples remaining = {len(samples)}')
      print('')

    time.sleep(5)
    """
}