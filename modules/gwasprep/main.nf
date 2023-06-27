/* Process 3 GWAS Model Preprocess + Fitting
 * -----------------------------------------
 * Get cohorts (study_arms) and subjects within each cohort to fit model
 * Compute PCs on a cohort level using Plink PC1 - PC10
 *
 * GWAS models available
 *  - Plink GLM for cross sectional GWAS
 *  - GALLOP for longitudinal GWAS
 *  - R Survival for cross-sectional & longitudinal survival GWAS
 */

// process GETPHENOS {
//   scratch true
//   label 'small'
  
//   input:
//     path covarfile, stageAs: 'covariates.tsv' //from "${params.covarfile}"
//   output:
//     path "phenos_list.txt", emit: allphenos//into p3_get_cohorts_processed
  
//   script:
//     """
//     #!/usr/bin/env python3
//     import pandas as pd
//     import time
    
//     study_id_colname = "${params.study_col}"

//     data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')
//     cohorts = data_df[study_id_colname].unique().tolist()
    
//     with open("phenos_list.txt", 'w') as f:
//       f.write('\\n'.join(cohorts))
    
//     time.sleep(5)
//     """
// }


// process REMOVEOUTLIERS {
//   scratch true
//   label 'medium'
//   storeDir "${STORE_DIR}/${params.dataset}/p3_COVARIATES_QC/${params.out}"

//   input:
//     path samplelist //from p2_qc_processed
//     path covarfile, stageAs: 'covariates.tsv' //from "${params.covarfile}"
//     each cohort //from p3_in_cohort_list 
//   output:
//     path "${params.ancestry}_${cohort}_filtered.tsv" //into gwas_samplelist

//   script:
//     //println "Separating IID based on ancestry (${ancestry})"
//     """
//     #!/usr/bin/env python3
//     import pandas as pd
//     import time

//     ancestry = "${params.ancestry}"
//     study_id_colname = "${params.study_col}"

//     ancestry_df = pd.read_hdf("${samplelist}", key="ancestry_keep")
//     outlier_df = pd.read_hdf("${samplelist}", key="outliers")
//     kin_df = pd.read_hdf("${samplelist}", key="kin")
//     data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')

//     cohorts = data_df[study_id_colname].unique().tolist()
//     cohorts = filter(lambda x: x == "${cohort}", cohorts)

//     kin_df = kin_df[kin_df.KINSHIP >= ${params.kinship}]
//     # TODO: address case when single cohort present
//     # TODO: currently does not address longitudinal covariates
//     for cohort in cohorts:
//       print(f'---- {cohort} ----')
//       samples = data_df[ (data_df.IID.isin(ancestry_df.IID)) &
//                          (data_df[study_id_colname] == cohort) &
//                         ~(data_df.IID.isin(outlier_df.IID)) ].copy(deep=True)

//       r = kin_df[(kin_df['#IID1'].isin(samples.IID)) & (kin_df.IID2.isin(samples.IID))].copy()
//       samples = samples[~samples.IID.isin(r.IID2)].copy()
//       samples.to_csv(f"{ancestry}_{cohort}_filtered.tsv", sep="\t", index=False)
//       print(f'Samples removed (outliers) = {data_df.IID.isin(outlier_df.IID).sum()}')
//       print(f'Samples removed (kinship) = {r.shape[0]}')
//       print(f'Samples remaining = {len(samples)}')
//       print('')

//     time.sleep(5)
//     """
// }


// process COMPUTE_PCA {
//   scratch true
//   label 'large_mem'
   
//   storeDir "${STORE_DIR}/${params.dataset}/p3_PCA_QC/"
//   publishDir "${OUTPUT_DIR}/${params.dataset}/MERGER_SPLITS_LOGS_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  
//   input:
//     each file(samplelist) //from gwas_samplelist.flatten()
//     file "*" //from input_p3_pca.collect()
    
//   output:
//     tuple file(samplelist), file("${cohort_prefix}.pca.eigenvec"), emit: eigenvec //into p3_cohort_pca_processed
//     //tuple file("*.log") file("*.eigenval") //into p3_cohort_pca_logs
//     //file "*.eigenval", emit: eigenvalues //into p3_cohort_pca_eval
//     //file "*.log", emit: logs

//   script:
//     def m = []
//     def cohort = ""
//     cohort = samplelist.getName()
//     m = cohort =~ /(.*)_filtered.tsv/
//     cohort = m[0][1]
//     cohort_prefix = "${cohort}"
    
//     """
//     plink2 \
//           --indep-pairwise 50 .2 \
//           --maf ${params.minor_allele_freq} \
//           --pfile "allchr_${params.dataset}_p2in" \
//           --out ${cohort}.ld

//     plink2 \
//           --keep ${samplelist} \
//           --out ${cohort}.pca \
//           --extract ${cohort}.ld.prune.in \
//           --pca 10 \
//           --threads ${task.cpus} \
//           --memory ${task.memory.toMega()} \
//           --pfile "allchr_${params.dataset}_p2in"
//     """
// }


// process MERGE_PCA {
//   scratch true
//   label 'small'
  
//   input:
//     tuple file(samplelist), file(cohort_pca) //from p3_cohort_pca_processed
  
//   output:
//     file "${params.ancestry}_*_filtered.pca.tsv" //into p3_merge_pca_processed
  
//   script:
//     """
//      #!/usr/bin/env python3
//      import pandas as pd
//      import time   
//      import os
     
//      print(os.listdir())
//      sample_fn = "${samplelist.getName()}"
//      cohort = sample_fn[:-(len('_filtered.tsv'))]
//      pc_fn = cohort + '.pca.eigenvec'

//      print(pc_fn, cohort)
//      pc_df = pd.read_csv(pc_fn, sep="\t")
//      samples_df = pd.read_csv(sample_fn, sep="\t")
//      pc_df.rename(columns={"#IID": "IID"}, inplace=True)
//      samples_df = samples_df.merge(pc_df, on="IID")
     
//      samples_df.to_csv(cohort + '_filtered.pca.tsv', sep="\t", index=False)
//      time.sleep(5)
//      """
// }

// /* todo
//  * collect input files and process GWAS with gallop or plink
//  */


// p3_merge_pca_processed
//   .into { gwas_samplelist_plink; gwas_samplelist_gallop }

// p1_run_processed_ch
//   .into{ p1_run_processed_ch_new; p1_run_processed_ch_old }

// p1_run_processed_ch_old
//   .groupTuple(by: 0)
//   .flatten()
//   .collate(5)
//   .into { p3_in_files_gallop; p3_in_files_plink }
  
// p1_run_processed_ch_new
//   .groupTuple(by: 0)
//   .flatten()
//   .filter( ~/.*pvar/ )
//   .map { it -> tuple( it.getSimpleName(), it ) }
//   .set{ p3_get_plink_chunks_pfile_input }


// process GALLOPCOX_INPUT {

//   scratch true
//   label 'small'
  
//   input:
//     tuple val(chrname), path(plink_input) //from p3_get_plink_chunks_pfile_input
//   output:
//     tuple val(chrname), path("allchr_${params.dataset}_p2in.txt") //into p3_plink_chunks
    
//    when:
//     params.longitudinal_flag || params.survival_flag
  
//   script:
//     """
//     #!/usr/bin/env python3
//     fn = "${plink_input}"
//     out_fn = "allchr_${params.dataset}_p2in.txt"
//     count = 0
//     id_pairs = []
//     start, end =  None,None
//     with open(fn, 'r') as f:
//       for l in iter(f.readline, ''):
//         if l[0] == '#':
//           continue
//         data = l.strip().split('\t')
//         vid = data[2]
//         count += 1
//         if start is None:
//           start = vid

//         if count >= ${params.plink_chunk_size}:
//           end = vid
//           id_pairs.append( (start, end) )
//           start = None
//           end = None
//           count = 0

//     if count > 0:
//       end = vid
//       id_pairs.append( (start, end) )

//     with open(out_fn, 'w') as f:
//       for start,end in id_pairs:
//         f.write( '\\t'.join([start,end]) + '\\n' )
//     """
// }

// p3_plink_chunks = p3_plink_chunks
//                     .splitText()
                    
// p3_plink_chunks = p3_in_files_gallop
//                   .cross(p3_plink_chunks)
//                   .flatten()
//                   .collate(7)
//                   .map{ it -> tuple(it[0], it[1], it[2], it[3], it[4], it[6]) } 


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
           --from \$from \
           --to \$to \
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

// gwas_rawfile
//   .into{ gwas_rawfile_gallop; gwas_rawfile_coxph }