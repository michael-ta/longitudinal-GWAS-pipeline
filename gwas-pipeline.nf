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

params.longitudinal_flag = false
params.survival_flag = false
params.chunk_flag = false
params.mh_plot = false

params.parallel_splits = 3
parallel_splits = params.parallel_splits > 5 ? 5 : params.parallel_splits
parallel_splits = params.parallel_splits < 1 ? 1 : params.parallel_splits


datetime = new java.util.Date()
params.datetime = new java.text.SimpleDateFormat("YYYY-MM-dd'T'HHMMSS").format(datetime)


/* Process 1 Run - Cache handling
 * ------------------------------
 * Attempt to map the parameters from `params.input` and `params.dataset` to stored cached files only re-process 
 * `params.input` files if a cache output does not exist
 * 
 * Instead of relying on `storeDir` we can reduce the need to cache intermediate files for each process
 * - sets `input_p1_run_ch` to empty if all outputs needed for the next step is already in the cache
 *   otherwise, re-run `p1_run` process given the missing outputs
 */

Channel
  .fromPath( "${STORE_DIR}/${params.dataset}/p1_run_cache/*" )
  .map{ f -> tuple(f.getSimpleName(), f) }
  .into{ cache_check_ch; input_cache_ch }

Channel
  .fromPath(params.input)
  .map{ f -> tuple(f.getSimpleName(), f) }
  .set{ input_check_ch}

input_check_ch
  .join(cache_check_ch, remainder: true)
  .filter{ fSimple, fOrig, fCache -> fCache == null }
  .map{ fSimple, fOrig, fCache -> tuple( fSimple, fOrig ) }
  .into{ input_basename_ch; input_p1_count_ch; input_split_ch}


input_p1_count_ch
  .count()
  .map{ val -> (1..val) }
  .flatten()
  .set{ input_idx_list }


input_idx_list
  .merge(input_split_ch)
  .branch {
     split_tasks_5: parallel_splits > 4 && it[0] % 5
     split_tasks_4: parallel_splits > 3 && it[0] % 4
     split_tasks_3: parallel_splits > 2 && it[0] % 3
     split_tasks_2: parallel_splits > 1 && it[0] % 2
     
     split_tasks_1: true
   }.set { input_p1 }

input_p1.split_tasks_5
  .map{ fIdx, fSimple, fOrig -> fOrig }
  .splitText(by: params.chunk_size, file: true, compress: true)
  .map{ fn -> tuple( fn.getSimpleName(), fn) }
  .set{ input_split_ch_5 }

input_p1.split_tasks_4
  .map{ fIdx, fSimple, fOrig -> fOrig }
  .splitText(by: params.chunk_size, file: true, compress: true)
  .map{ fn -> tuple( fn.getSimpleName(), fn) }
  .set{ input_split_ch_4 }

input_p1.split_tasks_3
  .map{ fIdx, fSimple, fOrig -> fOrig }
  .splitText(by: params.chunk_size, file: true, compress: true)
  .map{ fn -> tuple( fn.getSimpleName(), fn) }
  .set{ input_split_ch_3 }

input_p1.split_tasks_2
  .map{ fIdx, fSimple, fOrig -> fOrig }
  .splitText(by: params.chunk_size, file: true, compress: true)
  .map{ fn -> tuple( fn.getSimpleName(), fn) }
  .set{ input_split_ch_2 }

input_p1.split_tasks_1
  .map{ fIdx, fSimple, fOrig -> file(fOrig) }
  .splitText(by: params.chunk_size, file: true, compress: true)
  .map{ fn -> tuple( fn.getSimpleName(), fn) }
  .set{ input_split_ch_1 }

input_split_ch_1
  .mix(input_split_ch_2, 
       input_split_ch_3,
       input_split_ch_4,
       input_split_ch_5)
  .set{ input_chunks_ch }

input_basename_ch
  .cross(input_chunks_ch)
  .flatten()
  .collate(4, false)
  .map{ fSimple1, fOrig, fSimple2, fSplit -> 
    tuple( fSimple1, file(fOrig), fSplit) }
  .set{ input_p1_run_ch }


/* Process 1 Run - Definition
 * -------------
 * Genotype preprocessing
 * returns the SNPs after applying
 *  - pass (&R2) filter
 *  - split
 *  - left-normalize
 *  - autosmal-par
 *  - hg38 alignment (liftOver)
 *  - mac >= 2
 *  - geno < 0.05
 * 
 * Inputs:
 *  `input_p1_run_ch` -> tuple( vSimple: val, fOrig: file, fSplit: file )
 *    vSimple - simple name of input file (everything before first '.')
 *    fOrig   - original input file
 *    fSplit  - chunked file to process

 * Outputs:
 *  `p1_run_processed_ch` -> tuple( fPgen: file, fPvar: file, fPsam: file )
 *    process chunked plink files
 *  `p1_run_chunklist_ch` -> tuple( val, val )
 *    chunck prefix - group by simple name
 */
   
process p1_run {
  scratch true
  
  label 'small'

  input:
    set val(vSimple), file(fOrig), file(fSplit) from input_p1_run_ch
  output:
    tuple file("${output}.pgen"), file("${output}.pvar"), file("${output}.psam") into p1_run_processed_ch
    tuple val(vSimple), val("${output}") into p1_run_chunklist_ch

  script:
  def chrnum = ""
  def vPart = ""
  def prefix = ""

  // regex match the chromosome corresponding to the input file chunk from the simple filename - requires that
  // each input file have chr0-9 in the filename; also require the user not to include any periods in the
  // filename
  mChrom = vSimple =~ /(?i)(chr)([0-9]+)/
  chrnum = mChrom[0][2]

  vPart = fSplit.getBaseName()
  // basename should contain the filename without the compression (.gz) extension - requires the input file to
  // have some sort of .extension in the filename. We can regex match split part from the remaining filename
  mPart = vPart =~ /(.*)\.([0-9]+)\.(.*)$/
  vPart = mPart[0][2]

  prefix = "${vSimple}.${vPart}"
  output = "${prefix}_p1out"
  ext = fSplit.getExtension()

  """
  echo "Processing - ${fSplit}"
  echo "Assigned cpus: ${task.cpus}"
  echo "Assigned memory: ${task.memory}"
  
  set +x
  if [[ ${vPart} == 1 ]]; then
    cp $fSplit tmp_input.${ext}
  else
    bcftools view -h $fOrig | gzip > tmp_input.${ext}
    cat $fSplit >> tmp_input.${ext}
  fi

  /srv/GWAS-Pipeline/References/Scripts/process1.sh \
    ${task.cpus} \
    tmp_input.${ext} \
    ${params.r2thres} \
    ${params.assembly} \
    ${chrnum} \
    ${prefix}
  """
}

/* Process 1 Run - Chunk Merging
 * ------------------------------
 * Merge the processed chunked files from `p1_run` back into single corresponding plink output
 * 
 * Prepare data for merging all individual input files into single file for `p2_qc_pipeline`
 */

p1_run_chunklist_ch
  .collectFile( newLine: true ) { vSimple, prefix ->
    ["${vSimple}.mergelist.txt", prefix] }
  .into{ input_p1_merge_chunk_ch; input_p1_merge_nochunk_ch }
  
p1_run_processed_ch
  .into{ p1_run_processed_ch_chunk; p1_run_processed_ch_nochunk; tmp_view }


process p1_merge_chunk {
  scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p1_run_cache"
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  
  label 'small'

  input:
    file mergelist from input_p1_merge_chunk_ch
    file "*" from p1_run_processed_ch_chunk.collect()
  output:
    tuple file("${vSimple}.psam"), file("${vSimple}.pgen"), file("${vSimple}.pvar") into p1_merge_chunk_processed_ch
    file "*.log" into p1_merge_chunk_log
   when:
     chunk_flag

  script:
  vSimple = mergelist.getSimpleName()
  """
  # run merge command on tmp list
  set +x
  
  plink2 --pmerge-list ${mergelist} \
    --make-pgen \
    --sort-vars \
    --out ${vSimple}
  plink2 --make-bed \
    --pfile ${vSimple} \
    --out ${vSimple}
  """
}

process p1_merge_nochunk {
  scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p1_run_cache"
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  
  label 'small'

  input:
    file mergelist from input_p1_merge_nochunk_ch
    file "*" from p1_run_processed_ch_nochunk.collect()
  output:
    tuple file("${vSimple}.psam"), file("${vSimple}.pgen"), file("${vSimple}.pvar") into p1_merge_nochunk_processed_ch
    file "*.log" into p1_merge_nochunk_log
  when:
    ! chunk_flag

  script:
  vSimple = mergelist.getSimpleName()
  """
  # run merge command on tmp list
  set +x
  
  plink2 -pfile ${vSimple}.1_p1out \
    --make-pgen \
    --sort-vars \
    --out ${vSimple}
  """
}

p1_merge_chunk_processed_ch
  .mix(p1_merge_nochunk_processed_ch)
  .collect()
  .set{p1_merge_processed_ch}

p1_merge_processed_ch
  .flatten()
  .map{ fn -> tuple(fn.getSimpleName(), fn) }
  .concat(input_cache_ch)
  .into{ p1_run_processed_ch; tmp_input_p2_merge_list; tmp_input_p2_merge_list_files }


/* Process 2 QC Pipeline - Merge All
 * ------------------------------
 * Merge the input files (should be shredded by chromosome) from `p1_run` back into single plink1/2 output
 * for downstream consumption by `p2_qc_pipeline`
 * 
 */

tmp_input_p2_merge_list
  .collectFile() { vSimple, f ->
    [ "allchr.mergelist.txt", f.getBaseName() + '\n' ] }
  .set{ input_p2_merge_list }

tmp_input_p2_merge_list_files
  .map{ vSimple, f -> file(f) }
  .set{ input_p2_merge_list_files }

plink_prefix = "allchr_${params.dataset}_p2in"


process p2_merge_list {
  scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p2_merged_cache"
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  
  label 'large_mem'

  input:
    file mergelist from input_p2_merge_list
    file "*" from input_p2_merge_list_files.collect()
  output:
    file ("${plink_prefix}.{bed,fam,bim,pgen,pvar,psam}") into input_p2_merged_plink
    file "*.log" into p2_mergelist_log

  script:
    """
    set -x
    cat $mergelist | uniq > tmp_mergefile.txt
    plink2 --memory ${task.memory.toMega()} \
      --pmerge-list "tmp_mergefile.txt" \
      --make-bed \
      --out "${plink_prefix}"
    """
}

input_p2_merged_plink
  .flatten()
  .into{ input_p2_qc_pipeline; temp_p3_pca}

temp_p3_pca
  .filter{ fName -> ["pgen", "pvar", "psam"].contains( fName.getExtension() ) }
.set{ input_p3_pca }


/* Process 2 QC Pipeline - Definition
 * ------------------------------
 * Run GWAS QC Pipeline
 *  - ancestry assignment
 *  - outlier pruning
 *  - LD structure pruning
 *  - kinship relationship 
 * 
 */

process p2_qc_pipeline {
  //scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p2_qc_pipeline_cache"
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/plots", mode: 'copy', overwrite: true, pattern: "*.{html,png}"
  
  label 'large_mem'
  
  input:
    file "*" from input_p2_qc_pipeline.collect()
  output:
    file "${params.ancestry}_samplelist_p2out.h5" into p2_qc_processed
    file "*.{html,png}" into p2_out_figs
  
  script:
    """
    set +x

    
    python3 ${ADDI_QC_PIPELINE} \
      --geno "${plink_prefix}" \
      --ref "/srv/GWAS-Pipeline/References/ref_panel/1kg_ashkj_ref_panel_gp2_pruned_hg38_newids" \
      --ref_labels "/srv/GWAS-Pipeline/References/ref_panel/ancestry_ref_labels.txt" \
      --pop "${params.ancestry}" \
      --out "${params.ancestry}_samplelist_p2out"    
    """
}

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

process p3_get_cohorts {
  scratch true
  label 'small'
  
  input:
    path covarfile, stageAs: 'covariates.tsv' from "${params.covarfile}"
  output:
    file "cohort_list.txt" into p3_get_cohorts_processed
  
  script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import time
    
    study_id_colname = "${params.study_col}"

    data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')
    cohorts = data_df[study_id_colname].unique().tolist()
    
    with open("cohort_list.txt", 'w') as f:
      f.write('\\n'.join(cohorts))
      
    time.sleep(5)
    """
}

p3_get_cohorts_processed
  .readLines()
  .set{ p3_in_cohort_list }


process p3_cohort_samplelist {
  scratch true
  label 'medium'
  
  storeDir "${STORE_DIR}/${params.dataset}/p3_cohort_pca_cache/${params.out}"

  input:
    file samplelist from p2_qc_processed
    path covarfile, stageAs: 'covariates.tsv' from "${params.covarfile}"
    each cohort from p3_in_cohort_list 
  output:
    file "${params.ancestry}_${cohort}_filtered.tsv" into gwas_samplelist

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
      samples = data_df[ (data_df.IID.isin(ancestry_df.IID)) &
                         (data_df[study_id_colname] == cohort) &
                        ~(data_df.IID.isin(outlier_df.IID)) ].copy(deep=True)

      r = kin_df[(kin_df['#IID1'].isin(samples.IID)) & (kin_df.IID2.isin(samples.IID))].copy()
      samples = samples[~samples.IID.isin(r.IID2)].copy()
      samples.to_csv(f"{ancestry}_{cohort}_filtered.tsv", sep="\t", index=False)
      print(f'Samples removed (outliers) = {data_df.IID.isin(outlier_df.IID).sum()}')
      print(f'Samples removed (kinship) = {r.shape[0]}')
      print(f'Samples remaining = {len(samples)}')
      print('')

    time.sleep(5)
    """
}


process p3_cohort_pca {
  scratch true
  label 'large_mem'
   
  storeDir "${STORE_DIR}/${params.dataset}/p3_cohort_pca_cache/${params.out}"
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.{log,eigenval}"
   publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.{eigenval}"
  
  input:
    each file(samplelist) from gwas_samplelist.flatten()
    file "*" from input_p3_pca.collect()
    
  output:
    tuple file(samplelist), file("${cohort_prefix}.pca.eigenvec") into p3_cohort_pca_processed
    file "*.log" into p3_cohort_pca_logs
    file "*.eigenval" into p3_cohort_pca_eval
   
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
          --pfile "${plink_prefix}" \
          --out ${cohort}.ld

    plink2 \
          --keep ${samplelist} \
          --out ${cohort}.pca \
          --extract ${cohort}.ld.prune.in \
          --pca 10 \
          --threads ${task.cpus} \
          --memory ${task.memory.toMega()} \
          --pfile "${plink_prefix}"
    """
}

process p3_merge_pca {
  scratch true
  label 'small'
  
  input:
    
    set file(samplelist), file(cohort_pca) from p3_cohort_pca_processed
    
  output:
    file "${params.ancestry}_*_filtered.pca.tsv" into p3_merge_pca_processed
   
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

/* todo
 * collect input files and process GWAS with gallop or plink
 */


p3_merge_pca_processed
  .into { gwas_samplelist_plink; gwas_samplelist_gallop }

p1_run_processed_ch
  .into{ p1_run_processed_ch_new; p1_run_processed_ch_old }

p1_run_processed_ch_old
  .groupTuple(by: 0)
  .flatten()
  .collate(5)
  .into { p3_in_files_gallop; p3_in_files_plink }
  
p1_run_processed_ch_new
  .groupTuple(by: 0)
  .flatten()
  .filter( ~/.*pvar/ )
  .map { it -> tuple( it.getSimpleName(), it ) }
  .set{ p3_get_plink_chunks_pfile_input }


process p3_get_plink_chunks {
  scratch true
  label 'small'
  
  input:
    tuple val(plink_prefix), file(plink_input) from p3_get_plink_chunks_pfile_input
  output:
    tuple val(plink_prefix), file("${plink_prefix}.txt") into p3_plink_chunks
    
   when:
    params.longitudinal_flag || params.survival_flag
  
  script:
    """
#!/usr/bin/env python3
fn = "${plink_input}"
out_fn = "${plink_prefix}.txt"
count = 0
id_pairs = []
start, end =  None,None
with open(fn, 'r') as f:
  for l in iter(f.readline, ''):
    if l[0] == '#':
      continue
    data = l.strip().split('\t')
    vid = data[2]
    count += 1
    if start is None:
      start = vid

    if count >= ${params.plink_chunk_size}:
      end = vid
      id_pairs.append( (start, end) )
      start = None
      end = None
      count = 0

if count > 0:
  end = vid
  id_pairs.append( (start, end) )


with open(out_fn, 'w') as f:
  for start,end in id_pairs:
    f.write( '\\t'.join([start,end]) + '\\n' )
    """
}

p3_plink_chunks = p3_plink_chunks
                    .splitText()
                    
p3_plink_chunks = p3_in_files_gallop
                  .cross(p3_plink_chunks)
                  .flatten()
                  .collate(7)
                  .map{ it -> tuple(it[0], it[1], it[2], it[3], it[4], it[6]) } 


process p3_export_rawfile {
  scratch true
  label 'small'
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"

  input:
    set val(fSimple), file(plog), file(pgen), file(psam), file(pvar), file(plink_chunk) from p3_plink_chunks
    each file(samplelist) from gwas_samplelist_gallop.flatten()
  output:
    tuple val(fSimple), file(samplelist), file('*.raw') into gwas_rawfile 
    file "*.log" into p3_export_rawfile_log
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

gwas_rawfile
  .into{ gwas_rawfile_gallop; gwas_rawfile_coxph }

process p3_gwas_coxph {
  scratch true
  label 'small'

  //publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    set val(fSimple), file(samplelist), file(rawfile) from gwas_rawfile_coxph
    path x, stageAs: 'phenotypes.tsv' from "${params.phenofile}"
  output:
    tuple env(KEY), file("*.coxph") into coxph_results
  when:
    params.survival_flag && ! params.longitudinal_flag
  
  script:
    def m = []
    def cohort = rawfile.getName()
    m = cohort =~ /(.*).raw/
    outfile = "${m[0][1]}.${params.out}.coxph"

    def pheno_name = ""

    if (params.pheno_name != '') {
      pheno_name = "--pheno-name '${params.pheno_name}'"
    }

    """
    set -x
    KEY="${cohort}_${fSimple}"
    /srv/GWAS-Pipeline/References/Scripts/survival.R \
               --rawfile ${rawfile} \
               --pheno "phenotypes.tsv" \
               --covar ${samplelist} \
               --covar-name "${params.covariates}" \
               ${pheno_name} \
               --out ${outfile}
    """
}


process p3_gwas_gallop {
  scratch true
  label 'medium'

  //publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    set val(fSimple), file(samplelist), file(rawfile) from gwas_rawfile_gallop
    path x, stageAs: 'phenotypes.tsv' from "${params.phenofile}"
  output:
    tuple env(KEY), file("*.gallop") into gallop_results
  when:
    params.longitudinal_flag && ! params.survival_flag
  
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


// send the output to plink only if there are subjects to test <- determined by
// merging the covariates and outcome files on IID

process p3_format_gwas_plink {
  echo true
  scratch true
  label 'small'

  input:
    val samplelist from gwas_samplelist_plink.flatten()
    path x, stageAs: 'phenotypes.tsv' from "${params.phenofile}"
  output:
    file "*_analyzed.tsv" optional true into plink_samplelist
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


process p3_gwas_plink{
  scratch true
  label 'medium'
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"

  input:
    set val(fSimple), file(plog), file(pgen), file(psam), file(pvar) from p3_in_files_plink
    each samplelist from plink_samplelist
  output:
    file "*.linear" into gwas_results
    file "*.log" into p3_log
  when:
    ! params.longitudinal_flag

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

/* Process 3 collect and export results
 * ------------------------------
 * Collect shredded files and export GWAS results on chromosome level
 * 
 */

p3_processed = gallop_results.mix(coxph_results)
              .groupTuple()
              .map{ it[1] }
              .flatten()
              .collectFile(keepHeader: true)
              
process save_outputs {
  scratch true
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true
  
  input:
    file(sumstats) from p3_processed
    
  output:
    file(sumstats) into sumstats_out
    
  script:
    """
    bedtools sort -i ${sumstats} -header > "${sumstats}".tmp
    mv ${sumstats}.tmp ${sumstats}
    """
}

/* Process 3 manhattan plot (needs update)
 * ------------------------------
 * create manhattan plot and qq plot for gwas summary stats
 * 
 */
 
process p3_gwas_viz {
  scratch true
  label 'medium'

  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/plots", mode: 'copy', overwrite: true

  input:
    file x from gwas_results.mix(sumstats_out).collect()
  output:
    file "*.png" into gwas_plots
  when:
    params.mh_plot

  script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    from qmplot import manhattanplot
    from qmplot import qqplot
    from multiprocessing import Pool

    files = list(filter(lambda x: os.path.splitext(x)[-1] in ['.gallop', '.linear'], os.listdir()))
    dfs = []
  
    # group by phenotype
    lt_flag = "${params.longitudinal_flag}" == "true"
    suffix = "${params.out}"
    threads = int("${task.cpus}")

    def plot_summary_stats(data, cohort, outcome, lt_flag=False):
      xtick = set(['chr' + i for i in list(map(str, range(1, 10))) + ['11', '13', '15', '18', '21', 'X']])
      if lt_flag:
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="Pi", ax = ax, 
                      xtick_label_set=xtick)
        plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.gallop.png", dpi=300)
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["Pi"],
                marker="o",
                title=f"QQ Intercept {cohort} {outcome}",
                xlabel=r"Expected -log(P)",
                ylabel=r"Observed -log(P)",
                ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_intercept.gallop.png", dpi=300)
  
      else:
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="P", ax = ax,
                      xtick_label_set=xtick)
        plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.gallop.png", dpi=300)
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_intercept.gallop.png", dpi=300)

    plot_df = dict()
    for f in files:
      fc = f.split('.')
      pheno = fc[-2]
      cohort = '_'.join(fc[0].split('_')[:-1])
      try:
        plot_df[f'{cohort}+{pheno}'].append(f)
      except KeyError:
        plot_df[f'{cohort}+{pheno}'] = [f]


    def read_table(fn):
      'converts a filename to a pandas dataframe'
      df = None
      try:
        df = pd.read_table(fn, sep="\t")
      except pd.errors.EmptyDataError:
        df = pd.DataFrame()
      return df


    for cohort in plot_df.keys():
      df = None
      c, outcome = cohort.split('+')
      
      with Pool(processes=threads) as pool:

        # have your pool map the file names to dataframes
        df_list = pool.map(read_table, plot_df[cohort])

        # reduce the list of dataframes to a single dataframe
        df = pd.concat(df_list, ignore_index=True)
      
      df = df.dropna(how="any", axis=0)  # clean data
      
      if df.shape[0] == 0:
        continue

      df['chr_order'] = df['#CHROM'].str.replace('chr', '')
      df['chr_order'] = df['chr_order'].astype(int)
      df = df.sort_values(by=['chr_order', 'POS'])

      # generate manhattan plot and set an output file.
      plot_summary_stats(data=df, cohort=f'{c}.{suffix}', outcome=outcome, lt_flag=lt_flag)
    """
}
