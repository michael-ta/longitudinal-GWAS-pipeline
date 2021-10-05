#!/usr/bin/env nextflow

// Parameters needed for Process 1
params.input_vcf = "data/*.vcf.gz"
params.r2thres = -9
params.chr = ""
params.assembly = "hg38"
params.out = "cor"
params.dataset = ""

// Parameters needed for Process 2
// specify EUR or SAS
params.ancestry = "EUR"

// Parameters needed for Process 3
params.study_id = "study_arm"
params.covarfile = "data/covariates.tsv"
params.phenofile = "data/phenotypes.longitudinal.tsv"
params.kinship = "0.177"
params.missing_geno_rate = "0.05"

params.covariates = "SEX PC1 PC2 PC3"
params.model = ""
params.pheno_name = ""
params.pheno_name_file = ""

datetime = new java.util.Date()
params.datetime = new java.text.SimpleDateFormat("YYYY-MM-dd'T'HHMMSS").format(datetime)

params.longitudinal_flag = false

vfiles = Channel
              .fromPath(params.input_vcf)
              .map { file -> tuple( file.getName(), file.size(), file) }

// Define filepaths for resources needed for pipeline
//local_p1_script = "${GWAS_RESOURCE_DIR}/References/Scripts/process1.sh"
//local_assembly = "${GWAS_RESOURCE_DIR}/References/Genome/hg38.fa.gz"
//local_chain = "${GWAS_RESOURCE_DIR}/References/liftOver/${params.assembly}ToHg38.over.chain.gz"
//local_ref_labels = "${GWAS_RESOURCE_DIR}/References/ref_panel/ancestry_ref_labels.txt"
//local_ref_panel = "${GWAS_RESOURCE_DIR}/References/ref_panel/1kg_ashkj_ref_panel_gp2_pruned_hg38_newids.*"

//ref_panel_files = Channel
//                      .fromPath(local_ref_panel)

/* 
   --
   Process 1 - genotype preprocessing
   returns the pass (&R2) filtered, split, left-normalized, autosomal-par, hg38 aligned SNPs with mac >=2 & geno <0.05
   refer to process1.sh script for implementation details
   --
*/

process p1_run {
  scratch true

  memory { dataFile.getExtension() == 'gz' ? 4.GB * ((dataSize >> 30) + 1) :  2.GB * ((dataSize >> 30) + 1) }
  storeDir "${GWAS_STORE_DIR}/p1_run_cache"
  
  input:
    //path script, stageAs: 'script.sh' from local_p1_script
    //path assembly, stageAs: 'hg38.fa.gz' from local_assembly
    //path chainfile, stageAs: "${params.assembly}ToHg38.over.chain.gz" from local_chain
    set dataID, dataSize, file(dataFile) from vfiles
  output:
    // output filenames from process 1 will have the format "chr${chrnum}_${params.dataset}_p1out.*"
    file("${plink_prefix}.psam") into p1_out_psam
    file("${plink_prefix}.pgen") into p1_out_pgen
    file("${plink_prefix}.pvar") into p1_out_pvar
    val plink_prefix into mergelist_shred
    file "*.log" into p1_log
    
  script:
    def m = []
    def chrnum = ""
    def output = ""

    // case insensitive matching to get the chromosome number from input parameter or filename
    if (params.chr == "")
      m = dataID =~ /(?i)(chr)([0-9]+)/
    else
      m = params.chr =~ /(?i)(chr)([0-9]+)/

    // if matching fails - output error for user
    if (m.size() != 1 || m[0].size() < 3)
      error "Failed to identify chromosome from filename ${dataFile}, please rename or use the --chr argument"
    
    chrnum = m[0][2]
    output = "chr${chrnum}_${params.dataset}"
    plink_prefix = "${output}_p1out"    

    """
    echo "Processing - ${dataFile}"
    echo "Assigned cpus: ${task.cpus}"

    set +x
    /srv/GWAS-Pipeline/References/Scripts/process1.sh ${task.cpus} \
      ${dataFile} ${params.r2thres} ${params.assembly} ${chrnum} ${output}
    """
}

// Prepare nextflow input channels for process 2 & 3

mergelist_file = mergelist_shred
  .collectFile() { plink_prefix -> 
    [ "gwas_mergelist.txt", plink_prefix + '\n' ] 
  }

p1_out_pgen
  .into { p2_in_pgen; p3_in_pgen }
p1_out_pvar
  .into { p2_in_pvar; p3_in_pvar }
p1_out_psam
  .into { p2_in_psam; p3_in_psam }

p3_in_pgen
  .map { file ->
    tuple(file.getName().replaceFirst(/.pgen/, ""), file)
  }
  .set { p3_in_pgen_new }
p3_in_pvar
  .map { file ->
    tuple(file.getName().replaceFirst(/.pvar/, ""), file)
  }
  .set { p3_in_pvar_new }
p3_in_psam
  .map { file ->
    tuple(file.getName().replaceFirst(/.psam/, ""), file)
  }
  .set { p3_in_psam_new }

p3_in_pgen_new
  .join(p3_in_pvar_new)
  .join(p3_in_psam_new)
  .set { p3_in_files }

p3_in_files
  .into { p3_in_files_gallop; p3_in_files_plink }

// Merge plink files from process 1 for process 2
plink_prefix = "All_chr_${params.dataset}_p2in"

process p2_merge_list {
  scratch true

  label 'small'

  input:
    file mergelist from mergelist_file
    file "chr*_${params.dataset}_p1out.pgen" from p2_in_pgen.collect()
    file "chr*_${params.dataset}_p1out.pvar" from p2_in_pvar.collect()
    file "chr*_${params.dataset}_p1out.psam" from p2_in_psam.collect()
  output:
    tuple file("${plink_prefix}.bed"), file("${plink_prefix}.fam"), file("${plink_prefix}.bim") into p2_in_plink1
    tuple file("${plink_prefix}.pgen"), file("${plink_prefix}.pvar"), file("${plink_prefix}.psam") into p2_in_plink2
    file "*.log" into p2_mergelist_log

  script:
    """
    set -x
    ls *.pgen | sed 's/.pgen//g' > tmp_mergefile.txt
    plink2 --pmerge-list "tmp_mergefile.txt" --make-bed --out "${plink_prefix}"
    """
}

process p2_qc_pipeline {
  scratch true

  label 'medium'
  storeDir "${GWAS_STORE_DIR}/p2_qc_pipeline_cache"

  publishDir "${GWAS_OUTPUT_DIR}/${params.out}_${params.datetime}/Plots", mode: 'copy', overwrite: true

  input:
    set file("${plink_prefix}.bed"), file("${plink_prefix}.fam"), file("${plink_prefix}.bim") from p2_in_plink1
    set file("${plink_prefix}.pgen"), file("${plink_prefix}.pvar"), file("${plink_prefix}.psam") from p2_in_plink2
    //path ref_labels, stageAs: "ancestry_ref_labels.txt" from local_ref_labels
    //file "*" from ref_panel_files.collect()
  output:
    file "${params.ancestry}_samplelist_p2out.h5" into p2_out_file
    file "*.html" into p2_out_html
    file "*.png" into p2_out_png

  script:
    """
    set -x
    python3 ${ADDI_QC_PIPELINE} \
      --geno "${plink_prefix}" \
      --ref "/srv/GWAS-Pipeline/References/ref_panel/1kg_ashkj_ref_panel_gp2_pruned_hg38_newids" \
      --ref_labels "/srv/GWAS-Pipeline/References/ref_panel/ancestry_ref_labels.txt" \
      --pop "${params.ancestry}" \
      --out "${params.ancestry}_samplelist_p2out"
    """
}


process p3_cohort_samplelist {
  scratch true
  label 'small'

  input:
    file samplelist from p2_out_file
    path covarfile, stageAs: 'covariates.tsv' from "${params.covarfile}"
  output:
    file "${params.ancestry}_*_filtered.tsv" into gwas_samplelist

  script:
    //println "Separating IID based on ancestry (${ancestry})"
    """
    #!/usr/bin/env python3
    import pandas as pd
    import time

    ancestry = "${params.ancestry}"
    study_id_colname = "${params.study_id}"

    ancestry_df = pd.read_hdf("${samplelist}", key="ancestry_keep")
    outlier_df = pd.read_hdf("${samplelist}", key="outliers")
    pcs_df = pd.read_hdf("${samplelist}", key="pcs")
    kin_df = pd.read_hdf("${samplelist}", key="kin")
    data_df = pd.read_csv('covariates.tsv', sep="\\t", engine='c')

    cohorts = data_df[study_id_colname].unique().tolist()

    kin_df = kin_df[kin_df.KINSHIP >= ${params.kinship}]
    # TODO: address case when single cohort present
    # TODO: currently does not address longitudinal covariates
    for cohort in cohorts:
      print(f'---- {cohort} ----')
      samples = data_df[ (data_df.IID.isin(ancestry_df.IID)) &
                         (data_df[study_id_colname] == cohort) &
                        ~(data_df.IID.isin(outlier_df.IID)) ].copy(deep=True)
      samples = samples.merge(pcs_df, left_on="IID", right_on="#IID", how="inner")
      samples.drop(columns="#IID", inplace=True)
      r = kin_df[(kin_df['#IID1'].isin(samples.IID)) & (kin_df.IID2.isin(samples.IID))].copy()
      samples = samples[~samples.IID.isin(r.IID2)].copy()
      samples.to_csv(f"{ancestry}_{cohort}_filtered.tsv", sep="\t", index=False)
      print(f'Samples removed (outliers) = {data_df.IID.isin(outlier_df.IID).sum()}')
      print(f'Samples removed (kinship) = {r.shape[0]}')
      print(f'Samples remaining = {len(samples)}')
      print('')

    time.sleep(10)
    """
}

gwas_samplelist
  .into { gwas_samplelist_plink; gwas_samplelist_gallop }


process p3_export_rawfile {
  scratch true
  label 'small'

  input:
    set val(plink_prefix), file(psam), file(pgen), file(pvar) from p3_in_files_gallop
    each file(samplelist) from gwas_samplelist_gallop.flatten()
  output:
    tuple file(samplelist), file('*.raw') into gwas_rawfile 
  when:
    params.longitudinal_flag

  script:
    def m = []
    def cohort = ""
    cohort = samplelist.getName()
    m = cohort =~ /(.*)_filtered.tsv/
    cohort = m[0][1]

    m = []
    def chrnum = ""

    // case insensitive matching to get the chromosome number from input parameter or filename
    if (params.chr == "")
      m = plink_prefix  =~ /(?i)(chr)([0-9]+)/
    else
      m = params.chr =~ /(?i)(chr)([0-9]+)/

    // if matching fails - output error for user
    if (m.size() != 1 || m[0].size() < 3)
      error "Failed to identify chromosome from filename ${dataFile}, please rename or use the --chr argument"
    
    chrnum = "chr${m[0][2]}"
    def outfile = "${cohort}_${chrnum}"

    """
    set -x
    plink2 --pfile ${plink_prefix} \
           --keep ${samplelist} \
           --export A \
           --mac 20 \
           --update-sex ${samplelist} \
           --pheno ${samplelist} \
           --pheno-col-nums 4 \
           --hwe 1e-6 \
           --geno ${params.missing_geno_rate} \
           --out ${outfile}  \
           --threads ${task.cpus}
    """
}

process p3_gwas_gallop {
  scratch true
  label 'small'

  publishDir "${GWAS_OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    set file(samplelist), file(rawfile) from gwas_rawfile
    path x, stageAs: 'phenotypes.tsv' from "${params.phenofile}"
  output:
    file "*.gallop" into gallop_results
  when:
    params.longitudinal_flag
  
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
    gallop --gallop \
           --rawfile ${rawfile} \
           --pheno "phenotypes.tsv" \
           --covar ${samplelist} \
           --covar-name ${params.covariates} \
           ${model} ${pheno_name} ${pheno_name_file} \
           --time-name 'study_days' \
           --out "${outfile}"
    """
}


process p3_format_gwas_plink {
  echo true
  scratch true
  label 'small'

  input:
    val samplelist from gwas_samplelist_plink.flatten()
    path x, stageAs: 'phenotypes.tsv' from "${params.phenofile}"
  output:
    file "*_analyzed.tsv" into plink_samplelist
  when:
    ! params.longitudinal_flag

  script:
    //println "Samplelist: ${samplelist}"
    //println "Phenotypes: phenotypes.tsv"

    def m = []
    def cohort = ""
    cohort = samplelist.getName()
    m = cohort =~ /(.*)_filtered.tsv/
    outfile = "${m[0][1]}"

    def pheno_name = "y"
    if (params.pheno_name != '') {
      pheno_name = "${params.pheno_name}"
    }

    """
    #!/usr/bin/env python3
    import pandas as pd
    import time

    covars = "${params.covariates}"
    covars = covars.split(' ')

    d_pheno = pd.read_csv("phenotypes.tsv", sep="\t", engine='c')
    d_sample = pd.read_csv("${samplelist}", sep="\t", engine='c')

    d_result = pd.merge(d_pheno, d_sample, on='IID', how='inner')
    d_set = d_result.loc[:, ["#FID", "IID", "${pheno_name}"] + covars].copy()
    d_set.to_csv("${outfile}_analyzed.tsv", sep="\t", index=False)

    time.sleep(10)
    """
}


process p3_gwas_plink{
  scratch true
  label 'small'
  publishDir "${GWAS_OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    set val(plink_prefix), file(psam), file(pgen), file(pvar) from p3_in_files_plink
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

    m = []
    def chrnum = ""

    // case insensitive matching to get the chromosome number from input parameter or filename
    if (params.chr == "")
      m = plink_prefix  =~ /(?i)(chr)([0-9]+)/
    else
      m = params.chr =~ /(?i)(chr)([0-9]+)/

    // if matching fails - output error for user
    if (m.size() != 1 || m[0].size() < 3)
      error "Failed to identify chromosome from filename ${dataFile}, please rename or use the --chr argument"
    
    chrnum = "chr${m[0][2]}"
    def outfile = "${cohort}_${chrnum}.${params.out}"

    def pheno_name = "y"
    if (params.pheno_name != '') {
      pheno_name = "${params.pheno_name}"
    }

    """
    plink2 --pfile ${plink_prefix} \
           --glm hide-covar omit-ref cols=+beta,+a1freq \
           --pheno ${samplelist} \
           --pheno-name ${pheno_name} \
           --covar ${samplelist} \
           --covar-name ${covariates} \
           --covar-variance-standardize \
           --keep ${samplelist} \
           --output-chr chrM \
           --mac 20 \
           --hwe 1e-6 \
           --threads ${task.cpus} \
           --out ${outfile}
    """
}

process p3_gwas_viz {
  scratch true

  publishDir "${GWAS_OUTPUT_DIR}/${params.out}_${params.datetime}/Plots", mode: 'copy', overwrite: true

  input:
    file x from gwas_results.mix(gallop_results).collect()
  output:
    file "*.png" into gwas_plots

  script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    from qmplot import manhattanplot
    from qmplot import qqplot

    files = list(filter(lambda x: os.path.splitext(x)[-1] in ['.gallop', '.linear'], os.listdir()))
    dfs = []
  
    # group by phenotype
    lt_flag = "${params.longitudinal_flag}" == "true"
    suffix = "${params.out}"

    def plot_summary_stats(data, cohort, outcome, lt_flag=False):
      xtick = set(['chr' + i for i in list(map(str, range(1, 10))) + ['11', '13', '15', '18', '21', 'X']])
      if lt_flag:
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="Pi",
                      figname=f"{cohort}_{outcome}_manhattan_intercept.gallop.png",
                      xtick_label_set=xtick)
        
        manhattanplot(data=data,
                      title=f"Manhattan Slope {cohort} {outcome}",
                      pv="P",
                      figname=f"{cohort}_{outcome}_manhattan_slope.gallop.png",
                      xtick_label_set=xtick)
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["Pi"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               dpi=300,
               ax=ax,
               figname=f"{cohort}_{outcome}_qq_intercept.gallop.png")
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ Slope {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               dpi=300,
               ax=ax,
               figname=f"{cohort}_{outcome}_qq_slope.gallop.png")
      else:
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="P",
                      figname=f"{cohort}_{outcome}_manhattan.linear.png",
                      xtick_label_set=xtick)
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               dpi=300,
               ax=ax,
               figname=f"{cohort}_{outcome}_qq.linear.png")

    plot_df = dict()
    for f in files:
      fc = f.split('.')
      pheno = fc[-2]
      cohort = '_'.join(fc[0].split('_')[:-1])
      try:
        plot_df[f'{cohort}+{pheno}'].append(f)
      except KeyError:
        plot_df[f'{cohort}+{pheno}'] = [f]


    for cohort in plot_df.keys():
      dfs = []
      c, outcome = cohort.split('+')
      for f in plot_df[cohort]:
        dfs.append(pd.read_table(f, sep="\t"))
      
      df = pd.concat(dfs)
      df = df.dropna(how="any", axis=0)  # clean data
      df['chr_order'] = df['#CHROM'].str.replace('chr', '')
      df['chr_order'] = df['chr_order'].astype(int)
      df = df.sort_values(by=['chr_order', 'POS'])

      # generate manhattan plot and set an output file.
      plot_summary_stats(data=df, cohort=f'{c}.{suffix}', outcome=outcome, lt_flag=lt_flag)
    """
}

