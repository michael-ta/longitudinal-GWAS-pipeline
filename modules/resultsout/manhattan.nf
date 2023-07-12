process MANHATTAN {
  scratch true
  label 'medium'

  //publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/", mode: 'copy', overwrite: true
  publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/${model}_MANHATTAN_${params.datetime}", mode: 'copy', overwrite: true

  input:
    path x //from gwas_results.mix(sumstats_out).collect()
    val(model)
  output:
    path "*.png" //into gwas_plots

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
    
    if "${params.survival_flag}" == "true":
      gwasname = "CPH"
    else:
      gwasname = "GLM"

    suffix = "${params.out}"
    threads = int("${task.cpus}")

    def plot_summary_stats(data, cohort, outcome, lt_flag=False):
      xtick = set(['chr' + i for i in list(map(str, range(1, 10))) + ['11', '13', '15', '18', '21', '22']])
      #xtick = set(['chr' + i for i in list(map(str, range(1, 10))) + ['11', '13', '15', '18', '21', 'X']])
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
        plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.{gwasname}.png", dpi=300)
        f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_intercept.{gwasname}.png", dpi=300)

    plot_df = dict()
    
    for f in files:
      fc = f.split('.')
      pheno = fc[-2]
      cohort = fc[0][0:fc[0].find('chr')]
      if cohort[-1] == '_':
        cohort = cohort[:-1]
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