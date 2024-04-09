process MANHATTAN {
  scratch true
  label 'medium'

  publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/${model}_MANHATTAN_${params.datetime}", mode: 'copy', overwrite: true

  input:
    each x
    val(model)
  output:
    path "*.png"

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import matplotlib.pyplot as plt
  import os
  from qmplot import manhattanplot
  from qmplot import qqplot

  def plot_summary_stats(data, cohort, outcome, model):
    xtick = set(['chr' + i for i in list(map(str, range(1, 14))) + ['15', '17', '19', '22']])

    if model == "lmm_gallop":
      f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
      manhattanplot(data=data,
                    title=f"Manhattan Intercept {cohort} {outcome}",
                    pv="Pi", ax = ax, 
                    xtick_label_set=xtick)
      plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.{model}.png", dpi=300)
      f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
      qqplot(data=data["Pi"],
              marker="o",
              title=f"QQ Intercept {cohort} {outcome}",
              xlabel=r"Expected -log(P)",
              ylabel=r"Observed -log(P)",
              ax=ax)
      plt.savefig(f"{cohort}_{outcome}_qq_intercept.{model}.png", dpi=300)
      
      f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
      manhattanplot(data=data,
                    title=f"Manhattan Slope {cohort} {outcome}",
                    pv="P", ax = ax, 
                    xtick_label_set=xtick)
      plt.savefig(f"{cohort}_{outcome}_manhattan_slope.{model}.png", dpi=300)
      f, ax = plt.subplots(figsize=(10, 7), facecolor="w", edgecolor="k")
      qqplot(data=data["Pi"],
              marker="o",
              title=f"QQ Slope {cohort} {outcome}",
              xlabel=r"Expected -log(P)",
              ylabel=r"Observed -log(P)",
              ax=ax)
      plt.savefig(f"{cohort}_{outcome}_qq_slope.{model}.png", dpi=300)
    
    elif model in ["glm","cph"]:
      f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
      manhattanplot(data=data,
                    title=f"Manhattan {model} {cohort} {outcome}",
                    pv="P", ax = ax,
                    xtick_label_set=xtick)
      plt.savefig(f"{cohort}_{outcome}_manhattan.{model}.png", dpi=300)
      f, ax = plt.subplots(figsize=(10, 7), facecolor="w", edgecolor="k")
      qqplot(data=data["P"],
              marker="o",
              title=f"QQ {model} {cohort} {outcome}",
              xlabel=r"Expected -log(P)",
              ylabel=r"Observed -log(P)",
              ax=ax)
      plt.savefig(f"{cohort}_{outcome}_qq.{model}.png", dpi=300)
    
    else:
      print("Model flag not recognised")
      return
    
    print("Results plotting success") 
  
  # Get some metadata
  threads = int("${task.cpus}")
  suffix = "${params.out}"
  model = "${model}"
  
  data_path = "${x}"
  gwas_name = os.path.splitext(os.path.basename(data_path))[0].split('_')
  pheno = gwas_name[-2]
  cohort = gwas_name[1]
  cohort_suffix = f"{cohort}.{suffix}" if suffix != "" else cohort


  # Prep data for plotting
  df = pd.read_csv("${x}", sep="\\t", engine='c')
  df = df.dropna(how="any", axis=0)  # clean data
  df['chr_order'] = df['#CHROM'].str.replace('chr', '')
  df['chr_order'] = df['chr_order'].astype(int)
  df = df.sort_values(by=['chr_order', 'POS'])
  
  # generate manhattan plot and set an output file.
  plot_summary_stats(data=df, cohort=f'{cohort_suffix}', outcome=pheno, model=model)
  
  """
}