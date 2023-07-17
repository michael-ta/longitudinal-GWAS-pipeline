process GWASQC {
  //scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p2_qc_pipeline_cache"
  publishDir "${OUTPUT_DIR}/${params.dataset}/PLOTS/GWASQC_PLOTS_${params.datetime}/", mode: 'copy', overwrite: true, pattern: "*.{html,png}"
  label 'large_mem'
  
  input:
    file "*" //from input_p2_qc_pipeline.collect()
  output:
    path "${params.ancestry}_samplelist_p2out.h5", emit: gwasqc_h5_file //into p2_qc_processed
    path "*.{html,png}", emit: gwasqc_figures //into p2_out_figs
  
  script:
    """
    set +x
    
    python3 ${ADDI_QC_PIPELINE} \
      --geno "allchr_${params.dataset}_p2in" \
      --ref "/srv/GWAS-Pipeline/References/ref_panel/1kg_ashkj_ref_panel_gp2_pruned_hg38_newids" \
      --ref_labels "/srv/GWAS-Pipeline/References/ref_panel/ancestry_ref_labels.txt" \
      --pop "${params.ancestry}" \
      --out "${params.ancestry}_samplelist_p2out"    
    """
}