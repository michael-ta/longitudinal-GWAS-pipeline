process SAVEGWAS {
  scratch true
  publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/${model}_${params.datetime}", mode: 'copy', overwrite: true
  
  input:
    path(sumstats)
    val(model)
    
  output:
    path(sumstats) //into sumstats_out
    
  script:
    """
    bedtools sort -i ${sumstats} -header > "${sumstats}".tmp
    mv ${sumstats}.tmp ${sumstats}
    """
}