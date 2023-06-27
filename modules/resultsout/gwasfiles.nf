process SAVEGWAS {
  scratch true
  publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}", mode: 'copy', overwrite: true
  
  input:
    path(sumstats) 
    
  output:
    path(sumstats) //into sumstats_out
    
  script:
    """
    bedtools sort -i ${sumstats} -header > "${sumstats}".tmp
    mv ${sumstats}.tmp ${sumstats}
    """
}