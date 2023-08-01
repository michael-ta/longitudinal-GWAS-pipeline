process SAVEGWAS {
  scratch true
  publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/${model}_${params.datetime}", mode: 'copy', overwrite: true

  input:
    tuple val(pheno), path(sumstats)
    val(model)
  output:
    path(sumstats) , emit: res_split
    path "${pheno}_allresults.tsv", emit: res_all

  script:
    """
    echo $pheno
    COUNTER=0
    for file in $sumstats
    do
      COUNTER=\$((COUNTER+1))
      if [[ \$COUNTER -eq 1 ]]
      then
        cat \${file} > "${pheno}_allresults.tsv"
      else
        tail -n +2 \${file} >> "${pheno}_allresults.tsv"
      fi
    done
    bedtools sort -i "${pheno}_allresults.tsv" -header > "${pheno}_allresults.tsv.tmp"
    mv "${pheno}_allresults.tsv.tmp" "${pheno}_allresults.tsv"
    """
}
// process SAVEGWAS {
//   scratch true
//   publishDir "${OUTPUT_DIR}/${params.dataset}/RESULTS/${model}_${params.datetime}", mode: 'copy', overwrite: true
  
//   input:
//     path(sumstats)
//     val(model)
    
//   output:
//     path(sumstats) //into sumstats_out
    
//   script:
//     """
//     bedtools sort -i ${sumstats} -header > "${sumstats}".tmp
//     mv ${sumstats}.tmp ${sumstats}
//     """
// }
