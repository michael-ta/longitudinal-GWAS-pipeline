process MERGER_SPLITS {
  scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p1_run_cache"
  publishDir "${OUTPUT_DIR}/${params.dataset}/LOGS/MERGER_SPLITS_${params.datetime}/", mode: 'copy', overwrite: true, pattern: "*.log"
  label 'small'

  input:
    file mergelist
    file "*" //from p1_run_processed_ch_chunk.collect()
  output:
    tuple file("${vSimple}.psam"), file("${vSimple}.pgen"), file("${vSimple}.pvar"), file("${vSimple}.log"), emit: snpchunks_qc_merged //into p1_merge_chunk_processed_ch
    //file "*.log"
  script:
    vSimple = mergelist.getSimpleName()
    if (params.chunk_flag) {
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
    } else {
      """
      # run merge command on tmp list
      set +x
      
      plink2 -pfile ${vSimple}.1_p1out \
        --make-pgen \
        --sort-vars \
        --out ${vSimple}
      """
    }
}

process MERGER_CHRS {
  scratch true
  storeDir "${STORE_DIR}/${params.dataset}/p2_merged_cache"
  //publishDir "${OUTPUT_DIR}/${params.out}_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  publishDir "${OUTPUT_DIR}/${params.dataset}/MERGER_CHRS_LOGS_${params.datetime}/logs", mode: 'copy', overwrite: true, pattern: "*.log"
  
  label 'large_mem'

  input:
    file mergelist
    file "*" //from input_p2_merge_list_files.collect()
  output:
    //file (" ${plink_prefix}.{bed,fam,bim,pgen,pvar,psam}") //into input_p2_merged_plink
    file ("allchr_${params.dataset}_p2in.{bed,fam,bim,pgen,pvar,psam,log}") //into input_p2_merged_plink
    //file "*.log" into p2_mergelist_log
  script:
    """
    set -x
    cat $mergelist | uniq > tmp_mergefile.txt
    plink2 --memory ${task.memory.toMega()} \
      --pmerge-list "tmp_mergefile.txt" \
      --make-bed \
      --out "allchr_${params.dataset}_p2in"
    """
}