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

process GENETICQC {
  scratch true
  label 'medium'
  //storeDir "${STORE_DIR}/${params.dataset}/p1_run_cache"

  input:
    tuple val(vSimple), path(fOrig), path(fSplit) //from input_p1_run_ch
  output:
    tuple file("${output}.pgen"), file("${output}.pvar"), file("${output}.psam"), emit: snpchunks_qc //into p1_run_processed_ch
    tuple val(vSimple), val("${output}"), emit: snpchunks_names //into p1_run_chunklist_ch

  script:
  def chrnum = ""
  def vPart = ""
  def prefix = ""

  // regex match the chromosome corresponding to the input file chunk from the simple filename - requires that
  // each input file have chr0-9 in the filename; also require the user not to include any periods in the
  // filename
  def mChrom = vSimple =~ /(?i)(chr)([0-9]+)/
  chrnum = mChrom[0][2]
  vPart = fSplit.getBaseName()

  // basename should contain the filename without the compression (.gz) extension - requires the input file to
  // have some sort of .extension in the filename. We can regex match split part from the remaining filename
  def mPart = vPart =~ /(.*)\.([0-9]+)\.(.*)$/
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