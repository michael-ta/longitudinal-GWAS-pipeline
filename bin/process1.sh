#!/bin/bash
# Take the vcf file (genotyped or imputed)
# Return: pass(&R2)-filtered, split, left-normalized, autosomal-par, hg38-ref-alt-aligned SNPs with mac >=2. geno < 0.05

# Parameters
N=$1 # Threads to use
VFILE=$2
R2THRES=$3 # If imputed, give a number for R2 threshold (usually 0.3 - 0.8) -9 otherwise
ASSEMBLY=$4 # [hg18, hg19, hg38]. Define if the liftover is required or not
CHRNUM=$5 # [1..22] Needed for lift over. 
FILE=$6 # Base file name. can be anything as long as unique
## e.g. process1.sh 2 '/data/CARD/PD/imputed_data/CORIELL/chr21.dose.vcf.gz' 0.3 hg19 21 chr21_cor

# Resources (Make sure the paths are correct)
# hg38 FASTQ reference will be passed in by nextflow
FA=/srv/GWAS-Pipeline/References/Genome/hg38.fa.gz
LIFTOVERDATA=/srv/GWAS-Pipeline/References/liftOver/${ASSEMBLY}ToHg38.over.chain.gz


######## start processing ###############################
# Filter PASS (&R2) > Split > Create Plink binary (only keep autosome + par)
## Different pipeline for imputed and genotyped (R2THRES>0....imputed)
if [[ $R2THRES > 0 ]]
then 
    bcftools view -f '.,PASS' \
                  -i "INFO/R2>${R2THRES}" ${VFILE} \
                  -Oz -o ${FILE}.vcf.gz --threads ${N} # get PASS(or .) variant and R2>R2THRES

    bcftools norm -m-both ${FILE}.vcf.gz \
                  -Oz -o ${FILE}_split.vcf.gz --threads ${N} # split

    plink2 --threads ${N} \
           --vcf ${FILE}_split.vcf.gz dosage=DS \
           --make-pgen --allow-extra-chr --autosome-par --out ${FILE}_split # Import dosage
else
    bcftools view -f '.,PASS' ${VFILE} \
                  -Oz -o ${FILE}.vcf.gz --threads ${N} 
    bcftools norm -m-both ${FILE}.vcf.gz \
                  -Oz -o ${FILE}_split.vcf.gz --threads ${N} # split
    plink2 --threads ${N} \
           --vcf ${FILE}_split.vcf.gz --make-pgen --allow-extra-chr --autosome-par --out ${FILE}_split
fi

# left-normalize using fasta file (hg38)
if [[ $ASSEMBLY = hg38 ]]
then
    plink2 --threads ${N} \
           --pfile ${FILE}_split --make-pgen --fa $FA --normalize --sort-vars --out ${FILE}_split_hg38_normalized
else
    # need to liftover to hg38 before normalization
    ## giving temporary variant IDs (unique vID required for updating the map)
    plink2 --threads ${N} \
           --pfile ${FILE}_split --make-pgen \
           --set-all-var-ids "chr@:#[${ASSEMBLY}]:\$r:\$a" \
           --new-id-max-allele-len 999 truncate --out ${FILE}_split_temp_renamed
    ## remove dup
    plink2 --threads ${N} \
           --pfile ${FILE}_split_temp_renamed --make-pgen --rm-dup exclude-all --out ${FILE}_split_temp_renamed_uniq
    ## prepare input for liftover
    grep -v '#' ${FILE}_split_temp_renamed_uniq.pvar | awk '{print "chr"$1,$2,$2+1,$3}' > ${FILE}_split.liftoverInput 
    ### liftover
    liftOver ${FILE}_split.liftoverInput ${LIFTOVERDATA} ${FILE}_split.liftoverOutput ${FILE}_split.unMapped
    ### only keep variants in the same chromosome
    awk '/chr'${CHRNUM}'\t/{print}'  ${FILE}_split.liftoverOutput > ${FILE}_split.liftoverOutput.chr${CHRNUM}
    ### variant ID list to keep
    cut -f4 ${FILE}_split.liftoverOutput.chr${CHRNUM} > ${FILE}_split.liftoverOutput.chr${CHRNUM}.ID
    ## update the map
    plink2 --threads ${N} --pfile ${FILE}_split_temp_renamed_uniq --make-pgen --update-map ${FILE}_split.liftoverOutput.chr${CHRNUM} 2 4 --extract ${FILE}_split.liftoverOutput.chr${CHRNUM}.ID --sort-vars --out ${FILE}_split_hg38 
    ## noramlize with fa
    plink2 --threads ${N} --pfile ${FILE}_split_hg38 --make-pgen --fa $FA --normalize --sort-vars --out ${FILE}_split_hg38_normalized
fi

# select relevant snps with more than sigleton (For small cohorts, this process reduces a lot of variants)
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized --make-pgen --snps-only just-acgt --mac 2 --out ${FILE}_split_hg38_normalized_snps
# align ref alt (This will return provisional variants sometime)
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps --make-pgen --ref-from-fa force --fa $FA --out ${FILE}_split_hg38_normalized_snps_temp_aligned
# remove the "Provisional variants" because provisional variants prevents loading on plink2
grep 'PR$' ${FILE}_split_hg38_normalized_snps_temp_aligned.pvar | cut -f3 > ${FILE}_split_hg38_normalized_snps_temp_aligned_provisional.txt
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps --make-pgen --ref-from-fa force --fa $FA --exclude ${FILE}_split_hg38_normalized_snps_temp_aligned_provisional.txt --out ${FILE}_split_hg38_normalized_snps_aligned
# remame ID for standard chr:pos:ref:alt
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned --make-pgen --set-all-var-ids 'chr@:#:$r:$a' --out ${FILE}_split_hg38_normalized_snps_aligned_renamed
# remove dup
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned_renamed --make-pgen --rm-dup exclude-all --out ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq
# geno 0.05
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq --make-pgen --geno 0.05 dosage --out ${FILE}_p1out
