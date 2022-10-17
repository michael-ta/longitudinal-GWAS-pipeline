library('survival')
library('optparse')

option_list <- list(
  make_option( c("--covar"), type="character", default=NULL,
                 help="cross-sectional covariates"),
  make_option( c("--pheno"), type="character", default=NULL,
                 help="phenotype or outcome"),
  make_option( c("--rawfile"), type="character", default=NULL,
                 help="rawfile"),
  make_option( c("--covar-name"), type="character", default=NULL,
                 help="space delimited covariate list"),
  make_option( c("--pheno-name"), type="character", default="y",
                 help="phenotype / outcome column name"),
  make_option( c("--out"), type="character", 
                 default="survival.tbl",
                 help="output file")
) 

# the following are requirements for survival analysis
# tstart - time start for period, always 0 for cross-sectional analysis
# tend - time end for period


parser <- OptionParser(option_list=option_list)
arguments <- parse_args( parser, positional_arguments=TRUE )
opt <- arguments$options
args <- arguments$args

input.covariates <- strsplit(opt[['covar-name']], ' ')

print(input.covariates)


data.pheno = read.table(opt$pheno, header=TRUE)
data.covar = read.table(opt$covar, header=TRUE)
data.geno = read.table(opt$rawfile, header=TRUE)

data.merged = merge(data.covar, data.pheno)

# ----- Create SumStats Table ----
SNPs = colnames(data.geno)[grepl("^chr[0-9]", names(data.geno))]
print( paste("Found", as.character(as.numeric(length(SNPs))), "SNPs") )

tmp.data <- data.frame(
  sapply(SNPs,
         function(x) strsplit(x, '\\.')))

tmp.data = as.data.frame(t(tmp.data))
colnames(tmp.data) <- c("#CHROM", "POS", "REF", "REF_ALT")
tmp.data$ALT = do.call(rbind,
  sapply(tmp.data$REF_ALT,
         function(x) strsplit(x, '_')[1]))[,1]
tmp.data$ID = paste(tmp.data[["#CHROM"]],
                    tmp.data$POS,
                    tmp.data$REF,
                    tmp.data$ALT,
                    sep=":")
tmp.data$A1 = tmp.data$ALT
tmp.data$A1_FREQ = as.numeric(
  lapply(data.geno[SNPs], function(x) sum(x))) / as.numeric(
    colSums(!is.na(data.geno[SNPs])) * 2)
tmp.data$MISS_FREQ = as.numeric(colSums(is.na(data.geno[SNPs]))) / as.numeric(
  lapply(data.geno[SNPs], function(x) length(x)))
tmp.data$TEST = 'CoxPH'

# ------ Fit CoxPH model ----

## prepare results matrix
stats <- matrix(NA, nrow(tmp.data), 4)
idx = 1

basemod <- paste0("Surv(tstart,tend,", opt[['pheno-name']], ")~")
basemod <- paste0(basemod, paste(unlist(input.covariates), collapse="+"))
mod_cols = c('coef', 'exp(coef)', 'se(coef)', 'Pr(>|z|)')
print( paste("Base survival model", basemod) )


for (i in colnames(data.geno)[grepl( "^chr[0-9]", names(data.geno))]) {
    data.mtx = merge( data.merged, data.geno[, c("IID", i)], by='IID' )
    eq = paste0(basemod, "+", i)
    mdl = coxph(as.formula(eq), data=data.mtx)
    res = summary(mdl)
    stats[idx,] = res$coefficients[i,][mod_cols]
    idx = idx + 1
}

stats = as.data.frame(stats)
colnames(stats) = c('BETA', 'exp(BETA)', 'SE', 'P')

stats = cbind(tmp.data, stats)
write.csv(stats, opt$out, row.names=FALSE)

