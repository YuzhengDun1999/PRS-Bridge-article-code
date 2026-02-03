########### This script is used to run LASSOSUM #############
library(lassosum)
library(parallel)
setwd(system.file("data", package="lassosum")) # Directory where data and LD region files are stored
library(data.table)

trait = "Example"
ref_N = "1kg" # can be 1kg and individuals from UKBB data

for (chr in 1:22) {
  sum.raw = bigreadr::fread2('data/sumdat_Rcov.txt')
  sum.raw = as.data.frame(sum.raw)
  sumdat = sum.raw[sum.raw$CHR == chr, ]
  ref.bfile = paste0("chr", chr) # path to your LD reference file, has to be individual level data
  
  test.bfile <- paste0("chr", chr) # path to your test gentic data
  LDblocks <- "EUR.hg19"
  cl <- makeCluster(2, type="FORK")
  sumdat$PVAL[sumdat$PVAL < .Machine$double.xmin] = .Machine$double.xmin * 10
  cor <- p2cor(p = as.numeric(sumdat$PVAL), n = as.numeric(sumdat$N), sign=sumdat$BETA)
  cor[is.na(cor)] = max(abs(cor), na.rm = TRUE) * sign(sumdat$BETA[is.na(cor)])
  out <- lassosum.pipeline(cor=cor, chr=sumdat$CHR, snp=sumdat$SNP_ID, A1=sumdat$ALT, A2=sumdat$REF,
                           ref.bfile=ref.bfile, test.bfile=test.bfile, LDblocks = LDblocks, cluster=cl,
                           destandardize = TRUE)
  beta = out$sumstats[,c(2, 3, 4)]
  for (i in 1:length(out$beta)) {
    temp = out$beta[[i]]; temp[is.na(temp)] = 0
    beta = cbind(beta, temp)
  }
  beta[,1:3] = beta[, c(3, 2, 1)]
  names(beta) = c("snp", "A1", "A2", paste0("e", 1:80))
  system(paste0("mkdir -p ", trait, '/LASSOSUM/chr', chr, "/"))
  bigreadr::fwrite2(beta, paste0(trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'), sep = '\t')
  
  
  prscode = paste(paste0('plink2'),
                  paste0('--score ', trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-83 --bfile chr', chr),
                  paste0(' --out ', trait, '/LASSOSUM/chr', chr, "/ref_N", ref_N, "validation"))
  system(prscode)
  
  prscode = paste(paste0('plink2'),
                  paste0('--score ', trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-83 --bfile chr', chr),
                  paste0(' --out ', trait, '/LASSOSUM/chr', chr, "/ref_N", ref_N, "tuning"))
  system(prscode)
}
