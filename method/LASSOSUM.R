library(lassosum)
setwd(system.file("data", package="lassosum")) # Directory where data and LD region files are stored
library(data.table)
temp <- commandArgs(TRUE)
chr =  as.numeric(temp[1])
trait = temp[2]
ref_N =  temp[3]

sum.raw = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/GWAS_ukbb/', trait, '/sumdat_Rcov.txt'))
sum.raw = as.data.frame(sum.raw)
sumdat = sum.raw[sum.raw$CHR == chr, ]

if (ref_N == "5000") { ref.bfile <- paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/chr", chr, "_5000") }
if (ref_N == "1kg") { ref.bfile <- paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg/chr", chr, "/chr_", chr, "_merged") }

test.bfile <- paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation/chr", chr)
LDblocks <- "EUR.hg19"

library(lassosum)
library(parallel)
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
system(paste0("mkdir -p ", '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, "/"))
bigreadr::fwrite2(beta, paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'), sep = '\t')


prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-83 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, "/ref_N", ref_N, "validation"))
system(prscode)

prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, '/LASSOSUMeffect-hm3-EUR-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-83 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM/chr', chr, "/ref_N", ref_N, "tuning"))
system(prscode)
