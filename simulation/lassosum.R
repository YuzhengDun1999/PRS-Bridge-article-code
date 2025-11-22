library(lassosum)
setwd(system.file("data", package="lassosum")) # Directory where data and LD region files are stored
library(data.table)
temp <- commandArgs(TRUE)

chr =  as.numeric(temp[1])
rho = as.numeric(temp[2])
#size = as.numeric(temp[4])
GA = as.numeric(temp[3])
#rep = as.numeric(temp[6])
size = 4
rep = 1

sum.raw = bigreadr::fread2(paste0('/dcl01/chatterj/data/jin/prs/simulation/EUR/sumdata/megasum-rho',rho,'-size',size,'-rep',rep,'-GA',GA,'-chr',chr,'.txt'))
sum.raw = as.data.frame(sum.raw)
PRScs = bigreadr::fread2(paste0("/dcs04/nilanjan/data/ydun/PRScs/simulation_result/chr", chr, "-rho3_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt"))
sum.raw = sum.raw %>% filter(POS %in% PRScs[,3])
library(dplyr)
sumdat = sum.raw %>% filter(((BETA/SE)^2/N<=80)) %>% filter((REF_FRQ>0.01) & (REF_FRQ<0.99)) %>% filter(!(CHR==6 & POS>26e6 & POS<34e6))
ref.bfile <- paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg/chr", chr, "/chr_", chr, "_merged")

test.bfile <- paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/tuning/chr", chr)
LDblocks <- "EUR.hg19"

library(lassosum)
library(parallel)
cl <- makeCluster(2, type="FORK")
cor <- p2cor(p = as.numeric(sumdat$PVAL), n = as.numeric(sumdat$N), sign=sumdat$BETA)
cor[is.na(cor)] = max(abs(cor), na.rm = TRUE) * sign(sumdat$BETA[is.na(cor)])
out <- lassosum.pipeline(cor=cor, chr=sumdat$CHR, pos=sumdat$POS, A1=sumdat$ALT, A2=sumdat$REF,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, LDblocks = LDblocks, cluster=cl)
pgs = data.frame(IID = 100001:110000)
for (i in 1:length(out$pgs)) {
  temp = out$pgs[[i]]; temp[is.na(temp)] = 0
  pgs = cbind(pgs, temp)
}
names(pgs) = c("IID", paste0("e", 1:80))
system(paste0("mkdir -p ", '/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/LASSOSUM/rho', rho, '-size', size, '-rep', rep, '-GA',GA, '/chr', chr, "/"))
bigreadr::fwrite2(pgs, paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/LASSOSUM/rho', rho, '-size', size, '-rep', rep, '-GA',GA, '/chr', chr, "/tuning_prs.txt"), sep = '\t')
