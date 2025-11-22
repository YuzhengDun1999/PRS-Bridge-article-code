rm(list=ls())
setwd('/fastscratch/myscratch/ydun/')
library(bigsnpr)

temp <- commandArgs(TRUE)
race = "EUR"
chr =  as.numeric(temp[1])
ref_N = 5000
trait = temp[2]


ldr = 3/1000

sum.raw = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/GWAS/sumdat_Rcov.txt'))
sum.raw = as.data.frame(sum.raw)

# ------------------------ Run LASSOSUM2
#NCORES <- nb_cores()
# Read external summary statistics
temfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/chr',chr,'_', ref_N, '.bk')
if (!file.exists(temfile)){
  snp_readBed(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/chr',chr,'_', ref_N,'.bed'))
}
obj.bigSNP <- snp_attach(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/chr',chr,'_', ref_N,'.rds'))
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

sumstats = sum.raw[sum.raw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]

set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, strand_flip = T)
rownames(info_snp) = info_snp$rsid
## compute correlation
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/ydun/ldpred2/tmp_dir/EUR/'), ncores = 1)
## indices in info_snp
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
## indices in G
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 1, #size = ldradius)
                 infos.pos = POS2[ind.chr2], size = ldr) # default

system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, '/LASSOSUM2/chr', chr, "/tmp_ref", ref_N, "/"))
corr <- bigsparser::as_SFBM(corr0, tempfile(tmpdir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, '/LASSOSUM2/chr', chr, "/tmp_ref", ref_N, "/tmp-data")), compact = TRUE)
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = 1)
beta_grid = as.data.frame(beta_lassosum2)
beta_grid[is.na(beta_grid)] = 0
rownames(beta_grid) = info_snp$rsid
beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'), paste0('e',1:120))

system(paste0("mkdir -p ", '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/"))
bigreadr::fwrite2(beta_grid, paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'), sep = '\t')

prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_N", ref_N, "validation"))
system(prscode)

prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_N", ref_N, "tuning"))
system(prscode)
system(paste0("rm -r /dcs04/nilanjan/data/ydun/SinglePRS/", trait, '/LASSOSUM2/chr', chr, "/tmp_ref", ref_N, "/"))
