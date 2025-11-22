## ldpred2 real data analysis
rm(list=ls())
setwd('/fastscratch/myscratch/ydun/')
library(bigsnpr)

temp <- commandArgs(TRUE)
race = "EUR"
chr =  as.numeric(temp[1])
ref_N = "blk_provided"
trait = temp[2]

sum.raw = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/GWAS_ukbb/', trait, '/sumdat_Rcov.txt'))
sum.raw = as.data.frame(sum.raw)

map_ldref <- readRDS("/dcs04/nilanjan/data/ydun/ldpred2/ldref/map.rds")
sumstats = sum.raw[sum.raw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
info_snp <- snp_match(sumstats, map_ldref)
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
## indices in G
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
corr_chr <- readRDS(paste0("/dcs04/nilanjan/data/ydun/ldpred2/ldref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
corr <- as_SFBM(corr_chr, compact = TRUE)


beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = 4)
beta_grid = as.data.frame(beta_lassosum2)
beta_grid[is.na(beta_grid)] = 0
rownames(beta_grid) = info_snp$rsid
beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:120))

system(paste0("mkdir -p ", '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/"))
bigreadr::fwrite2(beta_grid, paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'), sep = '\t')

prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_N", ref_N, "validation"))
system(prscode)

prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-', race, '-ref_N', ref_N, '.txt'),
                paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr),
                paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_N", ref_N, "tuning"))
system(prscode)

