########### This script is used to run LDpred2
trait = "Example"
library(bigsnpr)
ref_N = "1kg" ###### your reference LD file, can be 1kg and individuals from UKBB

ldr = 3/1000
for (chr in 1:22) {
  sum.raw = bigreadr::fread2(('data/sumdat_Rcov.txt'))
  sum.raw = as.data.frame(sum.raw)
  
  # ------------------------ Run LDpred2
  temfile = paste0('ldpred_ref/chr',chr,'_', ref_N, '.bk')
  if (!file.exists(temfile)){
    snp_readBed(paste0('ldpred_ref/chr',chr,'_', ref_N,'.bed')) ####### path to your reference bed file
  }
  obj.bigSNP <- snp_attach(paste0('ldpred_ref/chr',chr,'_', ref_N,'.rds'))
  
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
  POS2 <- snp_asGeneticPos(CHR, POS, ncores = 2)
  ## indices in info_snp
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 4, infos.pos = POS2[ind.chr2], size = ldr) # default
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  
  # Automatic model
  ldsc <- snp_ldsc2(corr0, df_beta)
  h2_est <- abs(ldsc[["h2"]])
  
  # grid of models:
  H2_seq = signif(abs(h2_est) * c(0.7, 1, 1.4), 3)
  h2_seq <- abs(h2_est) 
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = H2_seq, sparse = c(FALSE))
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 4)
  beta_grid = as.data.frame(beta_grid)
  rownames(beta_grid) = info_snp$rsid
  beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
  colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:nrow(params)))
  
  bigreadr::fwrite2(beta_grid, paste0(trait, '/chr', chr,'/ldpred2effect-hm3-EUR-ref_N', ref_N, '.txt'), sep='\t')
}

############### run ldpred2 using blked################
for (chr in 1:22) {
  sum.raw = bigreadr::fread2('data/sumdat_Rcov.txt')
  sum.raw = as.data.frame(sum.raw)
  
  map_ldref <- readRDS("ldpred2/ldref/map.rds") ##### downloaded from LDpred2 tutorial
  sumstats = sum.raw[sum.raw$CHR == chr, c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  set.seed(2020)
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  info_snp <- snp_match(sumstats, map_ldref)
  (info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  corr_chr <- readRDS(paste0("ldref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3] ##### downloaded from LDpred2 tutorial
  corr <- as_SFBM(corr_chr, compact = TRUE)
  ld = Matrix::colSums(corr_chr^2)
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL)))
  h2_est <- ldsc[["h2"]]
  
  H2_seq = signif(abs(h2_est) * c(0.7, 1, 1.4), 3)
  h2_seq <- abs(h2_est) 
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = H2_seq, sparse = c(FALSE))
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 4)
  rownames(beta_grid) = info_snp$rsid
  beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
  colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:nrow(params)))
  bigreadr::fwrite2(beta_grid, paste0(trait, '/chr', chr,'/ldpred2effect-hm3-EUR-ref_Nblk_provided.txt'), sep='\t')
}


