############ We first preprocess summary statistics which could be used for running all methods ########
############ summary statistics, genotype, and phenotype data can be downloaded via https://doi.org/10.7910/DVN/COXHAP
############ Genotype data starts with "EUR" that have .bim, .bed, and .fam extension.
############ Phenotype data looks like EUR_pheno with different parameters for simulating effect size.
############ Summary statistics looks like EUR_summary with different parameters for simulating effect size.

### GA controls the genetic architecture.
### GA=1: strong negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.
### GA=4: no negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.
### GA=5: mild negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.

###rho controls causal SNPs proportion. rho is set to be 1, 2, 3 with the causal SNPs proportions as 0.01, 0.001, 0.0005.
temp <- commandArgs(TRUE)
rho = as.numeric(temp[1]) ###### taking values in 1, 2, 3
GA = as.numeric(temp[2]) ###### taking values in 1, 4, 5

sumdat = bigreadr::fread2(paste0("EUR_summary_rho_", rho, "_size_4_GA_4"))
for (chr in 1:22) {
  out_dir = paste0('rho', rho, 'GA', GA, '/sumdat/chr', chr, '/')
  sumdat_hm3 = sumdat[sumdat[,2]==chr, ]
  
  sumdat_hm3 = sumdat_hm3 %>% filter(((BETA/SE)^2/N <= 80)) %>%
    select(CHR, BP, SNP, A1, A2, P_rep_1, BETA_rep_1, SE_rep_1) %>% filter(!(CHR==6 & BP>26e6 & BP<34e6))#hg37
  names(sumdat_hm3) = c(CHR, POS, SNP_ID, REF, ALT, PVAL, BETA, SE)
  sumdat_hm3$N = 100000
  sumdat_hm3$PVAL = as.numeric(sumdat_hm3$PVAL)
  sumdat_ldsc = sumdat_hm3 %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  bigreadr::fwrite2(sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, SE, N), paste0(out_dir, '/hm3_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  sumdat_ldsc = sumdat_hm3 %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  bigreadr::fwrite2(sumdat_ldsc, paste0(out_dir, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  ldsc_sumdat_all_hm3 = rbind(ldsc_sumdat_all_hm3, sumdat_ldsc)
  sumdat_PRScs = sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, PVAL) %>% 
    rename(SNP=SNP_ID, A1=REF, A2=ALT, BETA=BETA, P=PVAL)
  bigreadr::fwrite2(sumdat_PRScs, paste0(out_dir, '/PRScs_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
}

