########## This script is used to preprocess disease summary statistics downloaded from public website to the required input pipeline.
########## Out put sumdat_Rcov.txt should have exactly same format and name as data/sumdat_Rcov.txt
library(dplyr)

########## pre-process Rheumatoid Arthritis summary statistics ##########
dat = bigreadr::fread2("RA_GWASmeta_European_v2.txt.gz")
names(dat) = c("SNP", "CHR", "BP", "A1", "A2", "OR", "CI_low", "CI_up", "P")
sumdat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 43290) %>%
  mutate(SE = (log(CI_up) - log(CI_low)) / 2 / 1.96) %>% mutate(FRQ = 0.1) %>%
  select(CHR, BP, SNP, A1, A2, FRQ, P, BETA, SE, N)
names(sumdat) = c("CHR","POS","SNP_ID","REF","ALT","REF_FRQ","PVAL","BETA","SE","N")
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = "\t")


########### pre-process Inflammatory Bowel disease summary statistics ########
dat = bigreadr::fread2("EUR.IBD.gwas_info03_filtered.assoc")
dat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 34652) %>% mutate(SE = SE/OR) %>%
  select(CHR, BP, SNP, A1, A2, FRQ_A_12882, P, BETA, SE, N)
names(dat) = c("CHR","POS","SNP_ID","REF","ALT","REF_FRQ","PVAL","BETA","SE","N")
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = "\t")


########### pre-process Coronary Artery disease summary statistics ########
dat = bigreadr::fread2('cad.add.160614.website.txt')
dat = dat %>% select(chr, bp_hg19, markername, effect_allele, noneffect_allele, effect_allele_freq, p_dgc, beta, se_dgc)
names(dat) = c("CHR", "POS", "SNP_ID", "REF", "ALT", "REF_FRQ", "PVAL", "BETA", "SE")
dat$N = 184305
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = "\t")


########### pre-process Depression summary statistics ########
dat = bigreadr::fread2("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB")
dat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 2 * Neff_half) %>% 
  select(CHR, BP, SNP, A1, A2, FRQ_A_45396, P, BETA, SE, N)
names(dat) = c("CHR", "POS", "SNP_ID", "REF", "ALT", "REF_FRQ", "PVAL", "BETA","SE", "N")
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = ",")


########## pre-process breast cancer summary statistics ##########
dat = bigreadr::fread2('oncoarray_bcac_public_release_oct17.txt', na.strings = "NULL",
                       select = c("chr", "position_b37", "phase3_1kg_id", "a0", "a1", "bcac_onco_icogs_gwas_eaf_controls",
                                  "bcac_onco_icogs_gwas_P1df", "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se"),
                       col.names = c("chr", "pos", "SNP_ID", "a0", "a1", "freq", "p", "beta", "beta_se"))
SNP_ID = sapply(dat[,2], function(x){return(strsplit(x,split=":")[[1]][1])})
names(SNP_ID) = NULL
dat$SNP_ID = SNP_ID
dat = dat %>% filter(!is.na(SNP_ID)) %>% filter(SNP_ID!="1")
names(dat) = c("CHR", "POS", "SNP_ID", "REF", "ALT", "REF_FRQ", "PVAL", "BETA", "SE", "PVAL")
dat = dat %>% select(CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE)
dat$N = 4 / (1 / 137045 + 1 / 119078)
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = ",")
