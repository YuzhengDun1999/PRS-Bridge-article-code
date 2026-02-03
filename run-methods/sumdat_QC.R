######### this script is for generating summary statistics required for input for every method after quality control.
######### Input summary statistics must have same format as in data/sumdat_Rcov.txt
sumdat = bigreadr::fread2('data/sumdat_Rcov.txt') ###### Change the path to your summary statistics, it must have the same format and name as the example shown in the data/sumdat_Rcov.txt
trait = "Example"
system(paste0("mkdir -p ", trait, "/data/chr{1..22}"))

################### preprocess summary level data ################
library(dplyr)
ldsc_sumdat_all_hm3 = data.frame()
for (chr in 1:22) {
  out_dir = paste0(trait, '/data/chr', chr, '/') ###### 
  sumdat_hm3 = sumdat[sumdat$CHR == chr, ]
  sumdat_hm3 = sumdat_hm3 %>% filter(((BETA/SE)^2/N <= 80)) %>%
    select(CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE, N) %>%
    filter((REF_FRQ>0.01) & (REF_FRQ<0.99)) %>% filter(!(CHR==6 & POS>26e6 & POS<34e6))#hg37
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
out_dir1 = paste0(trait, '/data/')
bigreadr::fwrite2(ldsc_sumdat_all_hm3, paste0(out_dir1, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
