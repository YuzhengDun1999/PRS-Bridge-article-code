########## This script is used to preprocess disease summary statistics downloaded from public website to the required input pipeline.
########## Out put sumdat_Rcov.txt should have exactly same format and name as data/sumdat_Rcov.txt
library(dplyr)

temp <- commandArgs(TRUE)
trait = temp[1] # phenotype name, same for the whole pipeline
pheno = temp[2] # path to your phenotype file pheno.rds
psam_path = temp[3] # path to your psam file of any chromosome, chr1.psam

system(paste0("mkdir ", trait))

########## load your phenotype data ##############
dat = readRDS(pheno) ##### path to your UK Biobank phenotype data
dat = dat[which(dat$f.22020 == "Yes"),] ##### only include independent EUR individuals
dat = dat[which(dat$f.21000.0.0 %in% c("White", "British", "Irish", "Any other white background")),] 

psam_chr1 <- bigreadr::fread2(psam_path) #### change to your path to genotype file
col_names <- colnames(dat)

# first ten principal components
top_ten_pc_col <- col_names[str_detect(col_names, "22009")][1:10]

# select the columns to include covariate file
all_columns_to_include <- c("f.eid", "f.31.0.0", "f.21022.0.0", top_ten_pc_col)
dat = dat %>% filter(f.eid %in% psam_chr1[,2])
dat$FID = dat$f.eid; dat$IID = dat$f.eid
dat_sub <- dat[, all_columns_to_include]

# rename the columns for readability
colnames(dat_sub) = c("IID", "sex", "age", paste0("PC",1:10))
cov = cbind(data.frame(FID = cov[,1]), cov)

get_incident = function(self_report_id, ICD_code, trait, self_report_field= "f.20002"){
  self_report = dat[, str_detect(names(dat), self_report_field)]
  selected_report = c()
  for (j in 1:ncol(self_report)) {
    selected_report = c(selected_report, which(self_report[, j] %in% self_report_id))
  }
  selected_report = unique(selected_report)
  
  ICD = dat[, str_detect(names(dat), "f.41270")]
  ICD_date = dat[, str_detect(names(dat), "f.41280")]
  recruit_date = dat$f.53.0.0
  selected_ICD = c()
  for (j in 1:ncol(ICD)) {
    selected_ICD = c(selected_ICD, which(ICD[, j] %in% ICD_code))
  }
  selected_ICD = unique(selected_ICD)
  
  case_index = intersect(selected_ICD, selected_report)
  incident_case = rep(1, length(case_index))
  for (i in 1:length(case_index)) {
    index = case_index[i]
    ICD_index = which(ICD[index,] %in% ICD_code)
    for (j in 1:length(ICD_index)) {
      if (ICD_date[index, ICD_index[j]] < recruit_date[index]){
        incident_case[i] = 0
      }
    }
  }
  
  set.seed(2023)
  case_incident_index = case_index[incident_case == 1]
  case_sample = sample(0:1, size = length(case_incident_index), replace = TRUE, prob = c(0.5,0.5))
  case_tuning_index = case_incident_index[case_sample == 0]
  case_validation_index = case_incident_index[case_sample == 1]
  case_tuning = dat[case_tuning_index,] %>% mutate(y = 1) %>% select(c("FID", "IID", y)) 
  case_validation = dat[case_validation_index,] %>% mutate(y = 1) %>% select(c("FID", "IID", y)) 
 
  control_index_C = unique(union(selected_ICD, selected_report))
  control_index = rep(TRUE, nrow(dat)); control_index[control_index_C] = FALSE
  control_all = dat[control_index,] %>% mutate(y = 0) %>% select(c("FID", "IID", y)) 
  
  control_sample = sample(0:1, size = nrow(control_all), replace = TRUE, prob = c(0.5,0.5))
  control_tuning = control_all[control_sample == 0,]; control_tuning = control_tuning[1:10000,]
  control_validation = control_all[control_sample == 1,]; control_validation = control_validation[1:10000,]
  
  tuning = rbind(case_tuning, control_tuning)
  validation = rbind(case_validation, control_validation)
  cov_tune = merge(cov_tune, tuning, by = c("FID", "IID"))
  readr::write_tsv(cov_tune, paste0('tuning/', trait, '_cov.txt'))
  cov_valid = merge(cov_valid, validation, by = c("FID", "IID"))
  readr::write_tsv(cov_valid, paste0('validation/', trait, '_cov.txt'))
}


########## pre-process Rheumatoid Arthritis summary statistics ##########
if (trait == "RA") {
  dat = bigreadr::fread2("RA_GWASmeta_European_v2.txt.gz")
  names(dat) = c("SNP", "CHR", "BP", "A1", "A2", "OR", "CI_low", "CI_up", "P")
  sumdat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 43290) %>%
    mutate(SE = (log(CI_up) - log(CI_low)) / 2 / 1.96) %>% mutate(FRQ = 0.1) %>%
    select(CHR, BP, SNP, A1, A2, FRQ, P, BETA, SE, N)
  names(sumdat) = c("CHR","POS","SNP_ID","REF","ALT","REF_FRQ","PVAL","BETA","SE","N")
  bigreadr::fwrite2(dat, file = paste0(trait, "/sumdat_Rcov.txt"), sep = "\t")
  get_incident(1464, paste0("M0", c(50, 51, 52, 53, 58, 59, 60, 61, 62, 63, 64, 68, 69)), trait)
}

########### pre-process Inflammatory Bowel disease summary statistics ########
if (trait == "IBD") {
  dat = bigreadr::fread2("EUR.IBD.gwas_info03_filtered.assoc")
  dat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 34652) %>% mutate(SE = SE/OR) %>%
    select(CHR, BP, SNP, A1, A2, FRQ_A_12882, P, BETA, SE, N)
  names(dat) = c("CHR","POS","SNP_ID","REF","ALT","REF_FRQ","PVAL","BETA","SE","N")
  bigreadr::fwrite2(dat, file = paste0(trait, "/sumdat_Rcov.txt"), sep = "\t")
  get_incident(1461, paste0("K", c(500, 501, 508, 509, 510, 511, 512, 513, 514, 515, 518, 519)), trait)
}

########### pre-process Coronary Artery disease summary statistics ########
if (trait == "CAD") {
  dat = bigreadr::fread2('cad.add.160614.website.txt')
  dat = dat %>% select(chr, bp_hg19, markername, effect_allele, noneffect_allele, effect_allele_freq, p_dgc, beta, se_dgc)
  names(dat) = c("CHR", "POS", "SNP_ID", "REF", "ALT", "REF_FRQ", "PVAL", "BETA", "SE")
  dat$N = 184305
  bigreadr::fwrite2(dat, file = paste0(trait, "/sumdat_Rcov.txt"), sep = "\t")
  get_incident(c(1074, 1075), paste0("I", c(210, 211, 212, 213, 214, 219, "21X", 220, 221, 228, 229, 231:238, 220, 240, 241, 248, 249, 250:259)), trait)
}

########### pre-process Depression summary statistics ########
if (trait == "Depression") {
  dat = bigreadr::fread2("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB")
  dat = dat %>% mutate(BETA = log(OR)) %>% mutate(N = 2 * Neff_half) %>% 
    select(CHR, BP, SNP, A1, A2, FRQ_A_45396, P, BETA, SE, N)
  names(dat) = c("CHR", "POS", "SNP_ID", "REF", "ALT", "REF_FRQ", "PVAL", "BETA","SE", "N")
  bigreadr::fwrite2(dat, file = paste0(trait, "/sumdat_Rcov.txt"), sep = ",")
  get_incident(1286, paste0("F", c(320, 321, 322, 323, 328, 329, 330, 331, 332, 333, 334, 338, 339)))
}

########## pre-process breast cancer summary statistics ##########
if (trait == "BC") {
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
  bigreadr::fwrite2(dat, file = paste0(trait, "/sumdat_Rcov.txt"), sep = ",")
  get_incident(1002, paste0("F", c(500:509)), trait, "f.20001")
}