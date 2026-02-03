########## This script is used to generate continuous summary statistics from UK Biobank.
########## Out put sumdat_Rcov.txt should have exactly same format and name as data/sumdat_Rcov.txt
library(dplyr)

########## load your phenotype data ##############
dat = readRDS('pheno.rds') ##### path to your UK Biobank phenotype data
dat = dat[which(dat$f.22020 == "Yes"),] ##### only include independent EUR individuals
dat = dat[which(dat$f.21000.0.0 %in% c("White", "British", "Irish", "Any other white background")),] 

psam_chr1 <- bigreadr::fread2("chr1.psam") #### change to your path to genotype file
col_names <- colnames(dat)

# first ten principal components
top_ten_pc_col <- col_names[str_detect(col_names, "22009")][1:10]

# select the columns to include covariate file
all_columns_to_include <- c("f.eid", "f.31.0.0", "f.21022.0.0", top_ten_pc_col)
dat_sub <- dat[, all_columns_to_include]

# rename the columns for readability
colnames(dat_sub) = c("IID", "sex", "age", paste0("PC",1:10))
cov = dat_sub %>% filter(IID %in% psam_chr1[,2])
cov = cbind(data.frame(FID = cov[,1]), cov)

#### split data into training, tuning, and validation dataset
set.seed(123)
n = nrow(cov)
idx = sample(n)
n_train = floor(0.90 * n)
n_tune = floor(0.05 * n)
train_idx = idx[1:n_train]
tune_idx = idx[(n_train + 1):(n_train + n_tune)]
valid_idx = idx[(n_train + n_tune + 1):n]
cov_train = cov[train_idx, ]
cov_tune = cov[tune_idx, ]
cov_valid = cov[valid_idx, ]

#### extract phenotype
trait = "BMI" ##### change your trait and corresponding field.id, can be searched on UKBB website
field_id = "21001.0.0"
trait_col = col_names[str_detect(col_names, field_id)]
pheno = data.frame(FID = dat[, "f.eid"], IID = dat[, "f.eid"], y = dat[, trait_col])
readr::write_tsv(pheno, paste0('data/pheno.txt'))

readr::write_tsv(cov_train, paste0('data/training_cov.txt'))
cov_tune = merge(cov_tune, pheno, by = c("FID", "IID"))
readr::write_tsv(cov_tune, paste0('tuning/', trait, '_cov.txt'))
cov_valid = merge(cov_valid, pheno, by = c("FID", "IID"))
readr::write_tsv(cov_valid, paste0('validation/', trait, '_cov.txt'))

### run GWAS
GWAS_code = paste(paste0('plink2'),
                  paste0(' --bfile allchr'), #### path to your genotype file
                  paste0(' --pheno data/pheno.txt '),
                  paste0(' --covar data/training_cov.txt  --snps-only '),
                  paste0(' --glm --hwe 1e-07 --geno 0.05  --mind 0.05 --maf 0.01 --mach-r2-filter 0.8 2 '),
                  paste0(' --out data/GWAS'))
system(prscode)
### calculate minor allele frequency
GWAS_code = paste(paste0('plink2'),
                  paste0(' --bfile allchr'), #### path to your genotype file
                  paste0(' --freq --hwe 1e-07 --geno 0.05  --mind 0.05 --maf 0.01 --mach-r2-filter 0.8 2 '),
                  paste0(' --out data/MAF'))
system(prscode)

### process GWAS to required input form
sumdat = bigreadr::fread2("data/GWAS.y.glm.linear")
freq = bigreadr::fread2(paste0('data/MAF.afreq'))
freq = freq %>% filter(ID %in% sumdat$ID)
names(sumdat) = c("CHR", "POS","SNP_ID","REF","ALT", "A1", "TEST", "N","BETA", "SE","T", "PVAL")
sumdat$REF_FREQ = 1-freq[,5]
sumdat = sumdat %>% select(CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE, N)
bigreadr::fwrite2(dat, file = "data/sumdat_Rcov.txt", sep = "\t")
