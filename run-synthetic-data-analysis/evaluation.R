### GA controls the genetic architecture.
### GA=1: strong negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.
### GA=4: no negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.
### GA=5: mild negative selection, common SNPs heritability fixed as 0.4 for each of the five ancestries, cross-ancestry genetic architecture as 0.8.

###rho controls causal SNPs proportion. rho is set to be 1, 2, 3 with the causal SNPs proportions as 0.01, 0.001, 0.0005.

cal_R2 = function(dat){
  model = lm(y ~ PRS, data = dat)
  return(summary(model)$r.squared)
}

temp <- commandArgs(TRUE)
rho = as.numeric(temp[1])
GA = as.numeric(temp[2])

library('dplyr')
############### PRS-Bridge ###############
percent_list = c(0.2, 0.4, 0.6, 0.8)
alpha_list = c(0.5, 0.25, 0.125)
LD_size = "small" # input what LD_size you use in run-methods/PRSBridge.R, can be "large", and "small"
ref = "1kg" # input what LD reference data you use in run-methods/PRSBridge.R
method = paste0("Bridge_", LD_size)
trait =  paste0("rho", rho, "GA", GA)
result_dir = paste0(trait, "/", ref, "/", method, "/result/")

result = data.frame(matrix(ncol = 5, nrow = 0))
names(result) = c("seed", "alpha", "percent", "tuning_R2", "validation_R2"); result_i = 1

cov = bigreadr::fread2(paste0('EUR_phen_rho_', rho, '_GA_',GA))
cov = cov[100001:120000, c(1,3)]; names(cov) = c("IID", "y")
for (alpha in alpha_list) {
  for (percent in percent_list) {
    cov$PRS = 0
    for (chr in 1:22) {
      bfile = paste0('validation/EUR_chr', chr) #### path to your validation simulated bfile
      alpha = alpha_list[alpha_i]; percent = percent_list[percent_i]
      output_dir = paste0(result_dir, "chr", chr, "/percent",percent,"/alpha",alpha)
      prs_output = paste0(output_dir, '/prsresult.txt')
      prsscore = paste0(output_dir, '/coef.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 3 4 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
      PRS = bigreadr::fread2(paste0(prs_output, '.profile'))
      cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
    }
    for (seed in 1:100) {
      set.seed(seed)
      idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
      tuning = cov[idx, ]; validation = cov[-idx, ]
      result[result_i,] = c(seed, alpha, percent, cal_R2(tuning), cal_R2(validation)); result_i = result_i + 1
    }
  }
}
readr::write_tsv(result, paste0(method, "-rho", rho, '-size4-rep1-GA', GA, '_all_result_std.txt'))

result_validation = data.frame(matrix(ncol = 4, nrow = 0))
names(result_validation) = c("seed", "alpha", "percent", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$alpha, result_tmp$percent, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
alpha_best = result_validation %>% group_by(alpha) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(alpha_best = alpha_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0(method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_validation_result_std.txt'))


############ Ldpred2 ##########
method = "ldpred2"
ref_N = "1kg"
result_dir = paste0(trait, "/", ref_N, "/", method, "/result/")
result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
cov = bigreadr::fread2(paste0('EUR_phen_rho_', rho, '_GA_',GA))
cov = cov[100001:120000, c(1,3)]; names(cov) = c("IID", "y")

for (x in 1:63) {
  cov$PRS = 0
  for (chr in 1:22) {
    bfile = paste0('validation/EUR_chr', chr) #### path to your validation simulated bfile
    prs_output = paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_prs_', x,'.txt')
    prscode = paste(paste0('plink'),
                    paste0('--score ', result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_x_', x, '.txt'),
                    paste0(' 1 3 4',' sum --bfile ',bfile),
                    paste0(' --out ', prs_output))
    system(prscode)
    PRS = bigreadr::fread2(paste0(prs_output, '.profile'))
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:100) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, cal_R2(tuning), cal_R2(validation)); result_i = result_i + 1
  }
}
result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0(method, "-rho", rho, '-size4-rep1-GA', GA, '_validation_result_std.txt'))

############## PRS-CS ###########
phi_list = c(1, 0.01, 0.0001, 0.000001)
name_list = c('1e+00', '1e-02', '1e-04', '1e-06')
method = "PRScs"
ref = "1kg"
result_dir = paste0(trait, "/", ref, "/", method, "/result/")

result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "phi", "tuning_R2", "validation_R2"); result_i = 1
cov = bigreadr::fread2(paste0('EUR_phen_rho_', rho, '_GA_',GA))
cov = cov[100001:120000, c(1,3)]; names(cov) = c("IID", "y")

for (phi_i in 1:length(phi_list)){  
  phi = phi_list[phi_i]; name = name_list[phi_i]
  cov$PRS = 0
  for (chr in 1:22) {
    bfile = paste0('validation/EUR_chr', chr) #### path to your validation simulated bfile
    phi = phi_list[phi_i]; name = name_list[phi_i]
    prs_output = paste0(result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '_prs.txt')
    prscode = paste(paste0('plink'),
                    paste0('--score ', result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '.txt'),
                    paste0(' 2 5 6 sum --bfile ', bfile),
                    paste0(' --out ', prs_output))
    system(prscode)
    PRS = bigreadr::fread2(paste0(prs_output, '.profile'))
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:100) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, phi, cal_R2(tuning), cal_R2(validation)); result_i = result_i + 1
  }
}
result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0(method, "-rho", rho, '-size4-rep1-GA', GA, '_validation_result_std.txt'))



