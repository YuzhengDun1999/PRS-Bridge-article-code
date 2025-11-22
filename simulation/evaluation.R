cal_R2 = function(dat){
  model = lm(y ~ PRS, data = dat)
  return(summary(model)$r.squared)
}
####### lassosum
temp <- commandArgs(TRUE)
rho = as.numeric(temp[1])
#size = as.numeric(temp[4])
GA = as.numeric(temp[2])
#rep = as.numeric(temp[6])
size = 4
rep = 1

library('dplyr')
args = commandArgs(trailingOnly = TRUE)
method = 'Lassosum'
result_dir = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/LASSOSUM/rho', rho, '-size', size, '-rep', rep, '-GA', GA)

result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1

PRS_list_tuning = lapply(1:22, function(chr) {
  file_path <- bigreadr::fread2(paste0(result_dir, '/chr', chr, "/tuning_prs.txt"))
})
cov = bigreadr::fread2(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/EUR/pheno_summary_out_GA/phenotypes_rho', rho, '_',GA, '.phen'))
cov = cov[, c(1,3)]; names(cov) = c("IID", "y")
cov = cov %>% filter(IID %in% PRS_list_tuning[[1]][,1])
for (x in 1:80) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS_tmp = PRS_list_tuning[[chr]][, c(1, x + 1)]; names(PRS_tmp) = c("IID", "SCORESUM")
    cov = merge(cov, PRS_tmp, by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:100) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, summary(lm(y~PRS, data = tuning))$r.squared, summary(lm(y~PRS, data = validation))$r.squared); result_i = result_i + 1
  }
}
#readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result/", method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_validation_result_std.txt'))


####### PRS-Bridge ########
percent_list = c(0, 0.2, 0.4, 0.6, 0.8)
alpha_list = c(0.5, 0.25, 0.125, 0.0625)
method = "Bridge"
result_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result_projection_prior_GA", GA, "/")

result = data.frame(matrix(ncol = 5, nrow = 0))
names(result) = c("seed", "alpha", "percent", "tuning_R2", "validation_R2"); result_i = 1

cov = bigreadr::fread2(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/EUR/pheno_summary_out_GA/phenotypes_rho', rho, '_',GA, '.phen'))
cov = cov[, c(1,3)]; names(cov) = c("IID", "y")
for (alpha in alpha_list) {
  for (percent in percent_list) {
    cov$PRS = 0
    for (chr in 1:22) {
      PRS = bigreadr::fread2(paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result_projection_prior_GA", GA, "/chr", chr, '/rho', rho, "/percent", percent,"/alpha",alpha, '/prsresult.txt.profile'))
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
#readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref,'_all_result.txt'))
readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result/", method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_all_result_std.txt'))

result_validation = data.frame(matrix(ncol = 4, nrow = 0))
names(result_validation) = c("seed", "alpha", "percent", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$alpha, result_tmp$percent, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
alpha_best = result_validation %>% group_by(alpha) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(alpha_best = alpha_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result/", method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_validation_result_std.txt'))


############ Ldpred2 ##########
method = "ldpred2"
result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
cov = bigreadr::fread2(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/EUR/pheno_summary_out_GA/phenotypes_rho', rho, '_',GA, '.phen'))
cov = cov[, c(1,3)]; names(cov) = c("IID", "y")

for (x in 1:63) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/ldpred2/simulation/result/chr', chr, '/ldpred2effect-hm3-EUR-rho',rho, '-size',size,'-rep',rep,'-GA',GA, '_prs_', x,'.txt.profile'))
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:100) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, cal_R2(tuning), cal_R2(validation)); result_i = result_i + 1
  }
}
#readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref,'_all_result.txt'))
result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result/", method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_validation_result_std.txt'))

############## PRS-CS ###########
phi_list = c(0.01, 0.0001, 0.000001)
name_list = c('1e-02', '1e-04', '1e-06')
method = "PRScs"

result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "phi", "tuning_R2", "validation_R2"); result_i = 1
cov = bigreadr::fread2(paste0('/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/EUR/pheno_summary_out_GA/phenotypes_rho', rho, '_',GA, '.phen'))
cov = cov[, c(1,3)]; names(cov) = c("IID", "y")

for (phi_i in 1:length(phi_list)){  
  phi = phi_list[phi_i]; name = name_list[phi_i]
  cov$PRS = 0
  for (chr in 1:22) {
    if (GA == 5) {PRS = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRScs/simulation_result/chr', chr, '-rho',rho, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '_prs.txt.profile'))}
    if (GA != 5) {PRS = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRScs/simulation_result/chr', chr, '-rho',rho, '-GA',GA, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '_prs.txt.profile'))}
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:100) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, phi, cal_R2(tuning), cal_R2(validation)); result_i = result_i + 1
  }
}
#readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref,'_all_result.txt'))
result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/simulation/result/", method, "-rho", rho, '-size', size, '-rep', rep, '-GA', GA, '_validation_result_std.txt'))



