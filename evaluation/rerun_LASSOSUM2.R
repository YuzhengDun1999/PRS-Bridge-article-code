run_sh_list = c()
for (trait in c("BMI", "HDL", "LDL", "APOEA", "APOEB", "APOEB", "resting_heart_rate", "IBD", "BC", "RA", "Depression", "CAD")) {
  for (ref in c("1kg", "5000", "blk_provided")) {
    if (ref == "1kg") {R_script = paste0("Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/LASSOSUM2_trait_1kg.R")}
    if (ref == "5000") {R_script = paste0("Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/LASSOSUM2_trait.R")}
    if (ref == "blk_provided") {R_script = paste0("Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/LASSOSUM2_trait_blk.R")}
    for (chr in 1:22) {
      filepath = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/LASSOSUM2/chr", chr, "/LASSOSUM2effect-hm3-EUR-ref_N", ref, ".txt")
      if (!file.exists(filepath)) {
        sink(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/rerun_LASSOSUM2_trait", trait, "_ref", ref, "_chr", chr, ".sh"))
        sent_run = paste0(R_script, " ", chr, " ", trait)
        cat('#!/bin/bash', "\n")
        cat("module load conda_R \n")
        cat(sent_run, "\n")
        sink()
        run_sh_list = c(run_sh_list, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/rerun_LASSOSUM2_trait", trait, "_ref", ref, "_chr", chr, ".sh"))
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/run_all.sh")
sink(run_all)
mem = 100
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --time=2-00 --cpus-per-task=2 --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-063,compute-062,compute-094,compute-065,compute-069 ', run_sh)
  cat(submit_sh, "\n")
}
sink()


run_sh_list = c()
for (trait in c("BMI", "HDL", "LDL", "APOEA", "APOEB", "APOEB", "resting_heart_rate", "IBD", "BC", "RA", "Depression", "CAD")) {
  for (ref in c("1kg", "5000")) {
    for (chr in 1:22) {
      filepath = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/LASSOSUM/chr", chr, "/LASSOSUMeffect-hm3-EUR-ref_N", ref, ".txt")
      if (!file.exists(filepath)) {
        sink(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/rerun_LASSOSUM_trait", trait, "_ref", ref, "_chr", chr, ".sh"))
        if (ref == "1kg") {sent_run = paste0("Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/LASSOSUM_trait.R ", chr, " ", trait, " ", ref)}
        if (ref == "5000") {sent_run = paste0("Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/LASSOSUM_trait.R ", chr, " ", trait, " ", ref)}
        cat('#!/bin/bash', "\n")
        cat("module load conda_R \n")
        cat(sent_run, "\n")
        sink()
        run_sh_list = c(run_sh_list, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/rerun_LASSOSUM_trait", trait, "_ref", ref, "_chr", chr, ".sh"))
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/code/rerun_Lassosum/run_all.sh")
sink(run_all)
mem = 15
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --cpus-per-task=2 --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-063,compute-062,compute-094,compute-065,compute-069 ', run_sh)
  cat(submit_sh, "\n")
}
sink()
