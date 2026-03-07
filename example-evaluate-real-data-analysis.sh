###### provide example code to run pipeline under the section evaluate Real Data Results

TRAIT="example"
OUTCOME="continuous"
PATH_TO_GENO="example/UKBB/allchr"

Rscript evaluate-real-data-analysis/get_sd_LASSOSUM.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_PRSBridge.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_PRScs.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_PRScs_proj.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_PRScs_regularized.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_PRScs_threshold.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}
Rscript evaluate-real-data-analysis/get_sd_ldpred2.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}

Rscript plot-result/plot_example.R