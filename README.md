


# PRS-Bridge-article-code

  

This is the viable code for the PRS-Bridge paper. 

## Organization:

### evaluation
 - RE_Lassosum.R: get the relative efficiency of each methods in Section S11.
 - get_sd.R: get $R^2$ and standard error of methods in Section 4.2.
 - get_sd_disease.R: get transformed AUC and standard error of methods in Section 4.3.
 
 ### plot-results
 - plot_auto.R: visualize the results of PRS-Bridge-auto in Figure S8 and S9. 
 - plot_PRScs_proj.R: visualize the results of PRS-CS-proj in Figure S7.
 - plot_RE.R: visualize relative efficiency of methods in Figure S10 and S11.
 - plot_tuning.R: visualize results of tuning in Figure S3 and S4.
 - plot.R: visualize results of each method in Figure 4 and 5.
 - plot_fig1.R: Reproduce and visualize posterior coefficients in Figure 1.
 
 ### run-methods
 - disease_pipeline.R: generate script files for running LDpred2, PRS-CS, and PRS-Bridge for disease traits. 
 - LASSOSUM_disease.R: generate script files for running Lassosum for disease traits.
 - LASSOSUM.R: generate script files for running Lassosum for continuous traits in UK Biobank.
 - UKBB_pipeline.R: generate script files for running LDpred2, PRS-CS, and PRS-Bridge for continuous traits in UK Biobank.

 ### run-simulation
 - evaluation.R: evaluation of each methods in simulation studies. 
 - lassosum.R: generate script files for running Lassosum in simulation studies.
 - plot.R: visualize results of each method in Figure 3.

 ### PRS-CS-proj
 This folder contains code revised from the original PRS-CS method
 (https://github.com/getian107/PRScs). The LD reference panel can be downloaded by following the instructions on the original PRS-CS repository. Two main files in this folder are used to reproduce the results shown in Figure 1. An example of how to run the code and generate the figure is provided in plot-result/plot_fig1.R.
 - PRScs_noconverge.py: Remove the ad hoc method implemented in the original PRScs.py 
 - PRScs_proj.R: Extends the original PRS-CS framework by incorporating our proposed projected summary statistics approach.

 ### data
 - PRScs_sumdat.txt: data for reproducing Figure 1.
 - chr22.bim: data for reproducing Figure 1.
