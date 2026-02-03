


# PRS-Bridge-article-code

  

This is the viable code for the PRS-Bridge paper. 

## Organization:

### evaluation
 - RE_Lassosum.R: get the relative efficiency of each methods in Section S11.
 - get_sd.R: get $R^2$ and standard error of methods in Section 4.2.
 - get_sd_disease.R: get transformed AUC and standard error of methods in Section 4.3.
 
 ### plot-results
 - plot.R: visualize results of each method in Figure 4 and 5.
 - plot_fig1.R: Reproduce and visualize posterior coefficients in Figure 1.
 - plot_auto.R: visualize the results of PRS-Bridge-auto in Figure S8 and S9. 
 - plot_RE.R: visualize relative efficiency of methods in Figure S10 and S11.
 - plot_tuning.R: visualize results of tuning in Figure S3 and S4.
 
 ### run-methods
 To reproduce the real data analysis results for disease traits, please first download the publicly available GWAS summary statistics and run process_disease_dat.R. For continuous traits, obtain the genotype and phenotype data from the UK Biobank and run process_continuous_dat.R. Next, run sumdat_QC.R to perform quality control and generate the required input summary statistics for each method. Finally, follow the guidance in each method’s tutorial to obtain the appropriate LD reference panel and run the corresponding R scripts for each method.
 - process_disease_dat.R: preprocess disease summary statistics downloaded from website.
 - process_continuous_dat.R: generate and preprocess GWAS summary statistics for continuous trait in UK Biobank.
 - sumdat_QC.R: Quality control for summary statistics.
 - PRScs.R: generate script files for running PRS-CS. 
 - PRScs_proj.R: generate script files for running PRS-CS-Projection. 
 - PRSBridge.R: generate script files for running PRS-Bridge.
 - PRSBridge.R: R script files for running LDpred2.
 - LASSOSUM.R: R script files for running LASSOSUM.
 
 
 ### run-simulation
 To reproduce the simulation results in Section 4.1, please first download the simulated GWAS summary statistics, phenotype data, and genotype files from https://doi.org/10.7910/DVN/COXHAP. Next, run data_proc.R to preprocess the downloaded summary statistics. Then, in the run-methods directory, run PRS-Bridge with the 1kg_small reference data, PRS-CS with the 1000G reference data, and LDpred2 with the 1000G reference data. Finally, run evaluation.R to compute PRS and $R^2$ values and plot.R to generate the corresponding figures.
 - data_proc.R: process downloaded simulated GWAS summary statistics to the required input form for each methods.
 - evaluation.R: evaluation of each methods in simulation studies. 
 - plot.R: visualize results of each method in Figure 3.

 ### PRS-CS-proj
 This folder contains code revised from the original PRS-CS method
 (https://github.com/getian107/PRScs). The LD reference panel can be downloaded by following the instructions on the original PRS-CS repository. Two main files in this folder are used to reproduce the results shown in Figure 1. An example of how to run the code and generate the figure is provided in plot-result/plot_fig1.R.
 - PRScs_noconverge.py: Remove the ad hoc method implemented in the original PRScs.py 
 - PRScs_proj.R: Extends the original PRS-CS framework by incorporating our proposed projected summary statistics approach.

 ### data
 - PRScs_sumdat.txt: data for reproducing Figure 1.
 - chr22.bim: data for reproducing Figure 1.
 - sumdat_Rcov.txt: Example of summary statistics for running pipeline.
