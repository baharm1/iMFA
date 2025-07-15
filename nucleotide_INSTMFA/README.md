# Quantification of nucleotide metabolism _in vivo_ using iMFA at metabolic steady state
We implemented an ordinary differential equation-based MFA model (INST-MFA) using time-course nucleotide mass isotopologue distribution profiles to to estimate the purine and pyrimidine pathway fluxes (shown in the figure below) in normal cortex and GBM38 tumor-bearing mice. 
The metabolic model and methodology is described in the supplementary methods of the manuscript.
Metabolites were assumed to be not accumulated in the treatment-naïve tissues. 
However, the enrichment of MIDs for mass-balanced metabolites were assumed to be accumulated over time based on our time course experiments. Hence, a metabolic steady state - isotopic non-steady steady state metabolic flux analysis (INST-MFA) was used.

<p align="center">
	<img width="100%" src="https://github.com/baharm1/iMFA/blob/main/readme_figs/nucleotide_INSTMFA.png">
</p>

## Purine INST-MFA Model
1. To estimate the flux profiles, change the user-defined inputs in `main_purines.m` and run the file. Output plots and excel files will be auto-generated and saved in a folder titled output_files. The files intervals.xlsx will contain the optimum flux value and confidence intervals. The folder titled plots will contain MID-time plots that compare the model simulation (blue lines) to the experimental data (orange markers). 95% confidence intervals of fluxes resulted from `main_purines.m` were used as lower and upper bounds of fluxes at t=0 in dynamic iMFA.
2. The user-defined inputs are described in the next section. The model inputs should be in the specified format. We have included the input files used to generate the flux profiles for the manuscript.
3. The file `MTHF_mid_estimation.m` was used to estimate the enrichment of one-carbon unit from the enrichment of serine and glycine. The input_data_file for this small code will be the same as the one used for the main model. The output from this code is included in the input files provided for the main model.
4. The results of purine INST-MFA model are shown in **Fig. 3b**.

Please note that the runtime for our machine with 12 cores, Intel core i-9 processor, and 64 GB RAM was 2-3 days. Runtime may differ depending on the machine's capability.

### Replication of Data in the Manuscript: Purines

The following user-defined inputs were used in all cases:
1. model_file: 'purine_model.xlsx'
2. min_sd: sqrt(10^(-5))
3. ncores: 1
4. alpha: 0.05
5. cl: 0.95

Fluxes in GBM:
1. input_data_file: 'Input-Data\Input_gbm_ctrl_norm.xlsx' (uncomment line 14 and comment line 16 in main_purines.m).
2. Set the bounds for GMP and IMP concentrations (uncomment lines 111-114 and comment lines 80, 90, 117-118 in main_purines.m).

Fluxes in cortex:
1. input_data_file: 'Input-Data\Input_brain_ctrl_norm.xlsx' (uncomment line 16 and comment line 14 in main_purines.m).
2. Set the bounds for IMP concentration (uncomment lines 80, 90, 117, 118 and comment 111-114 in main_purines.m).

## Pyrimidine INST-MFA Model

To estimate the flux profiles, change the user-defined inputs in `main_pyrimidines.m` and run the file. In pyrimidine model, we used a set of parameters for fractional contribution of isotopologues of aspartate that creates labeling patterns in UMP. Other instructions are similar to the purine model.
The results of pyrimidine INST-MFA model are shown in **Fig. 3c**.
   
### Replication of Data in the Manuscript: Pyrimidines

The following user-defined inputs were used in all cases:
1. model_file: 'pyrimidine_model.xlsx'
2. min_sd: sqrt(10^(-5))
3. ncores: 1
4. alpha: 0.05
5. cl: 0.95

Fluxes in GBM:
1. input_data_file: 'Input-Data\Input_gbm_ctrl_norm.xlsx'
2. Set lower and upper bounds of fractional contribution of aspartate in UMP labeling based on [1].

Fluxes in cortex:
1. input_data_file: 'Input-Data\Input_brain_ctrl_norm.xlsx'
2. Set lower and upper bounds of fractional contribution of aspartate in UMP labeling based on [1].

## INST-MFA General Information
User-defined inputs are as follows:
1. model_file: An Excel file with the list of reactions (purine_model.xlsx, pyrimidine_model.xlsx).
2. input_data_file: An excel file with the input data (see files in Input-Data).
3. min_sd: Minimum standard deviation of the isotope tracer data is set to a non-zero fininite value.
4. ncores: Number of optimizations to run in parallel. Possible when multiple Knitro lisences are available (this is diferent from using parallelization within Knitro finite difference implementation).
5. alpha: Threshold used for chi-square-goodness-of-fit test.
6. cl: Confidence level for flux bound estimation.

### Model File

model_file is an Excel file with the following data in separate sheets:
1. reactions: a list of model reactions formatted as in purine_model.xlsx
2. input_metabs: names of the input metabolites
3. balance_metabs: names of the balanced metabolites. Stoichiometry and isotopologue balance equations are solved for these metabolites.

### Input Data File

input_data_file is an Excel file with the following data in separate sheets:
1. MID: The percent MID enrichments of the metabolites in the model calculated based on normalized MIDs (`output_SEnorm`). Values are averaged for each time point and condition. The first column is the name of the metabolite, and the second column is the isotopologue. The subsequent columns contain the MID at the samples time points. The value of the first row from the third column onwards corresponds to the time (in hours). 
2. SD: Sample standard deviation of the percent MID enrichments (`output_SEnorm`) for each time point and condition with the same format as MID sheet. Normalized % mean enrichment and their standard deviations are shown in **Extended Data Fig. 6a,b**.
3. unlabeled_metabs: List of metabolites that are assumed to be unlabeled in the model based on our experimental data.
4. conc: Known metabolite concentrations. The first column has the concentrtaion names, the second has the concentration values and the third has the standard deviation. Concentraion unit is pmol/mg-tissue. Calculation of concentrations are shown in `norm_abundance_control.xlsx` and **Extended Data Fig. 6e**.

## References   
[1] Cai, F. et al. (2023) Comprehensive isotopomer analysis of glutamate and aspartate in small tissue samples. Cell Metab. [10.1016/j.cmet.2023.07.013](https://doi.org/10.1016/j.cmet.2023.07.013)