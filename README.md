# iMFA: _In Vivo_ Metabolic Flux Analysis

This code was used to estimate the purine and pyrimidine pathway fluxes in GBM and cortex in xenograft mice models. The metabolic model and methodology is described in the supplementary methods of the manuscript.
Metabolites were assumed to be not accumulated in the treatment-naïve tissues. 
However, the enrichment of MIDs for mass-balanced metabolites were assumed to be accumulated over time based on our time course experiments. Hence, a metabolic steady state - isotopic non-steady steady state metabolic flux analysis (INST-MFA) was used.
In the radiation treated tissues, we incorporated time-dependent changes of metabolite pool sizes and a dynamic MFA was performed.
Codes related to molecular subtype classification and differential expression analysis can be found in `other_analyses` folder. 

### Requirements
1. MATLAB with default installation (We used MATLAB R2021b on Windows 11 OS)
2. [Artelys Knitro Optimizer](https://www.artelys.com/solvers/knitro/) (MATLAB version)
3. MATLAB Parallel Processing toolkit (optional)

## Correction of mice inter-variability for U<sup>13</sup>C-glucose saturation (`saturation_enrichment`)

To remove the inter-variability of plasma glucose enrichment between mice, we fitted the plasma glucose mean enrichment to a two-compartment exponential decay function based on the tracer kinetics model [1]. The following parameters are needed to run `two_compartment.m`:

1. Rate of tracer infusion (mg/min/g): `r = 0.012`
2. Bolus dose (mg/min/g): `P = 0.4`
3. Plasma glucose mean enrichment at steady state to estimate parameters of fitting curve: `plasma_glucose_ME_4h.xlsx`
4. Plasma glucose mean enrichment for all samples to estimate saturation enrichment: `plasma_glucose_ME_all.xlsx`
5. Minimum steady state enrichment observed in the data: `min_ss = 0.35`
6. Maximum enrichment from bolus at t = 0: `t0_max = 0.6`
7. Minimum enrichment from bolus at t = 0: `t0_min = 0.3`

Output files of `two_compartment.m` are as follows:

1. `avgSE_2comp.xlsx`: asymptotes of fitted curve i.e. saturation enrichment for all samples
2. `avgSE_fit_2comp`: plots of plasma glucose enrichment with fitted curve

To remove the inter-variability in saturation enrichment of mice, a correction factor is calculated for each mouse by dividing the mouse saturation enrichment to the average of saturation enrichment values across mice (`SE_factors.xlsx`). 
Then the tissue labeled MIDs are divided by the correction factor in `SEnorm.R` which requires the following the input files:

1. `input_data_perform_SEnorm`: tissue metabolite MIDs
2. `SE_factors.xlsx`: correction factors for all samples

Output files of `SEnorm.R` are normalized MIDs (`output_SEnorm`).

**Data Exclusion Criteria**

The need for multiple timepoints to develop our hypothesis-generating metabolic flux models constrained sample numbers at each timepoint, and exclusions were necessary to minimize variance within groups and ensure data robustness. To ensure the reliability of our data for these models, we excluded data points (highlighted in red in excel files) based on study-specific criteria tailored to the limited sample size and proof-of-concept nature of this modeling. Samples were removed if: (1) Data disrupted the expected concave, increasing enrichment trend (Extended Data Fig. 6a,b), in particular at t=60; (2) technical replicates exhibited high variability (exceeding within-sample consistency thresholds); (3) MIDs of a replicate deviated significantly from the other replicates; or (4) tissue weight measurements used for abundance normalization (Extended Data Fig. 6h) were outliers (either excessively low or high).

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
2. Set lower and upper bounds of fractional contribution of aspartate in UMP labeling based on [2].

Fluxes in cortex:
1. input_data_file: 'Input-Data\Input_brain_ctrl_norm.xlsx'
2. Set lower and upper bounds of fractional contribution of aspartate in UMP labeling based on [2].

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

## Purine Dynamic MFA
This code was used to estimate time-dependent changes in purine pathway fluxes after radiation. We used B-splines to generate the transient flux profiles from time-course isotope tracing data. 
The metabolic model and methodology is described in the supplementary methods of the manuscript. The results of purine dynamic MFA are shown in **Extended Data Fig. 7**.

### Usage
1. To estimate the flux profiles, change the user-defined inputs in `main_purines_dmfa.m` and run the file. The flux-time plots will be generated and saved a folder titled output_files.
2. The user-defined inputs are described in the next section. The model inputs should be in the specified format. We have included the input files used to generate the flux profiles in the manuscript.
3. If the optimal spline knots are to be determined, then run Optimize_Breakpoints.m and choose the optimal knots according to the generated score plots.
4. To run the code for DMFA, we first need to run the INST-MFA method on steady-state data to get the flux profiles before treatment. The results are used to set the initial condition for the dynamic flux profiles. This helped start optimization from an informed parameter space. 95% confidence intervals of fluxes resulted from purine INST-MFA were used as lower and upper bounds of fluxes at t=0 in dynamic iMFA.

Please note that the runtime for our machine with 12 cores, Intel core i-9 processor, and 64 GB RAM was 3-4 days. Runtime may differ depending on the machine's capability.

### User-Defined Inputs
1. model_file: An Excel file with the list of reactions simialr to INST-MFA model_file (purine_model.xlsx).
2. input_data_file: An excel file with the input data (see files in Input-Data).
3. min_sd: Minimum standard deviation of the isotope tracer data is set to a non-zero fininite value.
4. niter: Number of random starting points for parameter optimization.
5. bspline_order: spline order
6. internal_knots: position of internal knots for spline
7. ub_flux: upper bound for the flux parameters
8. l2: weight for l2 normalization of flux parameters (optional). We set l2 to 0 for the fluxes reported in the manuscript.

### Input Data File

input_data_file is an Excel file with the following data in separate sheets:
1. MID: The percent MID enrichments of the metabolites in the model calculated based on normalized MIDs (`output_SEnorm`). Values are averaged for each time point and condition. The first column is the name of the metabolite, and the second column is the isotopologue. The subsequent columns contain the MID at the samples time points. The value of the first row from the third column onwards corresponds to the time (in hours).
2. SD: Sample standard deviation of the percent MID enrichments (`output_SEnorm`) for each time point and condition with the same format as MID sheet. MIDs and their standard deviations are shown in **Extended Data Fig. 6g**.
3. conc: The fold-change in metabolite concentrations relative to t=0 (_i.e._, ctrl condition). The first column corresponds to the metabolite name and the subsequent columns correspond to the values at sampled time points. The value of the first row from the third column onwards corresponds to the time (in hours).
4. conc SD: Standard deviation of the measured values in conc sheet with the same format as conc sheet. Calculation of concentrations are shown in `norm_abundance_rt.xlsx` and **Extended Data Fig 6h**.
5. unlabeled_metabs: List of metabolites that are assumed to be unlabeled in the model based on our experimental data.
6. c0: The initial concentrtaion of the mass-balanced metabolites. The values were set according to the results from the INST-MFA purine model. conc: concentration, sd: standard deviation, lb: lower bound, ub: upper bound, ss_opt: The optimum value determined from the INST-MFA code.
7. v0: The initial fluxes. The values were set according to the results from the INST-MFA purine model. lb: 95% lower confidence bound, ub: 95% upper confidence bound, ss_opt: The optimum value determined from the INST-MFA code.

### Replication of Data in the Manuscript

The following user-defined inputs were used in all cases:
1. model_file: 'purine_model.xlsx'
2. min_sd: sqrt(10^(-5))
3. niter: 100
4. bspline_order: 3
5. ub_flux: 2000
6. l2: 0

Fluxes in GBM after RT:
1. input_data_file: 'Input-Data\Input_gbm_rt_conc_abund.xlsx'
2. internal_knots: [0.3]

Fluxes in cortex after RT:
1. input_data_file: 'Input-Data\Input_brain_rt_conc_abund.xlsx'
2. internal_knots: [0.65]

## References   
[1] Tracer methods for in vivo kinetics, Chapter 9, Constant infusion of tracer, Shipley R.A. and Clark R.E. (1972)

[2] Cai, F. et al. (2023) Comprehensive isotopomer analysis of glutamate and aspartate in small tissue samples. Cell Metab. [10.1016/j.cmet.2023.07.013](https://doi.org/10.1016/j.cmet.2023.07.013)
   
