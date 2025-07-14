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