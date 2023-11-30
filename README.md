# iMFA

This code was used to estimate the purine and pyrimidine pathway fluxes in GBM and cortex in xenograft mice models. The metabolic model and methodology is described in the supplementary methods of the manuscript. The folder labeled 'Purines' has the model for the purine fluxes. The folder labeled 'Pyrimides' has the model for the pyrimidine fluxes, which is a modified version of the purine model.

### Requirements
1. MATLAB with default installation (We used MATLAB R2021b on Windows 11 OS)
2. [Artelys Knitro Optimizer](https://www.artelys.com/solvers/knitro/) (MATLAB version)
3. MATLAB Parallel Processing toolkit (optional)

## Usage: Purines
1. To estimate the flux profiles, change the user-defined inputs in Main_Purines.m and run the file. Output plots and excel files will be auto-generated and saved in a folder titled output_files. The files intervals.xlsx will contain the optimum flux value and confidence intervals. The folder titled plots will contain MID-time plots that compare the model simulation (lines) to the experimental data (markers).
2. The user-defined inputs are described in the next section. The model inputs should be in the specified format. We have included the input files used to generate the flux profiles for the manuscript.
3. The file 1c_mid_estimate was used to estimate the enrichment of one-carbon unit from the enrichment of serine and glycine. The input_data_file for this small code will be the same as the one used for the main model. The output from this code is included in the input files provided for the main model.

Please note that the runtime for our machine with 12 cores, Intel core i-9 processor, and 64 GB RAM was 2-3 days. Runtime may be longer depending on the machine's capability.

### Replication of Data in the Manuscript: Purines

The following user-defined inputs were used in all cases:
1. model_file: 'purine_model_v1.xlsx'
2. min_sd: sqrt(10^(-5))
3. ncores: 1
4. alpha: 0.05
5. cl: 0.05

Fluxes in GBM:
1. input_data_file: 'Input-Data\Input_gbm_rt_norm.xlsx'
2. Set the bounds for GMP and IMP concentrations (uncomment lines 93-96 in Main_Purines.m).

Fluxes in cortex:
1. input_data_file: 'Input-Data\Input_brain_rt_norm.xlsx'
2. Set the bounds for IMP concentration (uncomment lines 99-100 in Main_Purines.m).

## Usage: Pyrimidines

To estimate the flux profiles, change the user-defined inputs in Main_pyrimidine.m and run the file. Use the inputs specified below. Other instructions are similar to the purine model.
   
### Replication of Data in the Manuscript: Pyrimidines

The following user-defined inputs were used in all cases:
1. model_file: 'pyrimidine_model_v2.xlsx'
2. min_sd: sqrt(10^(-5))
3. ncores: 1
4. alpha: 0.05
5. cl: 0.05

Fluxes in GBM:
1. input_data_file: ??????????
2. ??????Any other specifications?????

Fluxes in cortex:
1. input_data_file: ??????????
2. ??????????

## User-Defined Inputs
1. model_file: An Excel file with the list of reactions (purine_model_v1.xlsx).
2. input_data_file: An excel file with the input data (see files in Input-Data).
3. min_sd: Minimum standard deviation of the isotope tracer data. The standard deviation needs to a non-zero fininite value.
4. ncores: Number of optimizations to run in parallel. Possible when multiple Knitro lisences are available (this is diferent from using parallelization within Knitro finite difference implementation).
5. alpha: Threshold used for chi-square-goodness-of-fit test.
6. cl: Confidence level for flux bound estimation.

### Model File

input_data_file is an Excel file with the following data in separate tabs:
1. reactions: a list of model reactions formatted as in purine_model_v1.xlsx
2. input_metabs: names of the input metabolites
3. balance_metabs: names of the balanced metabolites

### Input Data File

input_data_file is an Excel file with the following data in separate tabs:
1. MID: The MID enrichments of the metabolites in the model (in %). The first column is the name of the metabolite, and the second column is the isotopologue. The subsequent columns contain the MID at the samples time points. The value of the first row from the third column onwards corresponds to the time (in hours).
2. SD: Standard deviation of the measured MID values. Same format as MID.
3. conc: Known metabolite concentrations. The first column has the concentrtaion names, the second has the concentration values and the third has the standard deviation.
4. unlabeled_metabs: List of metabolites that are assumed to be unlabeled in the model.

   



   
