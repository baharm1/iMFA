# MFA of serine contribution in human tumors

This code was used to estimate the relative serine uptake and synthesis fluxes _in vivo_. A three-compartment model including enhancing, non-enhancing, and brain tissues is implemented to quantify fluxes in patient tumors (see figure below). The metabolic model and methodology is described in the supplementary methods of the manuscript. We used a previously developed IMM-MFA model by [Yang _et al._](https://pubmed.ncbi.nlm.nih.gov/27829138/). 

<p align="center">
	<img width="70%" src="https://github.com/baharm1/iMFA/blob/main/readme_figs/serine_MFA_3comp.png">
</p>

### Requirements
1. MATLAB with default installation (We used MATLAB R2021b on Windows 11 OS; The same files were also run on Linux in some cases)
2. [Artelys Knitro Optimizer](https://www.artelys.com/solvers/knitro/) (MATLAB version 12.4)
3. MATLAB Parallel Processing toolkit (optional)

## Usage
1. To estimate the relative fluxes, change the user-defined inputs in the main script and run the file. Output plots and excel files will be auto-generated and saved in a folder titled `output_files`. The file `flux_ratio.xlsx` will contain the value and confidence intervals for the ratio of de-novo serine synthesis to serine salvage. The file `flux_results.xlsx` will contain the relative flux values.
2. The user-defined inputs are described in the next section. The model inputs should be in the specified format. We have included the input files used to generate the results for the manuscript.
3. Since the three-compartment model (including brain, enhancing, and non-enhancing tumors) was run for seven patients individually (`main_3comp_human`), we saved the model parameters in `3compmodel.mat` and used it for patients 1-7. The code to generate the file is in `Generate_model_file.m`.
4. We use a two-compartment model for patient 8 (including brain and non-enhancing tumor) and for the xenograft mouse models (including brain and GBM). The main scripts for two-compartment models are `main_2comp_human.m` and `main_2comp_mice.m`, respectively.

Please note that the runtime on an HPC cluster with 24 cores and 180 GB RAM was 2-3 days. Runtime may be longer depending on the machine's capability.

### Replication of Data in the Manuscript

The following user-defined inputs were used in all cases:
1. unlabeled_metabs: {'CO2','SERu'}
2. min_sd: sqrt(10^(-5))
3. UBval: 20
4. P: 1
5. nmc: 1000
6. conf: 0.05
7. niter: 200

Fluxes for a 3-compartment model:
1. Load the model information with the load('3compmodel.mat') command.
2. Use the input data for the patient of interest for file2. For instance, file2 = 'patient_input_data/Patient1.xlsx'

Fluxes for a 2-compartment model:
1. file1: 'model_2comp_human.xlsx' or 'model_2comp_mice.xlsx'
2. file2 = 'patient_input_data/Patient8.xlsx' or 'mice_input_data'
3. Note: Running two-compartment model will overwrite f_NL_const.m for the three-compartment model. 

## User-Defined Inputs
1. file1: An excel file with the details of the metabolic model.
2. file2: An excel file with the input data.
3. min_sd: Minimum standard deviation of the isotope tracer data. The standard deviation needs to a non-zero fininite value.
4. unlabeled_metabs: The metabolites that are assumed to be unlabeled in the model.
5. UBval: Upper bound for flux values
6. P: A parameter used in the original IMM-MFA method. Not used for this model.
7. conf: Confidence level for flux bound estimation.
8. niter: Number of initial starting points for parameter optimization.
9. nmc: Number of Monte Carlo simulations for confidence interval estimation.

### Model File

file1 is an Excel file with the following tabs:
1. rxns: List of model reactions. The first column is the reaction number. Second column is reaction type (_E_ for a sink reaction and _I_ for other reactions). The third column has carbon atom transition information.
2. metab_size: The number of carbon atoms in the metabolites. Name and size are in the first and second column respectively.
3. input_metab: The metabolites that are not mass-balanced in the model
4. metab_to_remove: The metabolites that are not included in isotopomer mass balance

### Input Data File

file2 is an Excel file with the following data in a tab labeled 'MID':
1. Column 1: Metabolite name (_e_, _n_, _b_, and _p/x_ represent enhancing, non-enhancing, brain, and plasma)
2. Column 2: Metabolite isotopologue
3. Column 3: average MID value
4. Column 4: standard deviation of measured MID
