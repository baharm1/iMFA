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

