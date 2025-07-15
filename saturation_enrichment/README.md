## Correction of mice inter-variability for U<sup>13</sup>C-glucose saturation

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

## References   
[1] Tracer methods for in vivo kinetics, Chapter 9, Constant infusion of tracer, Shipley R.A. and Clark R.E. (1972)
