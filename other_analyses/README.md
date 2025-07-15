Requirements of each code were written in the scripts.
## Molecular subtype classification using scRNA-seq
We used scRNA-seq of our patient cohort (from our unpublished manuscript, _Meghdadi, B. et al._ “Digital twins for In Vivo Metabolic Flux Estimations in Cancer Patients”) to identify molecular subtypes of cells for each patient’s tumor tissue. 
A list of signature genes of three molecular subtypes has been provided in [1], which includes 50 genes per molecular subtype (`molecular_subtype_genes.gmt`). 
After integrating all patient samples and checking the quality control using Seurat (version 4.2.0), the UMI counts were used in AUCell (version 1.28.0) to calculate gene set enrichment scores based on the ranking of specified molecular subtype genes among all expressed genes in a cell [2]. 
Following AUCell guidelines, a threshold for each molecular subtype gene set was selected based on the global distribution of area under the curve of gene set rankings (by setting the selected thresholds to Global_k1; dashed grey line in exploreThresholds plots). 
These thresholds were used to assign cells to molecular subtypes. This analysis was performed in R (version 4.4.2).
The results of scRNA-seq subtype classification are shown in **Extended Data Fig. 1b-d**.

## Molecular subtype classification using bulk RNA-seq
Bulk RNA-seq data generated in this study are accessible through GEO Series accession number [GSE299102](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299102). We integrated our RNA-seq data (`read_counts`) with TCGA-GBM RNA-seq data. 
To download TCGA-GBM, we picked sample ids that contain bulk RNA-seq data from [GDC data portal](https://portal.gdc.cancer.gov/). The sample ids can be found in `cohort_Brain.2024-12-13.tsv`. 
Read counts were normalized to reads per kilobase per million mapped reads (RPKM) and transcript per million (TPM) using the human reference genome from [Ensembl](https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz) and the Ensemble transcripts with the maximum transcript length.
RPKM values were transformed to log2(1+RPKM) and Ensemble ids were converted to HGNC ids or gene names using biomaRt library (version 2.54.0) in R (4.2.2). 
Signature genes of molecular subtypes were defined previously in [1]. 
TPM values (`combined_tcga_tempus_tpm.gct`) were imported to GenePattern [3] single sample geneset enrichment analysis (ssGSEA) and sample normalization method was set to rank to calculate geneset enrichment of molecular subtypes (`ssgsea_gene_pattern_tcga_gbm_gdc_tempus`). 
Hierarchical clustering grouped integrated RNA-seq samples based on the similarities of their enrichment scores. 
Since the molecular subtypes of TCGA-GBM were identified previously in [1] (`tcga_subtypes_wang_2017.txt`), the samples from our patient cohort that clustered with the TCGA-GBM samples may have similar molecular subtypes as TCGA-GBM samples. 
To identify molecular subtypes of our patient cohort independent of TCGA-GBM, we imported log-normalized RPKM values (`combined_tcga_tempus_log2_rpkm.txt`) to GlioVis portal [30] and performed three available methods including support vector machine (SVM), k-nearest neighbor (KNN), and ssGSEA (`gliovis_combined_tcga_gbm_gdc_tempus`). A combined assessment of all methods was conducted to assign one subtype to each sample.
The results of bulk RNA-seq subtype classification are shown in **Extended Data Fig. 1a**.

## Transcriptional regulation of neurotransmitters and serine related enzymes in GBM vs cortex
To understand the transcriptional regulation of serine related enzymes and neurotransmitters, we performed differential expression analysis on two available RNA-seq datasets: GSE59612 which includes normal cortex, enhancing and non-enhancing glioma samples [4] and GSE165595 which includes normal cortex and GBM samples [5]. 
Raw counts were imported to DESeq function for differential expression analysis using Wald significance test and size factor type = poscount. Fold changes show the comparison of differential expressions in glioma samples over cortex samples with a false discovery rate less than 0.05 using Benjamini-Hochberg test. 
These calculations were performed on [GEO2R portal](https://www.ncbi.nlm.nih.gov/geo/geo2r/).
Related figures to differential expression analysis are shown in ** Extended Data Fig. 5** for neurotransmitters and serine related genes.

## References
[1] Wang, Q. _et al._ (2017) Tumor Evolution of Glioma-Intrinsic Gene Expression Subtypes Associates with Immunological Changes in the Microenvironment. _Cancer Cell_ 32, 42-56 e46. 10.1016/j.ccell.2017.06.003

[2] Aibar, S. _et al._ (2017) SCENIC: single-cell regulatory network inference and clustering. _Nat Methods_ 14, 1083-1086. 10.1038/nmeth.4463

[3] Reich, M. et al. (2006) GenePattern 2.0. Nat Genet 38, 500-501. 10.1038/ng0506-500

[4] Gill, B.J. et al. (2014) MRI-localized biopsies reveal subtype-specific differences in molecular and cellular composition at the margins of glioblastoma. Proc Natl Acad Sci U S A 111, 12550-12555. 10.1073/pnas.1405839111

[5] Hwang, T. et al. (2022) Genome-wide perturbations of Alu expression and Alu-associated post-transcriptional regulations distinguish oligodendroglioma from other gliomas. Commun Biol 5, 62. 10.1038/s42003-022-03011-w
