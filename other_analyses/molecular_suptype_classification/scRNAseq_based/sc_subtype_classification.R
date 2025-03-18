# R version 4.4.2
library(Seurat) # version 5.1.0
library(AUCell) # version 1.28.0
library(GSEABase) # version 1.68.0
library(tidyverse) # version 2.0.0
library(ggvenn) # version 0.1.10

# calculate molecular subtype enrichment scores for each cell using AUCell ----
# seurat object containing scRNA-seq data of our patients with cell annotations
# scRNA-seq data will be published in 
# Meghdadi et al., "Digital twins for In Vivo Metabolic Flux Estimations in Cancer Patients”
gbm = readRDS('gbm_w_metadata.rds')

# extract UMI counts for AUCell analysis
gbm_mtx = GetAssayData(gbm, slot = 'counts')

# read signature genes of molecular subtypes defined in 
# Wang et al., Cancer Cell, 2017
subtype_genes = getGmt('molecular_subtype_genes.gmt')

# check gene names to be consistent between gmt and gene expression matrix and
# rename them otherwise
available_genes = subsetGeneSets(subtype_genes, rownames(gbm_mtx))
cbind(nGenes(available_genes))

any("ARNTL" %in% rownames(gbm_mtx)) # BMAL1
any("ACPP" %in% rownames(gbm_mtx)) # ACP3

rownames(gbm_mtx) = replace(rownames(gbm_mtx), 
                            which(rownames(gbm_mtx) %in% c("ARNTL", "ACPP")),
                            c("BMAL1", "ACP3"))

subtype_genes = subsetGeneSets(subtype_genes, rownames(gbm_mtx))
cbind(nGenes(subtype_genes))

# AUCell workflow to calculate geneset scores and assign cells to genesets
cells_ranking = AUCell_buildRankings(gbm_mtx, nCores = 12, plotStats = T)
gc()
cells_AUC = AUCell_calcAUC(subtype_genes, cells_ranking)

set.seed(123)
par(mfrow = c(1, 3))
cells_assignment = AUCell_exploreThresholds(cells_AUC, plotHist = T, 
                                            assignCells = T)

warningMsg = sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg

saveRDS(cells_AUC, 'cells_AUC.rds')
saveRDS(cells_assignment, 'cells_assignment.rds')
# set new thresholds based on histograms
thr_PN = cells_assignment$Proneural$aucThr$thresholds['Global_k1', 'threshold']
thr_CL = cells_assignment$Classical$aucThr$thresholds['Global_k1', 'threshold']
thr_MS = cells_assignment$Mesenchymal$aucThr$thresholds['Global_k1', 'threshold']

new_PN = names(which(getAUC(cells_AUC)['Proneural', ] > thr_PN))
new_CL = names(which(getAUC(cells_AUC)['Classical', ] > thr_CL))
new_MS = names(which(getAUC(cells_AUC)['Mesenchymal', ] > thr_MS))
  
gc()
cellsAssigned = list('Proneural' = new_PN, 
                     'Classical' = new_CL, 
                     'Mesenchymal' = new_MS)
assignmentTable = reshape2::melt(cellsAssigned, value.name = "cell")
colnames(assignmentTable)[2] = "geneSet"
head(assignmentTable)

AUC_scores = cells_AUC@assays@data@listData[["AUC"]]
AUC_scores = t(AUC_scores)
AUC_scores = as.data.frame(AUC_scores)
colnames(AUC_scores)

all(rownames(AUC_scores) == colnames(gbm))

# add AUCell scores to seurat object for easier visualization purposes
gbm = AddMetaData(gbm, metadata = AUC_scores)

pdf(file = 'featurePlot_scwna_molecular_subtypes.pdf', height = 10, width = 12)
FeaturePlot(gbm, features = c('Proneural', 'Classical', 'Mesenchymal'))
dev.off()

# add AUCell cell assignments to seurat object for easier visualization purposes
gbm@meta.data[['AUC_PN']] = 0
gbm@meta.data[['AUC_CL']] = 0
gbm@meta.data[['AUC_MS']] = 0
gbm@meta.data[cellsAssigned$Proneural, 'AUC_PN'] = 1
gbm@meta.data[cellsAssigned$Classical, 'AUC_CL'] = 1
gbm@meta.data[cellsAssigned$Mesenchymal, 'AUC_MS'] = 1

pdf(file = 'dimPlot_AUC_CL_assignment.pdf', height = 4.5, width = 5)
DimPlot(gbm, reduction = 'umap', group.by = 'AUC_PN')
DimPlot(gbm, reduction = 'umap', group.by = 'AUC_CL')
DimPlot(gbm, reduction = 'umap', group.by = 'AUC_MS')
dev.off()

# add patient id and tissue site (Enhancing and NonEnhancing) to seurat object
gbm@meta.data[['patientSite']] = paste(gbm$patient, gbm$site, sep = '_')
levels(gbm@active.ident)

# heatmap of AUCell score for each molecular subtype --------------------------
plot_heatmap_molecular_subtype_group = function(group){
  my_plot = DotPlot(gbm, 
                    features = c('Proneural', 'Classical', 'Mesenchymal'), 
                    group.by = group)
  
  avg_exp_scaled = my_plot$data
  avg_exp_scaled = avg_exp_scaled[-c(2, 5)]
  
  avg_exp_scaled_w = reshape(avg_exp_scaled, idvar = "id", 
                             timevar = 'features.plot', direction = 'wide')
  rownames(avg_exp_scaled_w) = avg_exp_scaled_w$id
  avg_exp_scaled_w$id = NULL
  colnames(avg_exp_scaled_w) = c('Proneural', 'Classical', 'Mesenchymal')
  
  pdf(file = paste('pheatmap_avg_exp_scaled_molecular_subtype_',
                   group, '.pdf', sep = ''),
      height = 4, width = 3)
  pheatmap::pheatmap(as.matrix(avg_exp_scaled_w), 
                     cluster_cols = F, cluster_rows = F,
                     scale = "column", 
                     show_colnames = T, show_rownames = T, 
                     color = colorRampPalette(rev(brewer.pal(11, name = 'RdBu')))(100))
  dev.off()
}

plot_heatmap_molecular_subtype_group('patient')
plot_heatmap_molecular_subtype_group('sample')
plot_heatmap_molecular_subtype_group('patientSite')
gc()

# venn diagram of AUC assignments of molecular subtypes for each patient ------
gbm_split = SplitObject(gbm, split.by = 'patientSite')
for (i in 1:length(gbm_split)){
  print(names(gbm_split[i]))
  df_AUC = gbm_split[[i]]@meta.data[c('AUC_PN', 'AUC_CL', 'AUC_MS')]
  df_AUC = df_AUC %>%
    mutate(across(starts_with("AUC"), as.logical))
  p = ggplot(data = df_AUC) +
    geom_venn(mapping = aes(A = AUC_PN, B = AUC_CL, C = AUC_MS),
              fill_color = c("#7570B3", "#E7298A", "#66A61E"),
              fill_alpha = 0.7) +
    theme_classic() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  ggsave(filename = paste(names(gbm_split[i]), 'AUC_assignment.pdf', sep = "_"),
         device = 'pdf', width = 3, height = 3)
}
