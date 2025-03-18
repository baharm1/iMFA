# Differential expression analysis was performed on GEO2R portal.
setwd("from_GEO2R_portal/")
# R version 4.2.2
library(ggplot2) # version 3.4.2
library(EnhancedVolcano) # version 1.16.0
# GSE165595, expression of serine-related genes and neurotransmitters ----
# read norm counts
glioma = read.delim(file = gzfile('GSE165595_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz'),
                    header = T, sep = '\t', row.names = 1)
# convert ensemble IDs to gene symbols
annot_table = read.delim(gzfile('Human.GRCh38.p13.annot.tsv.gz'), 
                         header = T, sep = '\t')

all(rownames(glioma) == annot_table$GeneID)

serine_list = c('PHGDH', 'PSAT1', 'PSPH',
                'SLC1A4', 'SLC1A5', 'SLC3A2',
                'SLC6A9', 'SLC6A14', 'SLC7A8',
                'SLC7A10', 'SLC12A4', 'SLC25A1',
                'SLC25A15', 'SLC25A28', 'SLC36A1',
                'SLC38A1', 'SLC38A2', 'SLC38A4', 
                'SLC38A5', 'SLC38A7')
geneIDs = annot_table[annot_table$Symbol %in% serine_list, 'GeneID'] # comment if you'd like to visualize neurotransmitters

#https://maayanlab.cloud/Harmonizome/gene_set/neurotransmitter/GeneRIF+Biological+Term+Annotations
neurotransmitter_list = c('ADK', 'ADRA2A', 'APP', 'ATP1A1', 'CACNA1A', 'CACNB4',
                          'CCKBR', 'CCL2', 'CHAT', 'CHMP4A', 'CHMP4B', 'COMT',
                          'CPLX2', 'CYP19A1', 'DOC2A', 'DOC2B', 'EN1', 'FEV', 
                          'GABBR1', 'GABRA1', 'GABRB2', 'GAD1', 'GAD2', 'GAL', 
                          'GLS2', 'GNMT', 'HCRTR1', 'HTR5A', 'HTT', 
                          'IL1B', 'MAPT', 'MECP2', 'NPY', 'NR4A2', 'OTOF', 
                          'PRKN', 'PITX3', 'POU5F1', 'PPARG', 'PROM1', 'RAI1',
                          'RGN', 'RNLS', 'RPS6KA3', 'SCG3', 'SLC17A7', 'SLC29A4',
                          'SLC32A1', 'SLC44A1', 'SLC6A15', 'SLC6A2', 'SLC6A4',
                          'SLC9A9', 'SLURP1', 'SNAP29', 'SNCA', 'SNCB', 'SNCG',
                          'SYT1', 'TAC1', 'TP73', 'UNC13A', 'VIP')
# not found in annot_table:
a = annot_table[annot_table$Symbol %in% neurotransmitter_list, 'Symbol']
neurotransmitter_list[!neurotransmitter_list %in% a]
# "HTL", "PARK2" --> NA, PRKN

# comment the following line to visualize serine-related genes
geneIDs = annot_table[annot_table$Symbol %in% neurotransmitter_list, 'GeneID'] 

gene_symbol = annot_table[annot_table$GeneID %in% geneIDs, 'Symbol']

# read metadata and assign tissue status to samples (normal vs tumor)
glioma_metadata = read.delim(file = 'GSE165595_metadata.txt', header = T, 
                             sep = '\t')

gbm_normal = glioma_metadata[(glioma_metadata$Pathology == 'Glioblastoma') & 
                               (glioma_metadata$Tumor.status == 'Normal'), 
                              'Accession']
gbm_tumor = glioma_metadata[(glioma_metadata$Pathology == 'Glioblastoma') & 
                              (glioma_metadata$Tumor.status == 'Tumor'), 
                            'Accession']

gbm_counts = glioma[rownames(glioma) %in% geneIDs, c(gbm_normal, gbm_tumor)]
gbm_counts = log2(gbm_counts + 1) # converts to log-normalized TPM 
rownames(gbm_counts) = gene_symbol
gbm_counts = as.data.frame(t(gbm_counts))

gbm_counts[['group']] = c(rep('Normal', 15), rep('Tumor', 15))

long_gbm_counts = tidyr::pivot_longer(data = gbm_counts, 
                                      cols = 1:length(serine_list))

long_gbm_counts = tidyr::pivot_longer(data = gbm_counts, 
                                      cols = 1:length(neurotransmitter_list))
# visualize log-normalized TPM values
p <- ggplot(long_gbm_counts, aes(x = name, y = value, color = group)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16, position=position_jitterdodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
ggsave(filename = 'boxplot_GSE165595_neurotransmitter.pdf', plot = p, 
       device = 'pdf', width = 20, height = 5)

DEgenes = read.delim('GSE165595_tumor_normal.top.table.tsv', header = T, 
                     sep = '\t')

saveDElist = DEgenes[DEgenes$GeneID %in% geneIDs, ]
write.table(saveDElist, file = 'GSE165595_DE_neurotransmitter_tumor_normal.txt', 
            sep = '\t', col.names = T, row.names = F)

# GSE59612, expression of serine-related genes and neurotransmitters ----
glioma = read.delim(file = gzfile('GSE59612_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz'),
                    header = T, sep = '\t', row.names = 1)
annot_table = read.delim(gzfile('Human.GRCh38.p13.annot.tsv.gz'), 
                         header = T, sep = '\t')

all(rownames(glioma) == annot_table$GeneID)

geneIDs = annot_table[annot_table$Symbol %in% serine_list, 'GeneID']

# not found in annot_table:
a = annot_table[annot_table$Symbol %in% neurotransmitter_list, 'Symbol']
neurotransmitter_list[!neurotransmitter_list %in% a] # NA

# comment the following line to visualize serine-related genes
geneIDs = annot_table[annot_table$Symbol %in% neurotransmitter_list, 'GeneID']

gene_symbol = annot_table[annot_table$GeneID %in% geneIDs, 'Symbol']

# read metadata and assign tissue status to samples (normal vs tumor)
glioma_metadata = read.delim(file = 'GSE59612_metadata.txt', header = T, 
                             sep = '\t')

gbm_NE = glioma_metadata[(glioma_metadata$Tissue.type == 'glioma - non-enhancing FLAIR+ sample'), 
                         'Accession']
gbm_E = glioma_metadata[(glioma_metadata$Tissue.type == 'glioma - contrast-enhancing sample'),
                        'Accession']
gbm_C = glioma_metadata[(glioma_metadata$Tissue.type == 'non-neoplastic brain'),
                        'Accession']

gbm_counts = glioma[rownames(glioma) %in% geneIDs, c(gbm_C, gbm_NE, gbm_E)]
gbm_counts = log2(gbm_counts + 1) # converts to log-normalized TPM 
rownames(gbm_counts) = gene_symbol
gbm_counts = as.data.frame(t(gbm_counts))

gbm_counts[['group']] = c(rep('Cortex', length(gbm_C)), 
                          rep('NE', length(gbm_NE)),
                          rep('E', length(gbm_E)))

long_gbm_counts = tidyr::pivot_longer(data = gbm_counts, 
                                      cols = 1:length(serine_list))
long_gbm_counts = tidyr::pivot_longer(data = gbm_counts, 
                                      cols = 1:length(neurotransmitter_list))
# visualize log-normalized TPM values
p <- ggplot(long_gbm_counts, aes(x = name, y = value, color = group)) + 
  geom_boxplot(outlier.shape = NA, width = 1)+ 
  geom_jitter(shape=16, position=position_jitterdodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(filename = 'boxplot_GSE59612_neurotransmitter.pdf', plot = p, 
       device = 'pdf', width = 25, height = 5)

DEgenes = read.delim('GSE59612_nonenhancing_cortex.top.table.tsv', 
                     header = T, sep = '\t')

saveDElist = DEgenes[DEgenes$GeneID %in% geneIDs, ]
write.table(saveDElist, file = 'GSE59612_DE_neurotransmitter_nonenhancing_cortex.txt', 
            sep = '\t', col.names = T, row.names = F)

# volcano plots ----

# replace neurotransmitter with serine in the following lines to visualize serine-related genes in volcano plot
DElistE = read.delim(file = 'GSE59612_DE_neurotransmitter_enhancing_cortex.txt', 
                      sep = '\t', header = T)
DElistN = read.delim(file = 'GSE59612_DE_neurotransmitter_nonenhancing_cortex.txt', 
                      sep = '\t', header = T)
DElistG = read.delim(file = 'GSE165595_DE_neurotransmitter_tumor_normal.txt', 
                       sep = '\t', header = T)

DElistE[['color']] = rep('#D2691E', nrow(DElistE))
DElistN[['color']] = rep('black', nrow(DElistN))
DElistG[['color']] = rep('gold', nrow(DElistG))
DElist = rbind(DElistE, DElistN, DElistG)

keyvals = DElist$color
names(keyvals)[keyvals == '#D2691E'] <- 'GSE59612_Enhancing_vs_Cortex'
names(keyvals)[keyvals == 'black'] <- 'GSE59612_NonEnhancing_vs_Cortex'
names(keyvals)[keyvals == 'gold'] <- 'GSE165595_GBM_vs_Cortex'

pdf(file = 'volcano_DE_neurotransmitter_FC_1_GSE59612_GSE165595.pdf', 
    height = 20, width = 18)
EnhancedVolcano(DElist,
                lab = DElist$Symbol,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = DElist$Symbol,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                parseLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                arrowheads = FALSE,
                colConnectors = 'black',
                gridlines.minor = FALSE,
                gridlines.major = F) + coord_flip()
dev.off()
