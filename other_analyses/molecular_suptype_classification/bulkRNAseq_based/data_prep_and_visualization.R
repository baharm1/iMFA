# prep combined TCGA and Tempus data frame ------------------------------------
# R version 4.2.2
library(biomaRt) # version 2.54.0
library(rtracklayer) # version 1.58.0
library(dplyr) # version 1.1.2
## read Tempus data ----
P1 = read.delim(file = 'read_counts/htseq_P1.txt', header = FALSE, row.names = 1)
P3 = read.delim(file = 'read_counts/htseq_P3.txt', header = FALSE, row.names = 1)
P7 = read.delim(file = 'read_counts/htseq_P7.txt', header = FALSE, row.names = 1)

ensembl_ids = rownames(P1)
tail(ensembl_ids)
ensembl_ids = head(ensembl_ids, -5)

combined_P = cbind(P1, P3, P7)
colnames(combined_P) = c('P1', 'P3', 'P7')
combined_P['__ambiguous', ] / combined_P['__no_feature', ]
combined_P = combined_P[ensembl_ids, ]

## read TCGA downloaded from GDC portal (gdc.cancer.gov) -----------------------
tcga_path = 'tcga_gbm_gdc_rna/gdc_download_20241213_210040.770237.tar.gz'

# Create a temporary directory to extract the contents
temp_dir <- tempdir()

# Extract the contents of the tar.gz file to the temporary directory
untar(tcga_path, exdir = temp_dir)

# List the extracted files
extracted_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
#print(extracted_files)

# Keep tsv files only 170 files
extracted_files_2 <- extracted_files[grepl('_star_gene_counts.tsv', 
                                           extracted_files)]

# Extract file ids
split_file_path = stringr::str_split(extracted_files_2, pattern = '/')
file_ids = lapply(seq(1, length(split_file_path)), 
                  function(x) split_file_path[[x]][2])

file_ids = unlist(file_ids)

# create a data frame with gene names
read_file = read.table(extracted_files_2[1], sep = '\t', header = T)
read_file = read_file[-c(1, 2, 3, 4), ] 
tcga_ens_ids = read_file$gene_id
tcga_gene_names = read_file$gene_name

split_ens = stringr::str_split(tcga_ens_ids, pattern = '\\.')
split_ens = lapply(seq(1, length(split_ens)), function(x) split_ens[[x]][1])
split_ens = unlist(split_ens)

# common ensemble ids between our patients data and TCGA 
common_ens = intersect(split_ens, ensembl_ids)

# Read the individual files and save them all in another variable
com_files_unstranded = data.frame('gene_id' = tcga_ens_ids)

for (i in 1:length(extracted_files_2)) {
  file_id = file_ids[i]
  read_file = read.table(extracted_files_2[i], sep = '\t', header = T)
  # remove the first 4 rows: 
  #N_unmapped, N_multimapping, N_noFeature, N_ambiguous
  read_file = read_file[-c(1, 2, 3, 4), ] 
  com_files_unstranded[[file_id]] = read_file$unstranded
}  

rownames(com_files_unstranded) = com_files_unstranded$gene_id
com_files_unstranded$gene_id = NULL


com_files_unstr_common = com_files_unstranded[split_ens %in% common_ens, ]
combined_P_common = combined_P[rownames(combined_P) %in% common_ens, ]

a = stringr::str_split(rownames(com_files_unstr_common), pattern = '\\.')
a = lapply(seq(1, length(a)), function(x) a[[x]][1])
a = unlist(a)

rownames(com_files_unstr_common) = make.names(a, unique = T)
com_files_unstr_common = com_files_unstr_common[common_ens, ]

all(rownames(com_files_unstr_common) == rownames(combined_P_common))

# combine our data with TCGA
com_tcga_tempus = cbind(combined_P_common, com_files_unstr_common)

# import gtf file
gtf_file = import.gff(gzfile('Homo_sapiens.GRCh38.113.gtf.gz'))
gtf_df = data.frame('gene_id' = gtf_file$gene_id,
                    'transcript_length' = gtf_file@ranges@width)

# There are multiple transcript lengths for one ensemble id. Make it unique for
# Ensemble transcript with longest size
gtf_df_longest = gtf_df %>%
  group_by(gene_id) %>%
  dplyr::slice(which.max(transcript_length)) %>%
  ungroup()

gtf_df_longest = as.data.frame(gtf_df_longest)
# check order of genes in data and gtf
all(gtf_df_longest$gene_id == ensembl_ids)

# convert ensemble ids to hgnc ids
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = common_ens,
                   mart = ensembl)

b = unique(gene_names$ensembl_gene_id)

gene_hgnc = gene_names[!(gene_names$ensembl_gene_id %in% c('ENSG00000230417', 
                                                           'ENSG00000280739',
                                                           'ENSG00000286061')), ]

rownames(gene_hgnc) = gene_hgnc$ensembl_gene_id
com_tcga_tempus[["ensembl_id"]] = rownames(com_tcga_tempus)

res = merge(com_tcga_tempus, gene_hgnc, 
            by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE) # 59379, 175

res = na.omit(res) # 59376, 175
res[res$hgnc_symbol == "", 'hgnc_symbol'] = res[res$hgnc_symbol == "", 'ensembl_id']

rownames(gtf_df_longest) = gtf_df_longest$gene_id 
gtf_common_genes = gtf_df_longest[res$ensembl_id, ]

res$ensembl_id = NULL
rownames(res) = make.names(res$hgnc_symbol, unique = T)
res$hgnc_symbol = NULL

## convert read counts to tpm and rpkm ----
total_reads = colSums(res)

numinator = res / gtf_df_longest$transcript_length
denuminator = colSums(numinator)

# create input file for Gliovis
rpkm_P = sweep(numinator, 2, total_reads, `/`) * 10 ^ 9
log2_rpkm_P = log2(1 + rpkm_P)

write.table(log2_rpkm_P, file = 'combined_tcga_tempus_log2_rpkm.txt', 
            row.names = T, col.names = T, quote = F, sep = '\t')
			
# create input file for gene pattern ssGSEA analysis  
tpm_P = sweep(numinator, 2, denuminator, `/`) * 10 ^ 6
tpm_P[['Description']] = rownames(res)
tpm_P = tpm_P %>% relocate(Description)

write.table(tpm_P, file = 'combined_tcga_tempus_tpm.gct', 
            row.names = T, col.names = T, quote = F, sep = '\t')

# Fig S1A ---------------------------------------------------------------------
# read ssGSEA gene pattern output
tcga_out = read.table(file = 'ssgsea_gene_pattern_combined_tcga_gbm_gdc_tempus/combined_tcga_tempus_tpm.PROJ.txt', 
                      sep = '\t', header = T, row.names = 1)
t_tcga_out = as.data.frame(t(tcga_out))

hist(t_tcga_out$Proneural)
hist(t_tcga_out$Classical)
hist(t_tcga_out$Mesenchymal)

scaled_com_tcga = t(scale(t(as.data.frame(tcga_out)))) # scale over samples

patient_ids = read.delim(file = 'tcga_gbm_gdc_rna/gdc_sample_sheet.2024-12-13.tsv', 
                         header = T, quote = "")
tcga_subtypes = read.delim(file = 'tcga_gbm_gdc_rna/tcga_subtypes_wang_2017.txt', 
                           header = T, quote = "")

patient_ids$Sample.ID = gsub('-', '.', patient_ids$Sample.ID)
patient_ids$Sample.ID = gsub('.{1}$', '', patient_ids$Sample.ID)

rownames(patient_ids) = patient_ids$File.ID
patient_ids = patient_ids[colnames(com_files_unstr_common),]

rownames(tcga_subtypes) = tcga_subtypes$sampleId

tcga_subtypes_sorted = tcga_subtypes[patient_ids$Sample.ID, ]

# read output of log2rpkm run on gliovis
gliovis = read.delim(file = 'GlioVis_three_way_log2_rpkm_tcga_tempus.txt', 
                     header = T, quote = "")

gliovis$fileID[1:115] = gsub('^.', '', gliovis$fileID[1:115])
rownames(gliovis) = gliovis$fileID

patient_ids = read.delim(file = 'tcga_gbm_rna/gdc_sample_sheet.2024-12-13.tsv', 
                         header = T, quote = "")
rownames(patient_ids) = patient_ids$File.ID
patient_ids = patient_ids[colnames(com_files_unstr_common),]
rownames(patient_ids) = make.names(patient_ids$Sample.ID, unique = T)
 
gliovis[['sampleID']] = patient_ids$Sample.ID
gliovis[['subtype']] = tcga_subtypes_sorted$Group
rownames(gliovis) = rownames(patient_ids)

ann_colors = list(
  gsea.subtype.call = c(Proneural = "#7570B3", Classical = "#E7298A", 
                        Mesenchymal = "#66A61E"),
  knn.subtype.call = c(Proneural = "#7570B3", Classical = "#E7298A", 
                       Mesenchymal = "#66A61E"),
  svm.subtype.call = c(Proneural = "#7570B3", Classical = "#E7298A", 
                       Mesenchymal = "#66A61E"),
  subtype = c(PN = "#7570B3", CL = "#E7298A", 
              MS = "#66A61E"))
			  
tempus_gliovis = data.frame('fileID' = c('P1', 'P3', 'P7'),
                            'gsea.subtype.call' = c('Proneural', 'Proneural', 'Proneural'),
                            'knn.subtype.call' = c('Classical', 'Classical', 'Classical'),
                            'svm.subtype.call' = c('Proneural', 'Proneural', 'Classical'),
                            'sampleID' = c('P1', 'P3', 'P7'),
                            'subtype' = c(NA, NA, NA))

rownames(tempus_gliovis) = tempus_gliovis$sampleID
annotation_col = rbind(tempus_gliovis, gliovis)
annotation_col$fileID = NULL
annotation_col$sampleID = NULL
annotation_col$gsea.subtype.call = factor(annotation_col$gsea.subtype.call)
annotation_col$knn.subtype.call = factor(annotation_col$knn.subtype.call)
annotation_col$svm.subtype.call = factor(annotation_col$svm.subtype.call)
annotation_col$subtype = factor(annotation_col$subtype)

colnames(scaled_com_tcga) = rownames(annotation_col)

# heatmap with hierarchical clustering and Gliovis results 
pdf(file = 'scaled_combined_tcga_tempus_tpm_gliovis_subtypes_log2rpkm.pdf', 
    height = 5, width = 30)
pheatmap::pheatmap(as.matrix(scaled_com_tcga), cluster_rows = F, 
                   cluster_cols = T, show_colnames = T,
                   show_rownames = T, annotation_col = annotation_col,
                   annotation_colors = ann_colors)
dev.off()
