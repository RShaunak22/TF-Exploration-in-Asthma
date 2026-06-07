# ORIGINAL DATABASE SETUP

#dependencies and setup 
require(Seurat)
require(data.table)
library(dplyr)
library(R.utils)

#ChEA3 cross-referencing
#BiocManager::install("GSEABase")
library(GSEABase)
all_tissues_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/all_tissues.TFs.gmt")
View(all_tissues_TFs)
#Citation: Keenan AB, Torre D, Lachmann A, Leong AK, Wojciechowicz M, Utti V, Jagodnik K, 
#Kropiwnicki E, Wang Z, Ma'ayan A (2019) ChEA3: transcription factor enrichment 
#analysis by orthogonal omics integration. Nucleic Acids Research. doi: 10.1093/nar/gkz446
all_tissues_list <- setNames(geneIds(all_tissues_TFs), names(all_tissues_TFs))
#convert to dataframe
all_tissues_df <- stack(all_tissues_list)
colnames(all_tissues_df) <- c("gene","TF")
View(all_tissues_df)

#cross-referencing with diff exp genes
tf_crossover_chea3_all <- merge(sig_results_total, all_tissues_df, by="gene", all=FALSE)

#lungs TF
lung_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/lung.TFs.gmt")
View(lung_TFs)
lungs_TFs_list <- setNames(geneIds(lung_TFs), names(lung_TFs))
lung_df <- stack(lungs_TFs_list)
colnames(lung_df) <- c("gene", "TF")
tf_crossover_lung <- merge(sig_results_total, lung_df, by="gene", all=FALSE)

#coexpression TFs
coexp_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/ARCHS4_Coexpression.gmt")
View(coexp_TFs)
coexp_TFs_list <- setNames(geneIds(coexp_TFs), names(coexp_TFs))
coexp_df <- stack(coexp_TFs_list)
colnames(coexp_df) <- c("gene", "TF")
coexp_df$TF <- factor(gsub("_ARCHS4_PEARSON","",coexp_df$TF))
tf_crossover_coexp <- merge(sig_results_total, coexp_df, by="gene", all=FALSE)


# ENCODE
encode_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/ENCODE_ChIP-seq.gmt")
View(encode_TFs)
encode_TFs_list <- setNames(geneIds(encode_TFs), names(encode_TFs))
encode_df <- stack(encode_TFs_list)
colnames(encode_df) <- c("gene", "TF")
encode_df$TF <- factor(gsub("_.*","",encode_df$TF))
tf_crossover_encode <- merge(sig_results_total, encode_df, by="gene", all=FALSE)

#ENRICHR
enrichr_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/Enrichr_Queries.gmt")
View(enrichr_TFs)
enrichr_TFs_list <- setNames(geneIds(enrichr_TFs), names(enrichr_TFs))
enrichr_df <- stack(enrichr_TFs_list)
colnames(enrichr_df) <- c("gene", "TF")
tf_crossover_enrichr <- merge(sig_results_total, enrichr_df, by="gene", all=FALSE)

#GTEX
gtex_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/GTEx_Coexpression.gmt")
View(gtex_TFs)
gtex_TFs_list <- setNames(geneIds(gtex_TFs), names(gtex_TFs))
gtex_df <- stack(gtex_TFs_list)
colnames(gtex_df) <- c("gene", "TF")
tf_crossover_gtex <- merge(sig_results_total, gtex_df, by="gene", all=FALSE)

#Literature
lit_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/Literature_ChIP-seq.gmt")
View(lit_TFs)
lit_TFs_list <- setNames(geneIds(lit_TFs), names(lit_TFs))
lit_df <- stack(lit_TFs_list)
colnames(lit_df) <- c("gene", "TF")
lit_df$TF <- factor(gsub("_.*","",lit_df$TF))
tf_crossover_lit<- merge(sig_results_total, lit_df, by="gene", all=FALSE)

#Remap
remap_TFs <- getGmt("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/ChEA3/ReMap_ChIP-seq.gmt")
View(remap_TFs)
remap_TFs_list <- setNames(geneIds(remap_TFs), names(remap_TFs))
remap_df <- stack(remap_TFs_list)
colnames(remap_df) <- c("gene", "TF")
tf_crossover_remap <- merge(sig_results_total, remap_df, by="gene", all=FALSE)


#merging to final dataframe
final_ChEA3_df <- bind_rows(tf_crossover_chea3_all, tf_crossover_coexp, tf_crossover_encode,
                           tf_crossover_enrichr,tf_crossover_gtex, tf_crossover_lit, 
                         tf_crossover_lung, tf_crossover_remap)

#removing duplicates from dataframe
final_ChEA3_df_unique <- final_ChEA3_df[!duplicated(final_ChEA3_df[, c("TF", "gene","cell_type")]), ]
#assuming this has worked because observations for each dataframe have reduced.
#now looking for TFs that are present in original dataset
final_ChEA3_present <- final_ChEA3_df_unique[final_ChEA3_df_unique$TF %in% rownames(so_combined), , drop = FALSE]
#reduced but still rather large number of observations therefore sense check
sense_check <- final_ChEA3_present[final_ChEA3_present$gene %in% rownames(so_combined), , drop = FALSE]
#sense check passed. Genes in final_ChEA3_present list are also found in so_combined object
#i.e. only genes which are present in dataset are in this list
#it occurred to me that cell-type specific TFs in my seurat object may be reflected 
#across multiple cell-types in my dataset if they are present there. 
#However, I want TFs which reflect cell-type trends in my seurat object
#for each cell type, find cell-type specific TFs
ct_list <-unique(sense_check$cell_type)
final_ChEA3_specific <- do.call(rbind, lapply(ct_list, function(ct) {
  #extract cells for each cell type
  cells <- rownames(so_combined@meta.data[so_combined@meta.data$cell_types == ct,])
  #find cell types with expression levels > 0 
  present_TFs <- rownames(so_combined)[rowSums(so_combined@assays$RNA$counts[, cells, drop = FALSE]) > 0]
  #subset TFs with only present TFs
  subset_ct <- sense_check[sense_check$cell_type == ct, ]
  subset_ct[subset_ct$TF %in% present_TFs, ]
}))
#further reduced but still rather large list!
#inspired by JASPAR script, look for null values
final_ChEA3_notnull <- Filter(Negate(is.null), final_ChEA3_specific)
#no change
final_ChEA3_unique_TF <- final_ChEA3_notnull[!duplicated(final_ChEA3_notnull[,c("TF","cell_type")]), ]
#now want to group (using aggregate function) gene info and have unique TFs
#notes: gene ~ TF in many ~ one relationship
final_ChEA3_TFs <- aggregate(gene ~ TF + cell_type, data = final_ChEA3_unique_TF,FUN = function(x) paste(x, collapse = ","))

final_ChEA3_basal1 <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "Basal 1",]
final_ChEA3_basal2 <- final_ChEA3_TFs[final_ChEA3_unique_TF$cell_type == "Basal 2",]
final_ChEA3_mucociliated <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "mucociliated",]
final_ChEA3_ciliated1 <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "ciliated 1",]
final_ChEA3_ciliated2 <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "ciliated 2",]
final_ChEA3_goblet <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "goblet",]
final_ChEA3_club <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "club",]
final_ChEA3_ionocytes <- final_ChEA3_TFs[final_ChEA3_TFs$cell_type == "ionocytes",]
