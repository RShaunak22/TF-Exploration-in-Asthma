#this script is for determining the number of TFs conserved across cell-types and specific to cell-types
#dependencies
library(dplyr)
library(stats)
#JASPAR across all cell types with TFs recorded
conserved_JASPAR <- Reduce(intersect, list(
  final_JASPAR_ciliated1$symbol,
  final_JASPAR_ciliated2$symbol,
  final_JASPAR_club$symbol,
  final_JASPAR_goblet$symbol
))
#ChEA3 across all cell types with TFs recorded
conserved_ChEA3 <- Reduce(intersect, list(
  as.character(final_ChEA3_ciliated1$TF),
  as.character(final_ChEA3_ciliated2$TF),
  as.character(final_ChEA3_club$TF),
  as.character(final_ChEA3_goblet$TF),
  as.character(final_ChEA3_mucociliated$TF)
))
database_conserved <-intersect(conserved_ChEA3,conserved_JASPAR)
#alternative approach - since the original aim is to review cell-type TFs, want to have cell-type specific TFs
#want to go back a step actually, want to go to the non-aggregated version with a TF-gene observation in each row
#for JASPAR: final_JASPAR_specific, for ChEA3:since final_ChEA3_TFs and final_ChEA3_unique_TF and final_ChEA3_notnull all have 5649 obs use final_ChEA3_specifi
#JASPAR
ciliated1_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "ciliated 1",] %>% filter(score > 0)
colnames(ciliated1_gene_TF_JASPAR)[1] <- "gene"
ciliated1_KS <- ciliated1_gene_TF_JASPAR %>% left_join(sig_ciliated1, by = "gene") %>%
  dplyr::select("gene","symbol", "score", "p_val_adj", "avg_log2FC")
ciliated1_unique_TF <- unique(ciliated1_gene_TF_JASPAR$symbol)
#start with ciliated 1 cluster for JASPAR
ciliated1_KS_results <- data.frame(matrix(NA, nrow = nrow(final_JASPAR_ciliated1), ncol = 2))
colnames(ciliated1_KS_results) <- c("TF", "p_val")
for (i in seq_along(ciliated1_unique_TF)) {
  #subset target DEGs
  TF <- ciliated1_unique_TF[i]
  TF_target <- ciliated1_KS[ciliated1_KS$symbol == TF,]
  #subset non-target DEGs
  non_target_TF <- ciliated1_KS[ciliated1_KS$symbol != TF,]
  #adding results
  ciliated1_KS_results$TF[i] <- TF
  ciliated1_KS_results$p_val[i] <- ks.test(TF_target$avg_log2FC, non_target_TF$avg_log2FC)$p.val
  
}
#filtering for significant values
ciliated1_KS_results_filtered <- ciliated1_KS_results %>% filter(p_val < 0.05)

#ciliated 2 cluster for JASPAR
ciliated2_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "ciliated 2",] %>% filter(score > 0)
ciliated2_unique_TF <- unique(ciliated2_gene_TF_JASPAR$symbol)
colnames(ciliated2_gene_TF_JASPAR)[1] <- "gene"
ciliated2_KS <- ciliated2_gene_TF_JASPAR %>% left_join(sig_ciliated2, by = "gene") %>%
  dplyr::select("gene","symbol", "score", "p_val_adj", "avg_log2FC")
ciliated2_KS_results <- data.frame(matrix(NA, nrow = nrow(final_JASPAR_ciliated2), ncol = 2))
colnames(ciliated2_KS_results) <- c("TF", "p_val")
for (i in seq_along(ciliated2_unique_TF)) {
  #subset target DEGs
  TF <- ciliated2_unique_TF[i]
  TF_target <- ciliated2_KS[ciliated2_KS$symbol == TF,]
  #subset non-target DEGs
  non_target_TF <- ciliated2_KS[ciliated2_KS$symbol != TF,]
  #recording results
  ciliated2_KS_results$TF[i] <- TF
  ciliated2_KS_results$p_val[i] <- ks.test(TF_target$avg_log2FC, non_target_TF$avg_log2FC)$p.val
  
}
ciliated2_KS_results_filtered <- ciliated2_KS_results %>% filter(p_val < 0.05)

#club cluster for JASPAR
club_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "club",] %>% filter(score > 0)
club_unique_TF <- unique(club_gene_TF_JASPAR$symbol)
colnames(club_gene_TF_JASPAR)[1] <- "gene"
club_KS <- club_gene_TF_JASPAR %>% left_join(sig_club, by = "gene") %>%
  dplyr::select("gene","symbol", "score", "p_val_adj", "avg_log2FC")
club_KS_results <- data.frame(matrix(NA, nrow = nrow(final_JASPAR_club), ncol = 2))
colnames(club_KS_results) <- c("TF", "p_val")
for (i in seq_along(club_unique_TF)) {
  #subset target DEGs
  TF <- club_unique_TF[i]
  TF_target <- club_KS[club_KS$symbol == TF,]
  #subset non-target DEGs
  non_target_TF <- club_KS[club_KS$symbol != TF,]
  #recording results
  club_KS_results$TF[i] <- TF
  club_KS_results$p_val[i] <- ks.test(TF_target$avg_log2FC, non_target_TF$avg_log2FC)$p.val
  
}
club_KS_results_filtered <- club_KS_results %>% filter(p_val < 0.05)

#goblet cluster for JASPAR
goblet_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "goblet",] %>% filter(score > 0)
goblet_unique_TF <- unique(goblet_gene_TF_JASPAR$symbol)
colnames(goblet_gene_TF_JASPAR)[1] <- "gene"
goblet_KS <- goblet_gene_TF_JASPAR %>% left_join(sig_goblet, by = "gene") %>%
  dplyr::select("gene","symbol", "score", "p_val_adj", "avg_log2FC")
goblet_KS_results <- data.frame(matrix(NA, nrow = nrow(final_JASPAR_goblet), ncol = 2))
colnames(goblet_KS_results) <- c("TF", "p_val")
for (i in seq_along(goblet_unique_TF)) {
  #subset target DEGs
  TF <- goblet_unique_TF[i]
  TF_target <- goblet_KS[goblet_KS$symbol == TF,]
  #subset non-target DEGs
  non_target_TF <- goblet_KS[goblet_KS$symbol != TF,]
  #recording results
  goblet_KS_results$TF[i] <- TF
  goblet_KS_results$p_val[i] <- ks.test(TF_target$avg_log2FC, non_target_TF$avg_log2FC)$p.val
  
}
goblet_KS_results_filtered <- goblet_KS_results %>% filter(p_val < 0.05)

#ChEA3
#adapted from JASPAR loops
#slightly different as final_ChEA3_specific has log fold changes appended from ChEA3 script
#ciliated 1 cluster for ChEA3
ciliated1_gene_TF_ChEA3 <- final_ChEA3_specific[final_ChEA3_specific$cell_type == "ciliated 1",] %>%
  dplyr::select("gene","TF", "p_val_adj", "avg_log2FC")
#as.character addition because unique originally stored as factor vector rather than character vector
ciliated1_unique_ChEA3 <- as.character(unique(ciliated1_gene_TF_ChEA3$TF))
#start with ciliated 1 cluster for ChEA3
ciliated1_KS_results_ChEA3 <- data.frame(matrix(NA, nrow = nrow(final_ChEA3_ciliated1), ncol = 2))
colnames(ciliated1_KS_results_ChEA3) <- c("TF", "p_val")
for (i in seq_along(ciliated1_unique_ChEA3)) {
  #subset target DEGs
  TF <- ciliated1_unique_ChEA3[i]
  #col name is TF not symbol
  TF_target <- ciliated1_gene_TF_ChEA3[ciliated1_gene_TF_ChEA3$TF == TF,]
  #subset non-target DEGs
  non_target_TF <- ciliated1_gene_TF_ChEA3[ciliated1_gene_TF_ChEA3$TF != TF,]
  #adding results
  ciliated1_KS_results_ChEA3$TF[i] <- TF
  ciliated1_KS_results_ChEA3$p_val[i] <- ks.test(TF_target$avg_log2FC, non_target_TF$avg_log2FC)$p.val
  
}
ciliated1_ChEA3_KS_filtered <- ciliated1_KS_results_ChEA3 %>% filter(p_val < 0.05)