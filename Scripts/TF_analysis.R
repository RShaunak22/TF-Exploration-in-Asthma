#this script is for determining the number of TFs conserved across cell-types and specific to cell-types
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
#trialling with one TF
practice <- "TFAP2A"
practice_df <- data.frame()
# JASPAR practice
for (TF in practice){
  target_genes <- final_JASPAR_TFs$symbol_DEG[final_JASPAR_TFs$symbol == TF]
  #since DEGs stored as comma separated string
  # Get the comma-separated string and split it into a character vector
  raw_string <- final_JASPAR_TFs$DEG_symbol[final_JASPAR_TFs$symbol == TF]
  target_genes <- unlist(strsplit(raw_string, ","))
  target_genes <- trimws(target_genes)
  #get expression data for target genes
  exp_target <- FetchData(so_combined, vars = target_genes)
  practice_df <- rbind(practice_df, exp_target)
}
#alternative approach - since the original aim is to review cell-type TFs, want to have cell-type specific TFs
#want to go back a step actually, want to go to the non-aggregated version with a TF-gene observation in each row
#for JASPAR: final_JASPAR_specific, for ChEA3:since final_ChEA3_TFs and final_ChEA3_unique_TF and final_ChEA3_notnull all have 5649 obs use final_ChEA3_specifi
#JASPAR
ciliated1_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "ciliated 1",] %>% filter(score > 0)
ciliated1_unique_TF <- unique(ciliated1_gene_TF_JASPAR$symbol)
ciliated2_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "ciliated 2",] %>% filter(score > 0)
ciliated2_unique_TF <- unique(ciliated2_gene_TF_JASPAR$symbol)
club_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "club",] %>% filter(score > 0)
club_unique_TF <- unique(club_gene_TF_JASPAR$symbol)
goblet_gene_TF_JASPAR <- final_JASPAR_specific[final_JASPAR_specific$cell_type == "goblet",] %>% filter(score > 0)
goblet_unique_TF <- unique(goblet_gene_TF_JASPAR$symbol)
#start with ciliated 1 cluster
for (TF in ciliated1_unique_TF) {
  #subset target DEGs
  TF_target <- ciliated1_gene_TF_JASPAR[ciliated1_gene_TF_JASPAR$symbol == TF,]
  #get a unique target DEG list for the second loop to refer to
  target_unique <- unique(TF_target$symbol_DEG)
    for (DEG in target_unique){
    
    }
  #subset non-target DEGs
  non-target_TF <- final_JASPAR_specific[final_JASPAR_specific$symbol != TF,]
  #get a unique target DEG list for the second loop to refer to
  non_target_unique <- unique(non_target_TF$symbol_DEG)
    for (gene in non_target_unique) {
      
    }
}