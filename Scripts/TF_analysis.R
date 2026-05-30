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


#maybe need to readjust conserved_JASPAR and conserved_ChEA3 to a database with DEGs rather than a character list
for (TF in database_conserved) {
  Target <- subset(so_combined@meta.data$rownames == TF) #nope - needs to be genes
}
