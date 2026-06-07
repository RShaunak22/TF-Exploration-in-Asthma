#dependencies
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(reshape2)
library(SummarizedExperiment)
#setup
#reference genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#extracting genes
txdb_genes <- genes(txdb)
#pfm acquisition from jaspar
pfm_matrix <- getMatrixSet(JASPAR2022, 
                           opts = list(collection = "CORE", 
                                       tax_group = "vertebrates", 
                                       all_versions = FALSE))
#storing TF metadata from JASPAR to get TF names
tf_metadata <- data.frame(
  jaspar_id = sapply(pfm_matrix, ID),
  symbol = sapply(pfm_matrix, name)
)

#converting gene names to location names (ENTREZ ID)
loc_DEGs <- select(org.Hs.eg.db, 
                        keys=rownames(sig_results_total), 
                        keytype="SYMBOL", 
                        columns=c("ENTREZID"))
loc_DEGs$cell_type <- sig_results_total[match(loc_DEGs$SYMBOL, 
                                              rownames(sig_results_total)), 
                                        "cell_type"]
#check that genes exist
loc_DEGs <- loc_DEGs[!is.na(loc_DEGs$ENTREZID), ]
#find signficant DEGs in ENTREZ ID format which map to txdb genome
loc_DEGs_valid <- loc_DEGs[loc_DEGs$ENTREZID %in% names(txdb_genes), ]
valid_DEGs <- loc_DEGs_valid$ENTREZID
#if ENTREZ IDs do not intersect for a cell type
if (length(valid_DEGs) ==0) {
    return(NULL)
  }
#define promoter regions for DEGs in the txdb genome
promoters_DEGs <- promoters(genes(txdb, filter = list(gene_id = valid_DEGs)),
                            upstream = 2000, 
                            downstream = 200)
#find motifs in JASPAR which bind to this promoter region
motif_matches_DEGs <- matchMotifs(pfm_matrix, 
                              promoters_DEGs, 
                              genome = BSgenome.Hsapiens.UCSC.hg38,
                              out = "scores")
#extract prediction scores and convert them into a dataframe format
scores <- assay(motif_matches_DEGs)
scores_df <- as.data.frame(as.matrix(scores))
#rows = DEG, cols = TF
#so DEG = rownames
scores_df$DEG <- rownames(scores_df)
#need to convert object to long format to merge with final dataframe
sm_long <- melt(scores_df,
                id.var = "DEG",
                variable.name = "jaspar_id",
                value.name = "score")
#merge with metadata df to get relevant TF names
final_df <- merge(tf_metadata, sm_long, by="jaspar_id")
#convert back to gene symbols
genes_symbol <- select(org.Hs.eg.db,
                            keys = final_df$DEG,
                            keytype = "ENTREZID",
                            columns = "SYMBOL")
final_df$symbol_DEG <- genes_symbol$SYMBOL[match(final_df$DEG, genes_symbol$ENTREZID)]
#assign cell types to DEG even if there are repeats
cell_types <- data.frame(
  symbol_DEG = loc_DEGs_valid$SYMBOL,
  cell_type = loc_DEGs_valid$cell_type
)
#merge dataframes
final_df <- merge(final_df, cell_types, by = "symbol_DEG", all.x = TRUE)

#remove null entries
tf_results <- Filter(Negate(is.null), final_df)
#remove duplicate rownames
tf_results_unique <- final_df[!duplicated(final_df), ]
#no change - now try gene-TF interactions
tf_results_unique <- tf_results_unique[!duplicated(tf_results_unique[c("symbol","symbol_DEG","cell_type")]), ]
#check if these TFs are present in the original seurat object
final_JASPAR_present <- tf_results_unique[tf_results_unique$symbol %in% rownames(so_combined), , drop = FALSE]
#like ChEA3 reduced but still rather large number of observations therefore sense check
sense_check_JASPAR <- final_JASPAR_present[final_JASPAR_present$symbol_DEG %in% rownames(so_combined), , drop = FALSE]
#sense check passed: no DEG info added
#check that cell-type trends are retained
final_JASPAR_specific <- final_JASPAR_present %>%
  group_by(cell_type) %>%
  filter({
    cells <- rownames(so_combined@meta.data[so_combined@meta.data$cell_types == unique(cell_type), ])
    present_TFs <- rownames(so_combined)[rowSums(so_combined@assays$RNA$counts[, cells, drop = FALSE]) > 0]
    symbol %in% present_TFs
  }) %>%
  ungroup()

#now want to group (using aggregate function) gene info and have unique TFs
#notes: gene ~ TF in many ~ one relationship
final_JASPAR_TFs <- aggregate(symbol_DEG ~ symbol + cell_type, data = final_JASPAR_specific, FUN = function(x) paste(x, collapse = ","))

final_JASPAR_basal1 <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "Basal 1",]
final_JASPAR_basal2 <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "Basal 2",]
final_JASPAR_ciliated1 <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "ciliated 1",]
final_JASPAR_ciliated2 <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "ciliated 2",]
final_JASPAR_club <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "club",]
final_JASPAR_goblet <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "goblet",]
final_JASPAR_mucociliated <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "Mucociliated",]
final_JASPAR_ionocytes <- final_JASPAR_TFs[final_JASPAR_TFs$cell_type == "Ionocytes",]
