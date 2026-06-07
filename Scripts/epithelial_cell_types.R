#The workflow for this script follows the methodology as described by Satija et al as cited below:
#@Article{,
  #author = {Yuhan Hao and Tim Stuart and Madeline H Kowalski and Saket Choudhary and Paul Hoffman and Austin Hartman and Avi Srivastava and Gesmira Molla and Shaista Madad and Carlos Fernandez-Granda and Rahul Satija},
 # title = {Dictionary learning for integrative, multimodal and scalable single-cell analysis},
 # journal = {Nature Biotechnology},
  #year = {2023},
 # doi = {10.1038/s41587-023-01767-y},
  #url = {https://doi.org/10.1038/s41587-023-01767-y},
#}
#dependencies and setup 
require(Seurat)
require(data.table)
library(dplyr)
library(R.utils)
setwd("/Users/radhikashaunak/Downloads/")
mat_healthy <- fread("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/healthy_epi.tsv.gz")

meta_healthy <- read.table("/Users/radhikashaunak/Desktop/Asthma_TFs/Metadata/meta_healthy_epi.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes_healthy = mat_healthy[,1][[1]]
genes_healthy = gsub(".+[|]", "", genes_healthy)
mat_healthy = data.frame(mat_healthy[,-1], row.names=genes_healthy)
#HEALTHY RNA-SEQ
so_healthy <- CreateSeuratObject(counts = mat_healthy, project = "epithelial_cells", meta.data=meta_healthy)
so_healthy[["percent.mt"]] <- PercentageFeatureSet(so_healthy, pattern = "^MT-")
so_healthy.filtered <- subset(so_healthy, subset= nCount_RNA>800 & 
                        nFeature_RNA >200 &
                        nFeature_RNA < 5000 &
                        percent.mt <5)

#QC visuals for healthy RNA-seq
VlnPlot(so_healthy.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
scatterplot1 <- FeatureScatter(so_healthy.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatterplot2 <- FeatureScatter(so_healthy.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatterplot1 + scatterplot2

#ASTHMA RNA-SEQ
mat_asthma <- fread("/Users/radhikashaunak/Desktop/Asthma_TFs/Originals/epithelial_asthma.tsv.gz")
meta_asthma <- read.table("/Users/radhikashaunak/Desktop/Asthma_TFs/Metadata/meta_asthma.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes_asthma = mat_asthma[,1][[1]]
genes_asthma = gsub(".+[|]", "", genes_asthma)
mat_asthma = data.frame(mat_asthma[,-1], row.names=genes_asthma)
so_asthma <- CreateSeuratObject(counts = mat_asthma, project = "epithelial_cells", meta.data=meta_asthma)
so_asthma[["percent.mt"]] <- PercentageFeatureSet(so_asthma, pattern = "^MT-")
so_asthma.filtered <- subset(so_asthma, subset= nCount_RNA>800 & 
                        nFeature_RNA >200 &
                        nFeature_RNA < 5000 &
                        percent.mt <5 &
                          Loc == "Asthma")

#QC visuals for asthma RNA-seq
VlnPlot(so_asthma.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
scatterplot1 <- FeatureScatter(so_asthma.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatterplot2 <- FeatureScatter(so_asthma.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatterplot1 + scatterplot2
# Adding condition labels to each object
so_healthy.filtered$condition <- "healthy"
so_asthma.filtered$condition <- "asthma"

#merging datasets
so_combined <- merge(
  so_healthy.filtered,
  y = so_asthma.filtered,
  add.cell.ids = c("healthy", "asthma"),  # labelling for cell distinguishing
  project = "combined_UCSC_project"
)
#Normalization
so_combined <- NormalizeData(so_combined, normalization.method = "LogNormalize", scale.factor = 10000)
#Feature Selection
so_combined <- FindVariableFeatures(so_combined, selection.method = "vst", nfeatures = 2000)
#top 10 most highly variable genes (diff exp)
top10.UCSC <- head(VariableFeatures(so_combined), 10)
#plotting variable features without labels to consider dispersion
unlabelled.plt <- VariableFeaturePlot(so_combined)
#plotting labelled to identify top 10 features
labelled.plt <- LabelPoints(plot = unlabelled.plt, points=top10.UCSC, repel=TRUE)
unlabelled.plt + labelled.plt
#scaling data
#all.genes <- rownames(so_combined)
so_combined <- ScaleData(so_combined)
#linear dimension reduction (PCA)
so_combined <- RunPCA(so_combined, features = VariableFeatures(object = so_combined))
VizDimLoadings(so_combined, dims = 1:2, reduction = "pca")
DimPlot(so_combined, reduction = "pca") + NoLegend()
DimHeatmap(so_combined, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(so_combined)
#cell clustering
so_combined <- FindNeighbors(so_combined, dims = 1:7)
so_combined <- FindClusters(so_combined, resolution = 0.15,
                            algorithm = 1)
#non-linear dimensional reduction (umap/tsne)
so_combined <- RunUMAP(so_combined, dims = 1:7)
DimPlot(so_combined, reduction = "umap")
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
so_combined <- JoinLayers(so_combined)
epithelial.markers <- FindAllMarkers(so_combined, only.pos = TRUE)
epithelial.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#Basal cell identification
VlnPlot(so_combined, features =c("TP63", "KRT5", "KRT15"))
DotPlot(so_combined, features =c("TP63", "KRT5", "KRT15"))
# club cell identification
VlnPlot(so_combined, features ="SCGB1A1")
DotPlot(so_combined, features ="SCGB1A1")
VlnPlot(so_combined, features = c("SCGB3A2"))
DotPlot(so_combined, features = c("SCGB3A2"))
#Goblet cell identification
VlnPlot(so_combined, features = c("BPIFB2"))
DotPlot(so_combined, features = c("BPIFB2"))
VlnPlot(so_combined, features = c("ATP6V1B1"))
DotPlot(so_combined, features = c("ATP6V1B1"))
DotPlot(so_combined, features = c("AGR2"))
DotPlot(so_combined, features = c("MUC5AC"))
#DC identification
VlnPlot(so_combined, features = c("PTPRC","CCDC50","IRF8"))
DotPlot(so_combined, features = c("PTPRC","CCDC50","IRF8"))
#Overlapping Markers
VlnPlot(so_combined, features = c("DNAH5","SPEF2", "PIFO","FOXJ1"))
DotPlot(so_combined, features = c("DNAH5","SPEF2", "PIFO","FOXJ1"))
VlnPlot(so_combined, features = c("DNAI1", "CCDC40", "RSPH1", "TEKT1"))
DotPlot(so_combined, features = c("DNAI1", "CCDC40", "RSPH1", "TEKT1"))
VlnPlot(so_combined, features = c("FOXJ1"))
VlnPlot(so_combined, features = c("TUBB4B","SNTN"))
DotPlot(so_combined, features = c("TUBB4B","SNTN"))
#Club/Goblet cell markers
VlnPlot(so_combined, features = c("MUC5AC","MUC5B","SPDEF","CFTR"))
DotPlot(so_combined, features = c("MUC5AC","MUC5B","SPDEF","CFTR"))
VlnPlot(so_combined, features = c("SPDEF"))
VlnPlot(so_combined, features = c("CYP2F1"))
DotPlot(so_combined, features = c("CYP2F1"))
#deuterosomal markers
VlnPlot(so_combined, features = c("CDC20B", "DEUP1", "FOXN4"))
DotPlot(so_combined, features = c("CDC20B", "DEUP1", "FOXN4"))
#mast cell markers
VlnPlot(so_combined, features = c("TPSAB1", "TPSB2", "HPGDS","KIT"))
DotPlot(so_combined, features = c("TPSAB1", "TPSB2", "HPGDS","KIT"))
#T-cell marker
VlnPlot(so_combined, features = c("CD2"))
DotPlot(so_combined, features = c("CD2"))
#Ionocyte cell markers
VlnPlot(so_combined, features = c("LINC01187","ATP6V1G3","FOXI1"))
DotPlot(so_combined, features = c("LINC01187","ATP6V1G3","FOXI1"))
VlnPlot(so_combined, features = c("TMPRSS11E","BSND","SFTPB"))
DotPlot(so_combined, features = c("TMPRSS11E","BSND","SFTPB"))
#Macrophage cell markers
VlnPlot(so_combined, features = c("CD68","MSR1","PPARG","APOE","CD14","CD163"))
DotPlot(so_combined, features = c("CD68","MSR1","PPARG","APOE","CD14","CD163"))
VlnPlot(so_combined, features = c("CSF1R"))
DotPlot(so_combined, features = c("CSF1R"))
#alveolar cell markers
VlnPlot(so_combined, features = c("GABRP","P2X7"))
DotPlot(so_combined, features = c("GABRP","P2X7"))
VlnPlot(so_combined, features = c("SFTPC","AGER"))
DotPlot(so_combined, features = c("SFTPC","AGER"))
VlnPlot(so_combined, features = c("SFTPB","PDPN"))
DotPlot(so_combined, features = c("SFTPB","PDPN"))
VlnPlot(so_combined, features = c("CEBPD"))
DotPlot(so_combined, features = c("CEBPD"))
VlnPlot(so_combined, features = c("POU2F1"))
DotPlot(so_combined, features = c("POU2F1"))
VlnPlot(so_combined, features = c("FOXP1"))
DotPlot(so_combined, features = c("FOXP1"))
VlnPlot(so_combined, features = c("GRAMD2"))
DotPlot(so_combined, features = c("GRAMD2"))
VlnPlot(so_combined, features = c("SCNN1G"))
DotPlot(so_combined, features = c("SCNN1G"))

#checking
markers <- FindAllMarkers(so_combined, only.pos = TRUE, min.pct = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(so_combined, features = top_markers$gene)

# cluster assignment
cell_types <- c("0" = "basal 1", 
                "1" = "ciliated 1",
                "2" = "goblet",
                "3" = "club",
                "4" = "basal 2",
                "5" = "mucociliated",
                "6" = "alveolar type 2",
                "7" = "ciliated 2",
                "8" = "ionocytes")
so_combined$cell_types <- recode(
  as.character(so_combined$seurat_clusters),
  !!!cell_types)

#setting identities for further analysis
Idents(so_combined) <- "cell_type"

 
# UMAP confirmation
DimPlot(so_combined, reduction = "umap", label = TRUE, group.by = "cell_types")
so_combined$cell_types
#comparison with identities in dataset
DimPlot(so_combined, reduction = "umap", label = TRUE, group.by = "Cluster")

#Differential expression using mast

#citation = MAST: a flexible statistical framework for assessing transcriptional changes and
#characterizing heterogeneity in single-cell RNA sequencing data G Finak, A McDavid,
#M Yajima, J Deng, V Gersuk, AK Shalek, CK Slichter et al Genome biology 16 (1), 278

#installing MAST
#install.packages("BiocManager")
#BiocManager::install("MAST")
#ensure that cell identities are unique
cell_identity <- unique(so_combined$cell_types)
#cell types identified in asthma and healthy individuals differed
#therefore need to split analysis into common cell type features and condition specific
#object for asthma
asthma_celltypes  <- unique(so_combined$cell_types[so_combined$condition == "asthma"])
#object for healthy
healthy_celltypes <- unique(so_combined$cell_types[so_combined$condition == "healthy"])
#common cell types across conditions
shared <- intersect(asthma_celltypes, healthy_celltypes)
#object for only asthma
only_asthma  <- setdiff(asthma_celltypes, healthy_celltypes)
#object for only healthy
only_healthy <- setdiff(healthy_celltypes, asthma_celltypes)

#empty list for differential expression results
de_results <- list()
#differential expression for shared cell types
#for each cell type in the column of the seurat object
for (ct in shared) {
  #extract a cell type-specific subset of the data
  subset_so_combined <- subset(so_combined, cell_types == ct)
  Idents(subset_so_combined) <- "condition"  # using "healthy" and "asthma" as identifying factors
  #for each cell type find disease status-specific markers (i.e. DEGs)
  de_results[[ct]] <- FindMarkers(
    subset_so_combined,
    ident.1 = "asthma",
    ident.2 = "healthy",
    test.use = "MAST",        # using MAST
    latent.vars = "nCount_RNA", # Correct for sequencing depth
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
}
#filtering for significant results
sig_basal1 <- de_results[["basal 1"]]%>% 
  filter(p_val_adj < 0.05)
sig_basal1$cell_type <- "Basal 1"
sig_basal1$gene <- rownames(sig_basal1)
sig_basal2 <- de_results[["basal 2"]]%>% 
  filter(p_val_adj < 0.05)
sig_basal2$cell_type <- "Basal 2"
sig_basal2$gene <- rownames(sig_basal2)
sig_ciliated1 <- de_results[["ciliated 1"]]%>% 
  filter(p_val_adj < 0.05)
sig_ciliated1$cell_type <- "ciliated 1"
sig_ciliated1$gene <- rownames(sig_ciliated1)
sig_ciliated2 <- de_results[["ciliated 2"]]%>%
  filter(p_val_adj < 0.05)
sig_ciliated2$cell_type <- "ciliated 2"
sig_ciliated2$gene <- rownames(sig_ciliated2)
sig_mucociliated <- de_results[["mucociliated"]]%>%
  filter(p_val_adj < 0.05)
sig_mucociliated$cell_type <- "mucociliated"
sig_mucociliated$gene <- rownames(sig_mucociliated)
sig_club <- de_results[["club"]]%>% 
  filter(p_val_adj < 0.05)
sig_club$cell_type <- "club"
sig_club$gene <- rownames(sig_club)
sig_goblet <- de_results[["goblet"]]%>% 
  filter(p_val_adj < 0.05)
sig_goblet$cell_type <- "goblet"
sig_goblet$gene <- rownames(sig_goblet)
sig_ionocytes <- de_results[["ionocytes"]]%>% 
  filter(p_val_adj < 0.05)
sig_ionocytes$cell_type <- "ionocytes"
sig_ionocytes$gene <- rownames(sig_ionocytes)
#combining into one dataframe
sig_results_total <- bind_rows(sig_basal1, sig_basal2, sig_ciliated1, sig_ciliated2,
                               sig_club, sig_goblet, sig_ionocytes, sig_mucociliated)
sig_results_total$gene <- rownames(sig_results_total)
head(sig_results_total)
#note: upregulation and downregulation in this dataframe refers to upregulation in asthma relative to healthy and vice versa
View(sig_basal1)
rownames(sig_basal1)
#transcription factor analysis for basal 1 cells
basal1_total <- so_combined$cell_types["basal 1"]
DotPlot(so_combined, features = c("LTF"), group.by = "cell_types")
DotPlot(so_combined, features = c("PLSCR1"), group.by = "cell_types")
DotPlot(so_combined, features = c("ZNF581","KLF1"), group.by = "cell_types")
DotPlot(so_combined, features = c("ZNF469", "NFE4"), group.by = "cell_types")
DotPlot(so_combined, features = c("ELF3","CEBPE"), group.by = "cell_types")
DotPlot(so_combined, features = c("LYL1","GATA1"), group.by = "cell_types")
# differential gene expression analysis for basal 1 cells
DotPlot(so_combined, features = c("POSTN"), group.by = "cell_types")
DotPlot(so_combined, features = c("ALOX15"), group.by = "cell_types")
DotPlot(so_combined, features = c("IGHG1", "HLA-DRB5"), group.by = "cell_types")
DotPlot(so_combined, features = c("IFITM1", "MMP1", "LYZ"), group.by = "cell_types")
DotPlot(so_combined, features = c("NME2", "HBB", "NOS2", "HBA1"), group.by = "cell_types")
DotPlot(so_combined, features = c("IGHG1", "HLA-DRB5"), group.by = "cell_types")
DotPlot(so_combined, features = c("IGHA1", "IGLC3", "ZNF667-AS1"), group.by = "cell_types")
DotPlot(so_combined, features = c("CDH26", "HBA2","PRB4"), group.by = "cell_types")
DotPlot(so_combined, features = c("MIR4435-2HG", "PLS1", "PRB3"), group.by = "cell_types")
DotPlot(so_combined, features = c("RPP21", "PP7080", "FAM110C"), group.by = "cell_types")
DotPlot(so_combined, features = c("IGHG3", "SERPINB10", "UBE2V1"), group.by = "cell_types")
DotPlot(so_combined, features = c("AC004556.1", "FAM184B"), group.by = "cell_types")
