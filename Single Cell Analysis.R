
  
#Single-Cell mRNAseq Analysis - HNSC data - Puram et al., 2017
#Data download from NCBI GEO - Single cell transcriptomes and metadata (sample/ cell type):
  
#Import HNSc scRNAseq data - Use read.delim2 to import cell type annotations as "strings", does not show up as NA if use read_tsv

setwd("C:/Users/jinsu/OneDrive/Desktop/Jinsu//Research/Research with Dr. Bose/Dataset/single cell reference")

install.packages("Seurat")
install.packages("SeuratObject")
library(Seurat)
library(gdata)

all_data <- read.delim2("HNSCC_all_data.txt")
#Clean Gene names - drop '' around gene name

colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN23_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC20_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC12_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC10_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC6_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC5_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC8_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC7_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC13_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC16_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC18_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC22_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_24_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC24_", "", x)))




colnames(all_data)





colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P5", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P6", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P7", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P8", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P9", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P10", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P12", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P13", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P16", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P17", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P18", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P20", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P22", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P23", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P24", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P25", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P26", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P28", "LNM_YES_", x)))

colnames(all_data)
options(max.print = 5000)

all_data$X  <- as.list(sapply(all_data$X , function(x) gsub("\'", "", x)))
rownames(all_data) <- all_data$X
all_data$X <- NULL


index_TRUE <- startsWith(colnames(all_data), "LNM_YES_")
index_FALSE <- startsWith(colnames(all_data), "LNM_NO_")
LNM_PRESENT <-subset(all_data, select = ((startsWith(colnames(all_data), "LNM_YES_")== TRUE )))
LNM_ABSENT <-subset(all_data, select = ((startsWith(colnames(all_data), "LNM_NO_")== TRUE)))

#Create rna_meta data using the first 5 rows of the dataframe
rna_meta_PRESENT <- LNM_PRESENT[1:5,]
rna_meta_ABSENT <- LNM_ABSENT[1:5,]

rna_meta_PRESENT <- as.data.frame(t(rna_meta_PRESENT))
rna_meta_ABSENT <- as.data.frame(t(rna_meta_ABSENT))
rownames(rna_meta_PRESENT) <- gsub("_", "-", rownames(rna_meta_PRESENT))
rownames(rna_meta_ABSENT) <- gsub("_", "-", rownames(rna_meta_ABSENT))

index1 <- rna_meta_PRESENT$`non-cancer cell type` == 0
index2 <- rna_meta_ABSENT$`non-cancer cell type` == 0

rna_meta_PRESENT$`non-cancer cell type`[index1] <- "Cancer cell"
rna_meta_ABSENT$`non-cancer cell type`[index2] <- "Cancer cell"

index3 <- rna_meta_PRESENT$`non-cancer cell type` == "-Fibroblast"
index4 <- rna_meta_ABSENT$`non-cancer cell type` == "-Fibroblast"

rna_meta_PRESENT$`non-cancer cell type`[index3] <- "Fibroblast"
rna_meta_ABSENT$`non-cancer cell type`[index4] <- "Fibroblast"


#Create the counts matrix using the remaining data
rna_clean_PRESENT <- LNM_PRESENT[6:23691,]
rna_clean_ABSENT <- LNM_ABSENT[6:23691,]

names(rna_clean_PRESENT) <- gsub("_", "-", names(rna_clean_PRESENT))
names(rna_clean_ABSENT) <- gsub("_", "-", names(rna_clean_ABSENT))

#create seurat objects used for downstream analysis 
rna_seurat_LNM_POSITIVE <- CreateSeuratObject(rna_clean_PRESENT, project = "puram_data_LNM_POSITIVE")
rna_seurat_LNM_NEGATIVE <- CreateSeuratObject(rna_clean_ABSENT, project = "puram_data_LNM_ABSENT")

rna_seurat_LNM_POSITIVE <- AddMetaData(rna_seurat_LNM_POSITIVE, metadata = rna_meta_PRESENT)
rna_seurat_LNM_NEGATIVE <- AddMetaData(rna_seurat_LNM_NEGATIVE, metadata = rna_meta_ABSENT)


#Perform basic data validation - examine basic characteristics of data
counts_per_cell_P <- Matrix::colSums(rna_seurat_LNM_POSITIVE)
counts_per_gene_P <- Matrix::rowSums(rna_seurat_LNM_POSITIVE)
genes_per_cell_P <- Matrix::colSums(rna_seurat_LNM_POSITIVE)
cells_per_gene_P <- Matrix::rowSums(rna_seurat_LNM_POSITIVE)

counts_per_cell_A <- Matrix::colSums(rna_seurat_LNM_NEGATIVE)
counts_per_gene_A <- Matrix::rowSums(rna_seurat_LNM_NEGATIVE)
genes_per_cell_A <- Matrix::colSums(rna_seurat_LNM_NEGATIVE)
cells_per_gene_A <- Matrix::rowSums(rna_seurat_LNM_NEGATIVE)


#Graph the counts per cell, counts per gene, cells per gene, etc 
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')



#Normalize the RNA data in the seurat object:
rna_seurat_LNM_POSITIVE <- NormalizeData(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- NormalizeData(rna_seurat_LNM_NEGATIVE)

rna_seurat_LNM_POSITIVE <- FindVariableFeatures(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- FindVariableFeatures(rna_seurat_LNM_NEGATIVE)


top10_P <- head(VariableFeatures(rna_seurat_LNM_POSITIVE), 10)
top10_A <- head(VariableFeatures(rna_seurat_LNM_NEGATIVE), 10)

plot1_P <- VariableFeaturePlot(rna_seurat_LNM_POSITIVE)
LabelPoints(plot = plot1_P, points = top10_P, repel = TRUE, xnudge = 0, ynudge = 0)

plot1_A <- VariableFeaturePlot(rna_seurat_LNM_NEGATIVE)
LabelPoints(plot = plot1_A, points = top10_A, repel = TRUE, xnudge = 0, ynudge = 0)


rna_seurat_LNM_POSITIVE <- ScaleData(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- ScaleData(rna_seurat_LNM_NEGATIVE)


rna_seurat_LNM_POSITIVE <- RunPCA(rna_seurat_LNM_POSITIVE, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- RunPCA(rna_seurat_LNM_NEGATIVE, verbose = FALSE)


DimHeatmap(rna_seurat_LNM_POSITIVE, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat_LNM_POSITIVE)

DimHeatmap(rna_seurat_LNM_NEGATIVE, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat_LNM_NEGATIVE)

rna_seurat_LNM_POSITIVE <- FindNeighbors(rna_seurat_LNM_POSITIVE, dims = 1:10)
rna_seurat_LNM_NEGATIVE <- FindNeighbors(rna_seurat_LNM_NEGATIVE, dims = 1:10)

rna_seurat_LNM_POSITIVE <- FindClusters(rna_seurat_LNM_POSITIVE, resolution = 0.8, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- FindClusters(rna_seurat_LNM_NEGATIVE, resolution = 0.8, verbose = FALSE)

rna_seurat_LNM_POSITIVE <- RunUMAP(rna_seurat_LNM_POSITIVE, dims = 1:10, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- RunUMAP(rna_seurat_LNM_NEGATIVE, dims = 1:10, verbose = FALSE)

DimPlot(rna_seurat_LNM_POSITIVE, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE, DarkTheme()) 
DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE) 

rna_seurat_LNM_POSITIVE <- RunTSNE(rna_seurat_LNM_POSITIVE, dims = 1:10)
rna_seurat_LNM_NEGATIVE <- RunTSNE(rna_seurat_LNM_NEGATIVE, dims = 1:10)

plot_for_LNM_POS <- DimPlot(rna_seurat_LNM_POSITIVE, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)
plot_for_LNM_NEG <- DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)

#DEGs for Cancer cell type
Idents(rna_seurat_LNM_POSITIVE) <- 'non.cancer.cell.type'
Idents(rna_seurat_LNM_NEGATIVE) <- 'non.cancer.cell.type'

Cancer_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Cancer cell")
Cancer_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Cancer cell")

Idents(object = Cancer_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Cancer_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_Cancer_celltype <- merge(x = Cancer_Subsetted_LNM_PRESENT, y = Cancer_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_inCancer_LNMPOS <- FindMarkers(Merged_Cancer_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for T cell type
T_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "T cell")
T_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "T cell")


Idents(object = T_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = T_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_T_cell_celltype <- merge(x = T_Subsetted_LNM_PRESENT, y = T_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Tcells_LNMPOS <- FindMarkers(Merged_T_cell_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for Fibroblast cell type
Fibroblast_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Fibroblast")
Fibroblast_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Fibroblast")

Idents(object = Fibroblast_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Fibroblast_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Fibroblast_celltype <- merge(x = Fibroblast_Subsetted_LNM_PRESENT, y = Fibroblast_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Fibroblast_LNMPOS <- FindMarkers(Merged_Fibroblast_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

View(DEGs_inCancer_LNMPOS)
View(DEGs_in_Tcells_LNMPOS)
View(DEGs_in_Fibroblast_LNMPOS)

p_val_Cancer <- subset(DEGs_inCancer_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC > 0.58 )
logFold_lessthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC < -0.58)
Filtered_Cancer <- combine(logFold_greaterthan0.58_cancer, logFold_lessthan0.58_cancer)

p_val_T <- subset(DEGs_in_Tcells_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_T <- subset(p_val_T,  avg_log2FC > 0.58 )
logFold_lessthan0.58_T <- subset(p_val_T,  avg_log2FC < -0.58)
Filtered_T <- combine(logFold_greaterthan0.58_T, logFold_lessthan0.58_T)

p_val_Fibroblast <- subset(DEGs_in_Fibroblast_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC > 0.58 )
logFold_lessthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC < -0.58)
Filtered_Fibroblast <- combine(logFold_greaterthan0.58_Fibroblast, logFold_lessthan0.58_Fibroblast)
















#####
table(rna_seurat@meta.data$seurat_clusters)
DimPlot(rna_seurat, label.size = 4, repel = TRUE, label = TRUE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

rna_seurat <- CellCycleScoring(rna_seurat, s.features = s.genes, g2m.features = g2m.genes)
table(rna_seurat[[]]$Phase)

#FeaturePlot(rna_seurat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE) 

rna_seurat <- RunTSNE(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)

FeatureScatter(object = rna_seurat, feature1 = "RNF126", feature2 = "PC_1")

VlnPlot(object = rna_seurat, features = c("MIER2"))
