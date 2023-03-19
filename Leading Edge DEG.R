library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(Rfast2)
library(MAST)
library(DESeq2)

install.packages("cli")
install.packages('Seurat')

#SWD
setwd("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/DEGs/Leading Edge DEG")


#Loading the data
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/DEGs/Leading Edge DEG/comb.Robj")

Combined_samples <- comb
DimPlot(Combined_samples, group.by = "cluster_annotations")
View(Combined_samples@meta.data)




Idents(Combined_samples) <- "sample_id"
View(Combined_samples)

#DEGs
#Differentially Expressed Gene

Metastasis_Present_ver2 <- subset(x = Combined_samples, idents = c("sample_1", "sample_2", "sample_4", "sample_9", "sample_10", "sample_12"))
Metastasis_Absent_ver2 <- subset(x = Combined_samples, idents = c("sample_3", "sample_5", "sample_6", "sample_7", "sample_8", "sample_11"))

DimPlot(Metastasis_Present_ver2, group.by = "cluster_annotations")
DimPlot(Metastasis_Absent_ver2, group.by = "cluster_annotations")

Idents(Metastasis_Present_ver2) <- "cluster_annotations"
Idents(Metastasis_Absent_ver2) <- "cluster_annotations"

View(Metastasis_Present_ver2)
View(Metastasis_Absent_ver2)

LE_subsetted_Present <- subset(x = Metastasis_Present_ver2, idents = "edge")
LE_subsetted_Absent <- subset(x = Metastasis_Absent_ver2, idents = "edge") 

DimPlot(LE_subsetted_Present, group.by = "cluster_annotations")
DimPlot(LE_subsetted_Absent, group.by = "cluster_annotations")


Idents(object = LE_subsetted_Present) <- "Lymph_Node_Metatasis"
Idents(object = LE_subsetted_Absent) <- "Non_Lymph_Node_Metatasis"

View(LE_subsetted_Present)
View(LE_subsetted_Absent)

DEG_merged <- merge (x = LE_subsetted_Present, y = LE_subsetted_Absent, project = "merged_LE_seurat")

View(DEG_merged$orig.ident)
View(DEG_merged@meta.data)
View(DEG_merged)


save(DEG_merged, file = "DEG_merged_LE.RData")

#Performing DEGs Expression)
## Ident 1  = Lymph-Metastasis / Ident 2 = Non- Lymph Metatasis / We are comparing Ident 1 vs Ident 2 
DefaultAssay(DEG_merged) <- "Spatial"

DEGsMarkers_LE_USING_WILCOX <- FindMarkers(DEG_merged, ident.1 = "Lymph_Node_Metatasis", ident.2 = "Non_Lymph_Node_Metatasis", verbose = FALSE)


#Filtering for DEGs from Default (Wilcox)
p_less0.01 <- subset(DEGsMarkers_LE_USING_WILCOX, p_val_adj < 0.01)
Upregulated <- subset(p_less0.01,  avg_log2FC > 0.58 )
Downregulated <- subset(p_less0.01,  avg_log2FC < -0.58)

CombinedDEG_LE <- combine(Upregulated, Downregulated)



write.csv(DEGsMarkers_LE_USING_WILCOX, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/DEGsMarkers_LE_USING_WILCOX.csv")
write.csv(Upregulated, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/Upregulated.DEG_LE.csv")
write.csv(Downregulated, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/Downregulated_LE.csv")



