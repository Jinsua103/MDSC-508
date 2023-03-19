library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(Rfast2)
library(MAST)
library(DESeq2)

install.packages('Seurat')
install.packages('MAST')
install.packages('DESeq2')
install.packages("BiocManager") # Needed to install all Bioconductor packages
BiocManager::install("MAST")



 
#Loading the data
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_1.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_2.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_3.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_4.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_5.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_6.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_7.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_8.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_9.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_10.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_11.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_12.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_13.Robj")
load("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/obj_14.Robj")



#Assigning 
Sample_1 <- obj_1 #Patient with metastasis: T2N2bM1, lymph node and distant metastasis
Sample_1 <- AddMetaData(Sample_1, metadata = "Yes", col.name = "Metastasis")

Sample_2 <- obj_2 #Patient with metastasis: T4aN1M1, lymph node and distant metastasis
Sample_2 <- AddMetaData(Sample_2, metadata = "Yes", col.name = "Metastasis")

Sample_3 <- obj_3 #Patient without metastasis: T1N0M0
Sample_3 <- AddMetaData(Sample_3, metadata = "No", col.name = "Metastasis")

Sample_4 <- obj_4 #Patient with metastasis: T2N2bM1, lymph node and distant metastasis
Sample_4 <- AddMetaData(Sample_4, metadata = "Yes", col.name = "Metastasis")

Sample_5 <-obj_5 #Patient without metastasis: T2N0M0
Sample_5 <- AddMetaData(Sample_5, metadata = "No", col.name = "Metastasis")

Sample_6 <-obj_6 #Patient without metastasis: T2N0M0
Sample_6 <- AddMetaData(Sample_6, metadata = "No", col.name = "Metastasis")

Sample_7 <-obj_7 #Patient without metastasis: T4aN0M1
Sample_7 <- AddMetaData(Sample_7, metadata = "No", col.name = "Metastasis")

Sample_8 <-obj_8 #Patient without metastasis: T4aN0M1
Sample_8 <- AddMetaData(Sample_8, metadata = "No", col.name = "Metastasis")

Sample_9 <-obj_9 #Patient with metastasis: T4aN1M1, lymph node and distant metastasis
Sample_9 <- AddMetaData(Sample_9, metadata = "Yes", col.name = "Metastasis")

Sample_10 <-obj_10 #Patient with metastasis: T4aN1M1, lymph node and distant metastasis
Sample_10 <- AddMetaData(Sample_10, metadata = "Yes", col.name = "Metastasis")

Sample_11 <-obj_11 #Patient without metastasis: T2N0M0
Sample_11 <- AddMetaData(Sample_11, metadata = "No", col.name = "Metastasis")

Sample_12 <-obj_12 #Patient with metastasis: T4aN2bM, lymph node and distant metastasis
Sample_12 <- AddMetaData(Sample_12, metadata = "Yes", col.name = "Metastasis")

#DEGs
#Differentially Expressed Gene

##Merging the samples based on the metastasis
Metastasis_samples <- merge(Sample_1, y = list(Sample_2, Sample_4, Sample_9, Sample_10, Sample_12))
Non_Metastasis_samples <- merge(Sample_3, y = list(Sample_5, Sample_6, Sample_7, Sample_8, Sample_11))


View(Metastasis_samples)
View(Non_Metastasis_samples)

#Idents are now classifed based on pathologist annotation for Metastasis sample
Idents(Metastasis_samples) <- 'pathologist_anno'
View(Metastasis_samples)

#Idents are now classifed based on pathologist annotation for Non-Metastasis sample
Idents(Non_Metastasis_samples) <- 'pathologist_anno'
View(Non_Metastasis_samples)


#WhichCell Funtions (CHECLKING IF THIS IS NECESSASRY TO HAVE IT FOR NOW)
##Metastasis_samples <- WhichCells(Metastasis_samples, expression = pathologist_anno == 'SCC')
##Non_Metastasis_samples <- WhichCells(Non_Metastasis_samples, expression = pathologist_anno == 'SCC')

View(Metastasis_samples)
View(Non_Metastasis_samples)

#1st Subsetting based on "SCC" //This is to filter everything else but SCC

SCC_Subsetted_Metastasis_samples <- subset(x = Metastasis_samples, idents = "SCC")
SCC_Subsetted_Non_Metastasis_samples <- subset(x = Non_Metastasis_samples, idents = "SCC")


View(SCC_Subsetted_Metastasis_samples)
View(SCC_Subsetted_Metastasis_samples@meta.data)

View(SCC_Subsetted_Non_Metastasis_samples)
View(SCC_Subsetted_Non_Metastasis_samples@meta.data)


#Re-assigning the identity of Non-metastasis and Metastasis for the merged data
Idents(object = SCC_Subsetted_Metastasis_samples) <- "Metatasis"
Idents(object = SCC_Subsetted_Non_Metastasis_samples) <- "Non-Metatasis"


#Final Merged File with all SCC, but with annotated Non-Metastsis group and Metatasis group
Merged <- merge(x = SCC_Subsetted_Metastasis_samples, y = SCC_Subsetted_Non_Metastasis_samples, project = "merged_SCC_seurat")

View(Merged@meta.data$orig.ident)
View(Merged@meta.data)
View(Merged)

##
#WhichCells(Merged, idents = c("benign", "malignant"))


save(Merged, file = "Merged_AllSamplesCombined.RData")


#Performing DEGs Expression)
## Ident 1  = Metastasis / Ident 2 = Non- Metatasis / We are comparing Ident 1 vs Ident 2 

#DEGs Using Default (Wilcox)
DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX <- FindMarkers(Merged, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#FINAL PRINTOUT
DEGs_LNM_POSITIVE <- FindMarkers(Merged, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)
DEGs_LNM_NEGATIVE <- FindMarkers(Merged, ident.1 = "Non-Metatasis", ident.2 = "Metatasis", recorrect_umi = FALSE)



#DEGS using MAST 
DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST <- FindMarkers(Merged, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", test.use = "MAST", assay = "Spatial")

#DEGS using DESeq2 ##DOESNOT WORK
DEGsMarkers_Metastasis_vs_Non_metastasis_USING_DESeq2 <- FindMarkers(Merged, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", test.use = "DESeq2", assay = "Spatial")

##ERROR:Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
##every gene contains at least one zero, cannot compute log geometric means

###ERROR-SOLVING  ##DOESNOT WORK
Merged[["Spatial"]]@counts<-as.matrix(Merged[["Spatial"]]@counts)+1
DEGsMarkers_Metastasis_vs_Non_metastasis_USING_DESeq2 <- FindMarkers(Merged, ident.1 = "Metatasis", ident.2 = "Non-Metatasis",  test.use = "DESeq2", assay = "Spatial", slot = "counts")











#Filtering for DEGs from Default (Wilcox)
p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX <- subset(DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX <- subset(p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX <- subset(p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX,  avg_log2FC < -0.58)
Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX <- combine(p_adj_greaterthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX, p_adj_lessthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX)
list(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX)

#FINAL Fring for DEGs
p_less0.01_DEGs_LNM_POSITIVE <- subset(DEGs_LNM_POSITIVE, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGs_LNM_POSITIVE <- subset(p_less0.01_DEGs_LNM_POSITIVE,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGs_LNM_POSITIVE <- subset(p_less0.01_DEGs_LNM_POSITIVE,  avg_log2FC < -0.58)
Filtered_DEGs_LNM_POSITIVE <- combine(p_adj_greaterthan0.58_DEGs_LNM_POSITIVE, p_adj_lessthan0.58_DEGs_LNM_POSITIVE)

write.csv(Filtered_DEGs_LNM_POSITIVE, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Filtered_DEGs_LNM_POSITIVE.csv")

p_less0.01_DEGs_LNM_NEGATIVE <- subset(DEGs_LNM_NEGATIVE, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGs_LNM_NEGATIVE <- subset(p_less0.01_DEGs_LNM_NEGATIVE,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGs_LNM_NEGATIVE <- subset(p_less0.01_DEGs_LNM_NEGATIVE,  avg_log2FC < -0.58)
Filtered_DEGs_LNM_NEGATIVE <- combine(p_adj_greaterthan0.58_DEGs_LNM_NEGATIVE, p_adj_lessthan0.58_DEGs_LNM_NEGATIVE)
write.csv(Filtered_DEGs_LNM_NEGATIVE, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Filtered_DEGs_LNM_NEGATIVE.csv")





write.csv(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/Samples_reajusted_Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX.csv")

#Filtering for DEGs from MAST
p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST <- subset(DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGsMarker_vs_Metastasis_vs_Non_metastasis_USING_MAST <- subset(p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST <- subset(p_less0.01_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST,  avg_log2FC < -0.58)
Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST <- combine(p_adj_greaterthan0.58_DEGsMarker_vs_Metastasis_vs_Non_metastasis_USING_MAST, p_adj_lessthan0.58_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST)
list(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST)


write.csv(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Dataset/Samples_reajusted_Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_MAST.csv")


#Filtering for DEGs from DESeq2

##N/A FOR NOW DUE TO ERROR IN DESEQ2 


SpatialFeaturePlot(Sample_1, features = )



#markers_malignant_vs_malignant <- FindMarkers(merged, ident.1 = "malignant", ident.2 = "benign", recorrect_umi = FALSE)

write.csv(markers_benign_vs_malignant, 
          file = "~/Desktop/MDSC 508 Research/R obj sample data/markers_benign_vs_malignant_before_p_value_0.01adj.csv")
#write.csv(markers_malignant_vs_malignant, 
#file = "~/Desktop/MDSC 508 Research/R obj sample data/markers_malignant_vs_malignant_before_p_value_0.01adj.csv")

count(markers_benign_vs_malignant)

markers_benign_vs_malignant_Method_MAST_filtering <- subset(markers_benign_vs_malignant_Method_MAST, p_val_adj < 0.01)
markers_benign_vs_malignant <- subset(markers_benign_vs_malignant, p_val_adj < 0.01)
#markers_malignant_vs_malignant <- subset(markers_malignant_vs_malignant, p_val_adj < 0.01)
count(markers_benign_vs_malignant)
count(markers_benign_vs_malignant_Method_MAST)
count(markers_benign_vs_malignant_Method_MAST_filtering)
count(markers_benign_vs_malignant)


View(markers_benign_vs_malignant)
#View(markers_malignant_vs_malignant)

p_adj_greaterthan0.58_markers_benign_vs_malignant_Method_MAST <- subset(markers_benign_vs_malignant_Method_MAST_filtering,  avg_log2FC > 0.58 )
p_adj_lessthan0.58markers_markers_benign_vs_malignant_Method_MAST <- subset(markers_benign_vs_malignant_Method_MAST_filtering,  avg_log2FC < -0.58)
p_adj_markers_benign_vs_malignant_Method_MAST <- combine(p_adj_greaterthan0.58_markers_benign_vs_malignant_Method_MAST, p_adj_lessthan0.58markers_markers_benign_vs_malignant_Method_MAST)
count(p_adj_markers_benign_vs_malignant_Method_MAST)



write.csv(p_adj_markers_benign_vs_malignant_Method_MAST, 
          file = "~/Desktop/MDSC 508 Research/R obj sample data/compare.csv")

p_adj_greaterthan0.4_markers_benign_vs_malignant <- subset(markers_benign_vs_malignant,  avg_log2FC > 0.58 )
count(p_adj_greaterthan0.4_markers_benign_vs_malignant)

p_adj_lessthan0.4markers_benign_vs_malignant <- subset(markers_benign_vs_malignant,  avg_log2FC < -0.58)
count(p_adj_lessthan0.4markers_benign_vs_malignant)

p_adj_markers_benign_vs_malignant <- combine(p_adj_greaterthan0.4_markers_benign_vs_malignant, p_adj_lessthan0.4markers_benign_vs_malignant)

write.csv(p_adj_markers_benign_vs_malignant, 
          file = "~/Desktop/MDSC 508 Research/R obj sample data/DEGs betweeb conditions p-value and avg_log2FC adjusted .csv")
#write.csv(markers_malignant_vs_malignant, 
#file = "~/Desktop/MDSC 508 Research/R obj sample data/markers_malignant_vs_malignant_after_p_value_0.01adj.csv")

p_adj_markers_benign_vs_malignant_Method_MAST
