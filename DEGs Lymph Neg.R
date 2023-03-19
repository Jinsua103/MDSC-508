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

Lymph_Neg_Subsetted_Metastasis_samples <- subset(x = Metastasis_samples, idents = "Lymphocyte Negative Stroma")
Lymph_Neg_Subsetted_Non_Metastasis_samples <- subset(x = Non_Metastasis_samples, idents = "Lymphocyte Negative Stroma")


View(SCC_Subsetted_Metastasis_samples)
View(SCC_Subsetted_Metastasis_samples@meta.data)

View(SCC_Subsetted_Non_Metastasis_samples)
View(SCC_Subsetted_Non_Metastasis_samples@meta.data)


#Re-assigning the identity of Non-metastasis and Metastasis for the merged data
Idents(object = Lymph_Neg_Subsetted_Metastasis_samples) <- "Metatasis"
Idents(object = Lymph_Neg_Subsetted_Non_Metastasis_samples) <- "Non-Metatasis"


#Final Merged File with all SCC, but with annotated Non-Metastsis group and Metatasis group
Merged_LymNeg <- merge(x = Lymph_Neg_Subsetted_Metastasis_samples, y = Lymph_Neg_Subsetted_Non_Metastasis_samples, project = "merged_LymNeg_seurat")

View(Merged_LymNeg@meta.data$orig.ident)
View(Merged_LymNeg@meta.data)
View(Merged_LymNeg)

##
#WhichCells(Merged, idents = c("benign", "malignant"))


save(Merged, file = "Merged_AllSamplesCombined.RData")


#Performing DEGs Expression)
## Ident 1  = Metastasis / Ident 2 = Non- Metatasis / We are comparing Ident 1 vs Ident 2 

#DEGs Using Default (Wilcox)
DEGs_Lymph_NEG_POSITIVEvsNEGATIVE <- FindMarkers(Merged_LymNeg, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#FINAL Fring for DEGs
p_less0.01_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE <- subset(DEGs_Lymph_NEG_POSITIVEvsNEGATIVE, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE <- subset(p_less0.01_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE <- subset(p_less0.01_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE,  avg_log2FC < -0.58)
Filtered_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE <- combine(p_adj_greaterthan0.58_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE, p_adj_lessthan0.58_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE)

write.csv(Filtered_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Filtered_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE.csv")





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

Lymph_Pos_Subsetted_Metastasis_samples <- subset(x = Metastasis_samples, idents = "Lymphocyte Positive Stroma")
Lymph_Pos_Subsetted_Non_Metastasis_samples <- subset(x = Non_Metastasis_samples, idents = "Lymphocyte Positive Stroma")




#Re-assigning the identity of Non-metastasis and Metastasis for the merged data
Idents(object = Lymph_Pos_Subsetted_Metastasis_samples) <- "Metatasis"
Idents(object = Lymph_Pos_Subsetted_Non_Metastasis_samples) <- "Non-Metatasis"


#Final Merged File with all SCC, but with annotated Non-Metastsis group and Metatasis group
Merged_LymPos <- merge(x = Lymph_Pos_Subsetted_Metastasis_samples, y = Lymph_Pos_Subsetted_Non_Metastasis_samples, project = "merged_LymPos_seurat")

View(Merged_LymNeg@meta.data$orig.ident)
View(Merged_LymNeg@meta.data)
View(Merged_LymNeg)

##
#WhichCells(Merged, idents = c("benign", "malignant"))


save(Merged, file = "Merged_AllSamplesCombined.RData")


#Performing DEGs Expression)
## Ident 1  = Metastasis / Ident 2 = Non- Metatasis / We are comparing Ident 1 vs Ident 2 

#DEGs Using Default (Wilcox)
DEGs_Lymph_POS_POSITIVEvsNEGATIVE <- FindMarkers(Merged_LymPos, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#FINAL Fring for DEGs
p_less0.01_DEGs_Lymph_POS_POSITIVEvsNEGATIVE <- subset(DEGs_Lymph_POS_POSITIVEvsNEGATIVE, p_val_adj < 0.01)
p_adj_greaterthan0.58_DEGs_Lymph_POS_POSITIVEvsNEGATIVE <- subset(p_less0.01_DEGs_Lymph_POS_POSITIVEvsNEGATIVE,  avg_log2FC > 0.58 )
p_adj_lessthan0.58_DEGs_Lymph_POS_POSITIVEvsNEGATIVE <- subset(p_less0.01_DEGs_Lymph_POS_POSITIVEvsNEGATIVE,  avg_log2FC < -0.58)
Filtered_DEGs_Lymph_POS_POSITIVEvsNEGATIVE <- combine(p_adj_greaterthan0.58_DEGs_Lymph_POS_POSITIVEvsNEGATIVE, p_adj_lessthan0.58_DEGs_Lymph_POS_POSITIVEvsNEGATIVE)

write.csv(Filtered_DEGs_Lymph_POS_POSITIVEvsNEGATIVE, 
          file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/Research/Research with Dr. Bose/Filtered_DEGs_Lymph_POS_POSITIVEvsNEGATIVE.csv")
