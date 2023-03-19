if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)


DEGs_LNM_POSITIVE
DEGs_LNM_NEGATIVE




DEGs_LNM_POSITIVE_NO_FILTER <- EnhancedVolcano(DEGs_LNM_POSITIVE, 
                                               lab = rownames(DEGs_LNM_POSITIVE),
                                               x ="avg_log2FC", 
                                               y ="p_val_adj",
                                               title = 'DEGs in SCC of LNM+ group before Filtering ',
                                               pCutoff = 0.01,
                                               FCcutoff = 0.58,
                                               pointSize = 1.5,
                                               labSize = 4.0)

DEGs_LNM_POSTIVE_AFTER_FILTER <- EnhancedVolcano(Filtered_DEGs_LNM_POSITIVE , 
                                               lab = rownames(Filtered_DEGs_LNM_POSITIVE),
                                               x ="avg_log2FC", 
                                               y ="p_val_adj",
                                               title = 'DEGs in SCC of LNM+ group after Filtering ',
                                               pCutoff = 0.01,
                                               FCcutoff = 0.58,
                                               pointSize = 1.5,
                                               labSize = 4.0)



Fig_3_1 <- ggarrange(DEGs_LNM_POSITIVE_NO_FILTER, DEGs_LNM_POSTIVE_AFTER_FILTER + rremove("x.text"), 
                     
                     ncol = 2, nrow = 1)
 


DEGs_Lymph_NEG_POSITIVEvsNEGATIVE_NO_FILTER <- EnhancedVolcano(DEGs_Lymph_NEG_POSITIVEvsNEGATIVE , 
                                               lab = rownames(DEGs_Lymph_NEG_POSITIVEvsNEGATIVE),
                                               x ="avg_log2FC", 
                                               y ="p_val_adj",
                                               title = 'DEGs in Lymphocyte Negative Stroma of LNM+ before Filtering ',
                                               pCutoff = 0.01,
                                               FCcutoff = 0.58,
                                               pointSize = 1.5,
                                               labSize = 4.0)

DEGs_Lymph_NEG_POSITIVEvsNEGATIVE_AFTER_FILTER <- EnhancedVolcano(Filtered_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE , 
                                                  lab = rownames(Filtered_DEGs_Lymph_NEG_POSITIVEvsNEGATIVE),
                                                  x ="avg_log2FC", 
                                                  y ="p_val_adj",
                                                  title = 'DEGs in Lymphocyte Negative Stroma of LNM+ after Filtering ',
                                                  pCutoff = 0.01,
                                                  FCcutoff = 0.58,
                                                  pointSize = 1.5,
                                                  labSize = 4.0)



DEGs_Lymph_POS_POSITIVEvsNEGATIVE_NO_FILTER <- EnhancedVolcano(DEGs_Lymph_POS_POSITIVEvsNEGATIVE , 
                                                               lab = rownames(DEGs_Lymph_POS_POSITIVEvsNEGATIVE),
                                                               x ="avg_log2FC", 
                                                               y ="p_val_adj",
                                                               title = 'DEGs in Lymphocyte Positive Stroma of LNM+ before Filtering ',
                                                               pCutoff = 0.01,
                                                               FCcutoff = 0.58,
                                                               pointSize = 1.5,
                                                               labSize = 4.0)

DEGs_Lymph_POS_POSITIVEvsNEGATIVE_AFTER_FILTER <- EnhancedVolcano(Filtered_DEGs_Lymph_POS_POSITIVEvsNEGATIVE , 
                                                                  lab = rownames(Filtered_DEGs_Lymph_POS_POSITIVEvsNEGATIVE),
                                                                  x ="avg_log2FC", 
                                                                  y ="p_val_adj",
                                                                  title = 'DEGs in Lymphocyte Positive Stroma of LNM+ after Filtering ',
                                                                  pCutoff = 0.01,
                                                                  FCcutoff = 0.58,
                                                                  pointSize = 1.5,
                                                                  labSize = 4.0)



EnhancedVolcano(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX , 
                lab = rownames(Filtered_DEGsMarkers_Metastasis_vs_Non_metastasis_USING_WILCOX),
                x ="avg_log2FC", 
                y ="p_val_adj",
                title = 'DEGs in SCC Between Meta+ vs Meta - after Filtering ',
                pCutoff = 0.01,
                FCcutoff = 0.58,
                pointSize = 1.5,
                labSize = 4.0)

EnhancedVolcano(DEGsMarkers_LE_USING_WILCOX , 
                lab = rownames(DEGsMarkers_LE_USING_WILCOX),
                x ="avg_log2FC", 
                y ="p_val_adj",
                title = 'DEGs in LE Between Meta+ vs Meta - before Filtering ',
                pCutoff = 0.01,
                FCcutoff = 0.58,
                pointSize = 1.5,
                labSize = 4.0)


EnhancedVolcano(CombinedDEG_LE , 
                lab = rownames(CombinedDEG_LE),
                x ="avg_log2FC", 
                y ="p_val_adj",
                title = 'DEGs in LE Between Meta+ vs Meta - after Filtering ',
                pCutoff = 0.01,
                FCcutoff = 0.58,
                pointSize = 1.5,
                labSize = 4.0)



Figure1_Sample1 <- SpatialPlot(Sample_1, group.by = "pathologist_anno")
Figure1_Sample2 <- SpatialPlot(Sample_2, group.by = "pathologist_anno")
Figure1_Sample3 <- SpatialPlot(Sample_3, group.by = "pathologist_anno")
Figure1_Sample4 <- SpatialPlot(Sample_4, group.by = "pathologist_anno")
Figure1_Sample5 <- SpatialPlot(Sample_5, group.by = "pathologist_anno")
Figure1_Sample6 <- SpatialPlot(Sample_6, group.by = "pathologist_anno")
Figure1_Sample7 <- SpatialPlot(Sample_7, group.by = "pathologist_anno")
Figure1_Sample8 <- SpatialPlot(Sample_8, group.by = "pathologist_anno")
Figure1_Sample9 <- SpatialPlot(Sample_9, group.by = "pathologist_anno")
Figure1_Sample10 <- SpatialPlot(Sample_10, group.by = "pathologist_anno")
Figure1_Sample11 <- SpatialPlot(Sample_11, group.by = "pathologist_anno")
Figure1_Sample12 <- SpatialPlot(Sample_12, group.by = "pathologist_anno")


Figure2_Sample1 <- SpatialPlot(Sample_1, group.by = "core_edge_anno")
Figure2_Sample2 <- SpatialPlot(Sample_2, group.by = "core_edge_anno")
Figure2_Sample3 <- SpatialPlot(Sample_3, group.by = "core_edge_anno")
Figure2_Sample4 <- SpatialPlot(Sample_4, group.by = "core_edge_anno")
Figure2_Sample5 <- SpatialPlot(Sample_5, group.by = "core_edge_anno")
Figure2_Sample6 <- SpatialPlot(Sample_6, group.by = "core_edge_anno")
Figure2_Sample7 <- SpatialPlot(Sample_7, group.by = "core_edge_anno")
Figure2_Sample8 <- SpatialPlot(Sample_8, group.by = "core_edge_anno")
Figure2_Sample9 <- SpatialPlot(Sample_9, group.by = "core_edge_anno")
Figure2_Sample10 <- SpatialPlot(Sample_10, group.by = "core_edge_anno")
Figure2_Sample11 <- SpatialPlot(Sample_11, group.by = "core_edge_anno")
Figure2_Sample12 <- SpatialPlot(Sample_12, group.by = "core_edge_anno")

arrange()
SpatialPlot(Sample_1, group.by = "core_edge_anno")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)

ggarrange(Figure1_Sample1, Figure1_Sample2, Figure1_Sample3, Figure1_Sample4, Figure1_Sample5, Figure1_Sample6 
         , Figure1_Sample7, Figure1_Sample8, Figure1_Sample9 
          , Figure1_Sample10 , Figure1_Sample11 , Figure1_Sample12 + rremove("x.text"), 
          labels = c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10","S11","S12"),
          ncol = 2, nrow = 6)


Fig_1_1 <- ggarrange(Figure1_Sample1, Figure1_Sample2, Figure1_Sample3, Figure1_Sample4  + rremove("x.text"), 
          labels = c("S1 LNM+ ", "S2 LNM+", "S3 LNM-", "S4 LNM+"),
          ncol = 2, nrow = 2)

Fig_1_2 <- ggarrange(Figure1_Sample5, Figure1_Sample6, Figure1_Sample7, Figure1_Sample8  + rremove("x.text"), 
                     labels = c("S5 LNM- ", "S6 LNM-", "S7 LNM-", "S8 LNM-"),
                     ncol = 2, nrow = 2)

Fig_1_3 <- ggarrange(Figure1_Sample9, Figure1_Sample10, Figure1_Sample11, Figure1_Sample12  + rremove("x.text"), 
                     labels = c("S9 LNM+ ", "S10 LNM+", "S11 LNM-", "S12 LNM+"),
                     ncol = 2, nrow = 2)

Fig_2_1 <- ggarrange(Figure2_Sample1, Figure2_Sample2, Figure2_Sample3, Figure2_Sample4  + rremove("x.text"), 
                  
                     ncol = 2, nrow = 2)

Fig_2_2 <- ggarrange(Figure2_Sample5, Figure2_Sample6, Figure2_Sample7, Figure2_Sample8  + rremove("x.text"), 
                     
                     ncol = 2, nrow = 2)

Fig_2_3 <- ggarrange(Figure2_Sample9, Figure2_Sample10, Figure2_Sample11, Figure2_Sample12  + rremove("x.text"), 
                    
                     ncol = 2, nrow = 2)



DEGs_scRNA_Fibroblast <- EnhancedVolcano(DEGs_in_Fibroblast_LNMPOS, 
                                     lab = rownames(DEGs_in_Fibroblast_LNMPOS),
                                     x ="avg_log2FC", 
                                     y ="p_val_adj",
                                     title = 'DEGs in Fibroblast celltype of LNM+ before Filtering ',
                                     pCutoff = 0.01,
                                     FCcutoff = 0.58,
                                     pointSize = 1.5,
                                     labSize = 4.0)

DEGs_scRNA_Tcell <- EnhancedVolcano(DEGs_in_Tcells_LNMPOS, 
                                         lab = rownames(DEGs_in_Tcells_LNMPOS),
                                         x ="avg_log2FC", 
                                         y ="p_val_adj",
                                         title = 'DEGs in T celltype of LNM+ before Filtering ',
                                         pCutoff = 0.01,
                                         FCcutoff = 0.58,
                                         pointSize = 1.5,
                                         labSize = 4.0)


DEGs_scRNA_Cancer <- EnhancedVolcano(DEGs_inCancer_LNMPOS, 
                                         lab = rownames(DEGs_inCancer_LNMPOS),
                                         x ="avg_log2FC", 
                                         y ="p_val_adj",
                                         title = 'DEGs in Cancer celltype of LNM+ before Filtering ',
                                         pCutoff = 0.01,
                                         FCcutoff = 0.58,
                                         pointSize = 1.5,
                                         labSize = 4.0)

DEGs_scRNA_Cancer_afterFiltering <- EnhancedVolcano(Filtered_Cancer, 
                                     lab = rownames(Filtered_Cancer),
                                     x ="avg_log2FC", 
                                     y ="p_val_adj",
                                     title = 'DEGs in cancer celltype of LNM+ after Filtering ',
                                     pCutoff = 0.01,
                                     FCcutoff = 0.58,
                                     pointSize = 1.5,
                                     labSize = 4.0)
