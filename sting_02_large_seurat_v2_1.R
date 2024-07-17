setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat_v2")

library(Seurat)

combo_seurat <- readRDS("seurat_full_SCT.RDS")

# Look at MT percentages
VlnPlot(combo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 783276 cells pre filt

# Subset based on nFeature_RNA and MT%
combo_seurat <- subset(combo_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 783276 cells post filt. Should be 777188. What is happening here???

# Post filter vlnplot
VlnPlot(combo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Remove IL15 data
combo_seurat <- subset(combo_seurat, subset = Treatment == "IL15", invert = TRUE)

# Use SCTransform to normalize - got this to work with integration, this is the preferred method
combo_seurat <- SCTransform(combo_seurat)
