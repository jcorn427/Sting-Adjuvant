setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis")

library(Seurat)
library(kBET)

# Read in the data
bridging_01.data <- Read10X(data.dir = "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_3/outs/per_sample_outs/MK434_048/count/sample_filtered_feature_bc_matrix")
bridging_01 <- CreateSeuratObject(counts = bridging_01.data, project = "bridging_01")

bridging_02.data <- Read10X(data.dir = "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_G/outs/per_sample_outs/MK434_160/count/sample_filtered_feature_bc_matrix/")
bridging_02 <- CreateSeuratObject(counts = bridging_02.data, project = "bridging_02")

# merge bridging samples

bridging_combo <- merge(bridging_01, y = bridging_02, add.cell.ids = c("01", "02"), project = "bridging_combo")

#standard analysis pipeline from here
bridging_combo[["percent.mt"]] <- PercentageFeatureSet(bridging_combo, pattern = "^MT-")

VlnPlot(bridging_combo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(bridging_combo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bridging_combo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

bridging_combo <- subset(bridging_combo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Go back and recheck the previous plots

bridging_combo <- NormalizeData(bridging_combo)
bridging_combo <- FindVariableFeatures(bridging_combo)

top10 <- head(VariableFeatures(bridging_combo), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bridging_combo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

bridging_combo <- ScaleData(bridging_combo)
bridging_combo <- RunPCA(bridging_combo)

# Take a look at PCA
print(bridging_combo[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(bridging_combo, dims = 1:2, reduction = "pca")
DimPlot(bridging_combo, reduction = "pca", shuffle = TRUE)
DimPlot(bridging_combo, reduction = "pca", split.by = "orig.ident")
DimHeatmap(bridging_combo, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(bridging_combo, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(bridging_combo) # going with 20 dims

bridging_combo <- FindNeighbors(bridging_combo, dims = 1:20)
bridging_combo <- FindClusters(bridging_combo, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(bridging_combo), 5)

# UMAP time
bridging_combo <- RunUMAP(bridging_combo, dims = 1:20)
p1 <- DimPlot(bridging_combo, reduction = "umap", group.by = "orig.ident")

# Calculate kBET batch effect score
k_mtx <- GetAssayData(bridging_combo, layer = "scale.data")
batch.estimate <- kBET(t(k_mtx), bridging_combo@meta.data$orig.ident)

# calculate Principle component regression
pca.data <- prcomp(t(k_mtx), center = TRUE)
batch.silhouette <- batch_sil(pca.data, bridging_combo@meta.data$orig.ident)
batch.pca <- pcRegression(pca.data, bridging_combo@meta.data$orig.ident)
