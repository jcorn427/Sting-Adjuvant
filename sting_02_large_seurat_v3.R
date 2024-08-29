setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat_v3")

library(Seurat)
library(BPCells)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(ggplot2)
library(tm)
library(stringr)
library(clustermole)
library(dittoSeq)
library(ggrepel)
library(SeuratData)

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

# Read in target file
target <- read.csv("../target_file.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Crash recovery:
pbmc_seurat <- readRDS("pbmc_seurat_mt_filt.RDS")



# Load data from h5 file that was created using cellranger aggr
pbmc.data <- open_matrix_10x_hdf5(path = "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_full_h5/outs/count/filtered_feature_bc_matrix.h5")

# write the matrix to a directory
write_matrix_dir(mat = pbmc.data, dir = "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat_v3/sample_files/")

# Load matrix from disk, convert Ensembl to gene names
pbmc.mat <- open_matrix_dir(dir = "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat_v3/sample_files/")
pbmc.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = pbmc.mat, species = "human")

# Add some metadata columns to the target file
target$Genotype <- removeNumbers(target$Donor)
target$Sample <- row.names(target)

#Create Seurate pbmc_seuratect
pbmc_seurat <- CreateSeuratObject(counts = pbmc.mat)

# Generate metadata data frame to add to seurat pbmc_seuratect
meta_data <- pbmc_seurat@meta.data

barcodes <- strsplit(row.names(meta_data), split = "-")

n <- 1
l <- data.frame()

for (i in row.names(meta_data)){
  #print(i)
  if (str_extract(i, '\\b\\w+$') == n){
    for (j in row.names(target)){
      if (as.numeric(str_split(j, '_', simplify = TRUE)[,2]) == n){
        k <- cbind(meta_data[i,],target[j,])
        l <- rbind(l, k)
      }
    
    }
  }
  
  if (str_extract(i, '\\b\\w+$') != n){
    n <- n+1
  }
}

# Add metadata to seurat object
cell_data <- l %>% select(Treatment, Timepoint, Genotype, Sample)
pbmc_seurat <- AddMetaData(SeuratObject = pbmc_seurat, metadata = cell_data)


# Get MT percentages
pbmc_seurat <- PercentageFeatureSet(pbmc_seurat, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(pbmc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(pbmc_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Subset based on nFeature_RNA and MT%
pbmc_seurat <- subset(pbmc_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

VlnPlot(pbmc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Remove IL15 data
pbmc_seurat <- subset(pbmc_seurat, subset = Treatment == "IL15", invert = TRUE)

saveRDS(pbmc_seurat, file = "pbmc_seurat_mt_filt.RDS")

# Generate bridging pbmc_seuratect for individual analysis
bridging_seurat <- subset(x = pbmc_seurat, subset = Treatment == "BRIDGING")

## Subset, no bridging ##
pbmc_seurat <- subset(x = pbmc_seurat, subset = Treatment != "BRIDGING")

# Subset, no HAQ for initial analysis
pbmc_seurat <- subset(x = pbmc_seurat, subset = Genotype != "HAQ")

# Normalize data
pbmc_seurat <- NormalizeData(pbmc_seurat)

# Split into individual samples
pbmc_seurat[["RNA"]] <- split(pbmc_seurat[["RNA"]], f = pbmc_seurat$Sample)


# Find variable features for each sample
pbmc_seurat <- FindVariableFeatures(pbmc_seurat)

# Sketch 500 cells from each sample
pbmc_seurat <- SketchData(object = pbmc_seurat, ncells = 500, method = "LeverageScore", sketched.assay = "sketch")

saveRDS(pbmc_seurat, file = "pbmc_seurat_160_sketch.RDS")

# Make sure the active assaly is "sketch" then run standard analysis pipeline
DefaultAssay(pbmc_seurat) <- "sketch"
pbmc_seurat <- FindVariableFeatures(pbmc_seurat)
pbmc_seurat <- ScaleData(pbmc_seurat)
pbmc_seurat <- RunPCA(pbmc_seurat)

ElbowPlot(pbmc_seurat, ndims = 50)

# integrate the datasets
pbmc_seurat <- IntegrateLayers(pbmc_seurat, method = HarmonyIntegration, orig = "pca", new.reduction = "integrated.harmony")


pbmc_seurat <- FindNeighbors(pbmc_seurat, reduction = "integrated.harmony", dims = 1:50)

#Trying different resolution, 0.2 was too few clusters. Trying 0.4, trying 0.3, 0.5, 0.6, 1.0, 1.4
pbmc_seurat <- FindClusters(pbmc_seurat, resolution = 1.4)

pbmc_seurat <- RunUMAP(pbmc_seurat, dims = 1:50, return.model = T)
# DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap")

# Integrate the full dataset
pbmc_seurat <- ProjectIntegration(object = pbmc_seurat, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.harmony")





# Extend sketch results to full dataset
pbmc_seurat <- ProjectData(
  object = pbmc_seurat,
  assay = "RNA",
  full.reduction = "integrated.harmony.full",
  sketched.assay = "sketch",
  sketched.reduction = "integrated.harmony.full",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)

# Run UMAP on full dataset
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "integrated.harmony.full", dims = 1:50, reduction.name = "umap.full",
                  reduction.key = "UMAP_full_")

# Look at full dataset:
DefaultAssay(pbmc_seurat) <- "RNA"

DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap.full", group.by = "cluster_full", alpha = 0.1, shuffle = TRUE)
DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "full.umap", group.by = "cluster_full", alpha = 0.1, shuffle = TRUE)

DefaultAssay(pbmc_seurat) <- "sketch"
DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap.full", group.by = "cluster_full", alpha = 0.1, shuffle = TRUE)

Idents(pbmc_seurat) <- "cluster_full"
DimPlot(pbmc_seurat, split.by = "cluster_full", shuffle = TRUE, ncol = 4, reduction = "umap.full", alpha = 0.1)
DimPlot(pbmc_seurat, split.by = "cluster_full", shuffle = TRUE, ncol = 4, reduction = "umap.full")

# Saving clustered seurat object pre-cell types
saveRDS(pbmc_seurat, file = "pbmc_seurat_clustered.RDS")
pbmc_seurat <- readRDS("pbmc_seurat_clustered.RDS")

# Running cluster mole

agg_exp_mat <- AggregateExpression(pbmc_seurat)
agg_exp_mat_RNA <- as.matrix(agg_exp_mat$RNA)

agg_exp_mat_RNA <- log1p(agg_exp_mat_RNA)

enrich_tbl <- clustermole_enrichment(expr_mat = agg_exp_mat_RNA, species = "hs")

write.csv(enrich_tbl, file = "cluster_mole_enrichments.csv")

# cluster mole cell assignments
# reorder levels of seurat object
levels(pbmc_seurat) <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27')

new.cluster.ids <- c("CD4+", "Naive CD4+", "Macrophage", "Tfh T", "Naive B", "Naive CD4+", "CD4+", "CD8+", "NK", "NK", "Naive B", "Memory B", "Progenitor", "NK", "DC", "CD8+", "Endothelial", "CD4+", "Progenitor", "DC", "T", "Treg", "CD4+", "NK", "Plasmacytoid DC", "Naive B", "Naive B", "Megakaryocytes")
names(new.cluster.ids) <- levels(pbmc_seurat)
pbmc_seurat <- RenameIdents(pbmc_seurat, new.cluster.ids)
pbmc_seurat@meta.data$clustermole_ids <- Idents(pbmc_seurat)

DimPlot(pbmc_seurat, reduction = "umap.full", split.by = "clustermole_ids", alpha = 0.1, shuffle = T, label = T, ncol = 4)

DimPlot(pbmc_seurat, reduction = "umap.full", alpha = 0.1, label = TRUE)

# Rejoin layers
pbmc_seurat[["sketch"]] <- JoinLayers(pbmc_seurat[["sketch"]])
pbmc_seurat[["RNA"]] <- JoinLayers(pbmc_seurat[["RNA"]])

# Run azimuth
pbmc_seurat <- RunAzimuth(pbmc_seurat, reference = "pbmcref", assay = "RNA")

Idents(pbmc_seurat) <- "predicted.celltype.l1"
DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap.full", group.by = "predicted.celltype.l1", alpha = 0.1, shuffle = TRUE)
DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap.full", group.by = "predicted.celltype.l2", alpha = 0.1, shuffle = TRUE, repel = T)
DimPlot(pbmc_seurat, label = T, label.size = 3, reduction = "umap.full", group.by = "predicted.celltype.l3", alpha = 0.1, shuffle = TRUE, repel = T)


FeaturePlot(pbmc_seurat, "CD3D", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "CD14", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "NCAM1", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "FCGR3A", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "CD8A", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "CD4", reduction = "umap.full")
FeaturePlot(pbmc_seurat, "CD19", reduction = "umap.full")

Idents(pbmc_seurat) <- "Treatment"
DimPlot(pbmc_seurat, reduction = "pca", shuffle = TRUE)

#### Verify predicted cell type assignments using canonical marker genes ####
tcells <- FeaturePlot(pbmc_seurat, features = "CD3D", reduction = "umap.full")
tcells

cd8_t <- FeaturePlot(pbmc_seurat, features = c("CD8A", "CD8B"), reduction = "umap.full")
cd8_t

cd4_t <- FeaturePlot(pbmc_seurat, features = c("CD4"), reduction = "umap.full")
cd4_t

cd8_effector <- FeaturePlot(pbmc_seurat, features = c("GZMK", "GZMH", "PRF1", "CCL5"), reduction = "umap.full")
cd8_effector

cd8_memory <- FeaturePlot(pbmc_seurat, features = c("ITGB1"), reduction = "umap.full")
cd8_memory

cd8_naive <- FeaturePlot(pbmc_seurat, features = c("CCR7"), reduction = "umap.full")
cd8_naive

cd4_naive <- FeaturePlot(pbmc_seurat, features = c("IL7R", "CCR7"), reduction = "umap.full")
cd4_naive

cd4_memory <- FeaturePlot(pbmc_seurat, features = c("IL7R", "S100A4"), reduction = "umap.full")
cd4_memory

treg <- FeaturePlot(pbmc_seurat, features = c("FOXP3", "IL2RA"), reduction = "umap.full")
treg

MAIT <- FeaturePlot(pbmc_seurat, features = c("SLC4A10", "TCRAV7S2"), reduction = "umap.full")
MAIT

tfh_t <- FeaturePlot(pbmc_seurat, features = c("CXCR5", "PDCD1"), reduction = "umap.full")
tfh_t

#gamma_delta <- FeaturePlot(pbmc_seurat, features = c("TRGV9", "TRDV2"), reduction = "umap.full")
#gamma_delta

b_cells <- FeaturePlot(pbmc_seurat, features = c("CD79A", "CD79B"), reduction = "umap.full")
b_cells

b_naive <- FeaturePlot(pbmc_seurat, features = c("CD27"), reduction = "umap.full")
b_naive

b_plasma <- FeaturePlot(pbmc_seurat, features = c("SDC1", "MZB1", "XBP1"), reduction = "umap.full")
b_plasma

mono_classic <- FeaturePlot(pbmc_seurat, features = c("CD14", "LYZ"), reduction = "umap.full")
mono_classic

mono_nonclassic <- FeaturePlot(pbmc_seurat, features = c("FCGR3A", "MS4A7"), reduction = "umap.full")
mono_nonclassic

dcs <- FeaturePlot(pbmc_seurat, features = c("IL3RA", "CLEC4C", "CST3", "NRP1"), reduction = "umap.full")
dcs

dc_myeloid <- FeaturePlot(pbmc_seurat, features = c("FCER1A", "CD1C"), reduction = "umap.full")
dc_myeloid

dc_plasmacytoid <- FeaturePlot(pbmc_seurat, features = c("FCER1A", "LILRA4"), reduction = "umap.full")
dc_plasmacytoid

nk_cells <- FeaturePlot(pbmc_seurat, features = c("KLRB1", "KLRC1", "KLRD1", "GNLY", "NKG7", "NCAM1"), reduction = "umap.full")
nk_cells

platelets <- FeaturePlot(pbmc_seurat, features = c("PPBP"), reduction = "umap.full")
platelets

# Generating figures for update
dittoDimHex(pbmc_seurat, reduction.use = "umap.full", "cluster_full")
dittoDimHex(pbmc_seurat, reduction.use = "umap.full")

# Mapping to Seurat PBMC CITE-seq reference

reference <- readRDS("pbmc_multimodal_2023.rds")

anchor <- FindTransferAnchors(
  reference = reference,
  query = pbmc_seurat,
  reference.reduction = "spca",
  normalization.method = "SCT",
  dims = 1:50
)
pbmc_seurat <- MapQuery(
  anchorset = anchor,
  query = pbmc_seurat,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)

# Look at mapping results and compare to reference

plot1 <- DimPlot(pbmc_seurat, reduction = "ref.umap", group.by = "predicted.celltype.l2", alpha = 0.1, label = TRUE, repel = T)
DimPlot(pbmc_seurat, reduction = "umap.full", group.by = "predicted.celltype.l2", alpha = 0.1, label = TRUE, repel = T)
DimPlot(pbmc_seurat, reduction = "umap.full", group.by = "predicted.celltype.l1", alpha = 0.1, label = TRUE, repel = T)
DimPlot(pbmc_seurat, reduction = "umap.full", group.by = "predicted.celltype.l3", alpha = 0.1, label = TRUE, repel = T)

plot2 <- DimPlot(reference, group.by = "celltype.l2", reduction = "wnn.umap", label = T, alpha = 0.1, repel = T)

plot1 + plot2

# Make final cell calls
Idents(pbmc_seurat) <- "cluster_full"

levels(pbmc_seurat) <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27')
new.cluster.ids <- c("CD4 Naive", "CD4 Naive", "CD14 Mono", "CD4", "B Plasma", "CD4 Naive", "CD4 Memory", "CD8 Naive", "CD8 Memory", "CD8", "B Intermediate", "B Naive", "T", "NK", "cDC2", "CD8 Effector", "CD16 Mono", "CD4", "Progenitor", "cDC2", "MAIT", "Treg", "T", "NK", "cDC1", "B Naive", "B Naive", "Platelet")


names(new.cluster.ids) <- levels(pbmc_seurat)
pbmc_seurat <- RenameIdents(pbmc_seurat, new.cluster.ids)
pbmc_seurat@meta.data$final_cell_ids <- Idents(pbmc_seurat)
DimPlot(pbmc_seurat, reduction = "umap.full", label = TRUE, pt.size = 0.5)
DimPlot(pbmc_seurat, reduction = "umap.full", group.by = "final_cell_ids", alpha = 0.1, label = TRUE)

#save seurat object with final cell calls
saveRDS(pbmc_seurat, "pbmc_seurat_final_ids.RDS")

####Perform DE between groups in wt only####

dmso.mock.markers <- FindMarkers(pbmc_seurat, ident.1 = "DMSO", ident.2 = "Mock", group.by = "Treatment")
DA5nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "DiABZi_5nM", ident.2 = "DMSO", group.by = "Treatment")
DA25nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "DiABZi_25nM", ident.2 = "DMSO", group.by = "Treatment")
ST002_25nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST002_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST002_50nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST002_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_25nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST012_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_50nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST012_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_25nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST020_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_50nm.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "ST020_50uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_25uM.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "INI3069_25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50uM.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "INI3069_50uM", ident.2 = "DMSO", group.by = "Treatment")
treat_1148_10uM.dmso.markers <- FindMarkers(pbmc_seurat, ident.1 = "1148_10uM", ident.2 = "DMSO", group.by = "Treatment")

write.csv(dmso.mock.markers, "./dmso_vs_mock_markers.csv")
write.csv(DA5nm.dmso.markers, "./DA5nm_vs_dmso_markers.csv")
write.csv(DA25nm.dmso.markers, "./DA25nm_vs_dmso_markers.csv")
write.csv(ST002_25nm.dmso.markers, "ST002_25nm.dmso_markers.csv")
write.csv(ST002_50nm.dmso.markers, "ST002_50nm.dmso_markers.csv")
write.csv(ST012_25nm.dmso.markers, "ST012_25nm.dmso_markers.csv")
write.csv(ST012_50nm.dmso.markers, "ST012_50nm.dmso_markers.csv")
write.csv(ST020_25nm.dmso.markers, "ST020_25nm.dmso_markers.csv")
write.csv(ST020_50nm.dmso.markers, "ST020_50nm.dmso_markers.csv")
write.csv(INI3069_25uM.dmso.markers, "INI3069_25uM.dmso_markers.csv")
write.csv(INI3069_50uM.dmso.markers, "INI3069_50uM.dmso_markers.csv")
write.csv(treat_1148_10uM.dmso.markers, "treat_1148_10uM.dmso_markers.csv")


# INI3069 analysis
DoHeatmap(pbmc_seurat, features = VariableFeatures(pbmc_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "final_cell_ids", label = F)

pbmc_seurat <- FindVariableFeatures(pbmc_seurat)
pbmc_seurat <- ScaleData(pbmc_seurat)

dittoDimPlot(pbmc_seurat, "final_cell_ids", reduction.use = "umap.full")
dittoBarPlot(pbmc_seurat, "final_cell_ids", group.by = "Treatment")
dittoBarPlot(pbmc_seurat, "final_cell_ids", group.by = "Treatment", scale = "count")
dittoBarPlot(pbmc_seurat, "final_cell_ids", group.by = "Timepoint")
dittoHeatmap(pbmc_seurat, genes = VariableFeatures(pbmc_seurat)[1:50], cells.use = 1:500)

Idents(pbmc_seurat) <- "Treatment"

INI3069_25uM_seurat <- subset(x = pbmc_seurat, idents = "INI3069_25uM")

DoHeatmap(INI3069_25uM_seurat, features = VariableFeatures(pbmc_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "final_cell_ids") + NoLegend()

INI3069_50uM_seurat <- subset(x = pbmc_seurat, idents = "INI3069_50uM")

DiABZI_5nM_seurat <- subset(x = pbmc_seurat, idents = "DiABZi_5nM")

DiABZI_25nM_seurat <- subset(x = pbmc_seurat, idents = "DiABZi_25nM")

Idents(INI3069_25uM_seurat) <- "Timepoint"
Idents(INI3069_50uM_seurat) <- "Timepoint"
Idents(DiABZI_25nM_seurat) <- "Timepoint"
Idents(DiABZI_5nM_seurat) <- "Timepoint"


INI3069_25uM_seurat_6 <- subset(x = INI3069_25uM_seurat, idents = "6hr")
INI3069_50uM_seurat_6 <- subset(x = INI3069_50uM_seurat, idents = "6hr")
DiABZI_5nM_seurat_6 <- subset(x = DiABZI_5nM_seurat, idents = "6hr")
DiABZI_25nM_seurat_6 <- subset(x = DiABZI_25nM_seurat, idents = "6hr")

INI3069_25uM_seurat_24 <- subset(x = INI3069_25uM_seurat, idents = "24hr")
INI3069_50uM_seurat_24 <- subset(x = INI3069_50uM_seurat, idents = "24hr")
DiABZI_5nM_seurat_24 <- subset(x = DiABZI_5nM_seurat, idents = "24hr")
DiABZI_25nM_seurat_24 <- subset(x = DiABZI_25nM_seurat, idents = "24hr")

DoHeatmap(INI3069_25uM_seurat_6, features = VariableFeatures(INI3069_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "final_cell_ids", label = F)
DoHeatmap(INI3069_50uM_seurat_6, features = VariableFeatures(INI3069_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_5nM_seurat_6, features = VariableFeatures(DiABZI_5nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_25nM_seurat_6, features = VariableFeatures(DiABZI_25nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DoHeatmap(INI3069_25uM_seurat_24, features = VariableFeatures(INI3069_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(INI3069_50uM_seurat_24, features = VariableFeatures(INI3069_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_5nM_seurat_24, features = VariableFeatures(DiABZI_5nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_25nM_seurat_24, features = VariableFeatures(DiABZI_25nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DimPlot(pbmc_seurat, reduction = "umap.full", group.by = "cluster_full", alpha = 0.1, label = TRUE, repel = F)
