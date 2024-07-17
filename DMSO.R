setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis")

library(Seurat)
library(scCustomize)
library(ggplot2)
library(sctransform)
library(harmony)
library(SeuratWrappers)
library(batchelor)
library(reticulate)
library(stringr)
library(Azimuth)
library(SeuratDisk)
library(clustermole)
library(HGNChelper)
library(dplyr)

# Source scType functions and prepare gene sets
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


options(future.globals.maxSize = 3e+09)

sample_directories <- c("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/MK434_001/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/MK434_012/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_2/outs/per_sample_outs/MK434_024/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_3/outs/per_sample_outs/MK434_035/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_A/outs/per_sample_outs/MK434_049/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_A/outs/per_sample_outs/MK434_060/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_B/outs/per_sample_outs/MK434_071/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_C/outs/per_sample_outs/MK434_082/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_C/outs/per_sample_outs/MK434_093/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_D/outs/per_sample_outs/MK434_104/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_E/outs/per_sample_outs/MK434_115/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_E/outs/per_sample_outs/MK434_126/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_F/outs/per_sample_outs/MK434_137/count/sample_filtered_feature_bc_matrix",
                        "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_G/outs/per_sample_outs/MK434_148/count/sample_filtered_feature_bc_matrix")

sample_names <- c("haq2_DMSO_6hr",
                  "wt5_DMSO_6hr",
                  "haq2_DMSO_24hr",
                  "wt5_DMSO_24hr",
                  "wt1_DMSO_6hr",
                  "wt2_DMSO_6hr",
                  "haq1_DMSO_6hr",
                  "wt3_DMSO_6hr",
                  "wt4_DMSO_6hr",
                  "wt1_DMSO_24hr",
                  "wt2_DMSO_24hr",
                  "haq1_DMSO_24hr",
                  "wt3_DMSO_24hr",
                  "wt4_DMSO_24hr")

# Read in data, make seurat objects, rename cells to eliminate non-unique barcodes, add them to a list

x <- 0
data.list <- c()
for (i in sample_directories){
  print(i)
  temp.data <- Read10X(data.dir = i)
  x <- x+1
  print(paste(x, sample_names[x]))
  temp.seurat <- CreateSeuratObject(counts = temp.data, project = sample_names[x])
  temp.seurat <- Replace_Suffix(temp.seurat, current_suffix = "-1", new_suffix = paste0("-", sample_names[x]))
  data.list[[sample_names[x]]] <- temp.seurat
}

# Merge seurat objects
combo_seurat <- Merge_Seurat_List(data.list)


# Rename layers
combo_seurat[['RNA']] <- JoinLayers(object = combo_seurat[['RNA']])
combo_seurat[['RNA']] <- split(combo_seurat[['RNA']], f = combo_seurat$orig.ident)

# Get MT percentages
combo_seurat <- PercentageFeatureSet(combo_seurat, pattern = "^MT-", col.name = "percent.mt")

# VlnPlot(combo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# Subset based on nFeature_RNA and MT%
combo_seurat <- subset(combo_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Run standard processing steps - only run this if having issues with SCT
combo_seurat <- NormalizeData(combo_seurat)
combo_seurat <- FindVariableFeatures(combo_seurat)
combo_seurat <- ScaleData(combo_seurat)

# Use SCTransform to normalize - got this to work with integration, this is the preferred method
combo_seurat <- SCTransform(combo_seurat)




# Run PCA
combo_seurat <- RunPCA(combo_seurat)

# Take a look at PCA
# print(combo_seurat[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(combo_seurat, dims = 1:2, reduction = "pca")
# DimPlot(combo_seurat, reduction = "pca")
# DimHeatmap(combo_seurat, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(combo_seurat, dims = 1:15, cells = 500, balanced = TRUE)

# ElbowPlot(combo_seurat) # going with 20 dims

# Run UMAP
combo_seurat <- RunUMAP(combo_seurat, dims = 1:20)

# Find neigbors and clusters
combo_seurat <- FindNeighbors(combo_seurat, dims = 1:20)
combo_seurat <- FindClusters(combo_seurat)
# DimPlot(combo_seurat, label = TRUE)
# DimPlot(combo_seurat, group.by = "orig.ident")
# DimPlot(combo_seurat, group.by = "orig.ident", split.by = "orig.ident", ncol = 3)

# There's clear separation between the different samples. Moving on to integration to get cleaner results
# Going to look at several different integration methods to determine which would be best

# CCA
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca", k.weight = 91
)

combo_seurat <- FindNeighbors(combo_seurat, reduction = "integrated.cca", dims = 1:20)
combo_seurat <- FindClusters(combo_seurat, resolution = 2, cluster.name = "cca_clusters")

combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

p1 <- DimPlot(
  combo_seurat,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  combine = FALSE, label.size = 2
)

# RPCA
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT", assay = "SCT", k.weight = 84
)

combo_seurat <- FindNeighbors(combo_seurat, reduction = "integrated.rpca", dims = 1:20)
combo_seurat <- FindClusters(combo_seurat, resolution = 2, cluster.name = "rpca_clusters")

combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")

DimPlot(
  combo_seurat,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters"),
  combine = FALSE, label.size = 2
)

# Harmony
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony", verbose = TRUE
)

combo_seurat <- FindNeighbors(combo_seurat, reduction = "harmony", dims = 1:20)
combo_seurat <- FindClusters(combo_seurat, resolution = 2, cluster.name = "harmony_clusters")

combo_seurat <- RunUMAP(combo_seurat, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

p3 <- DimPlot(
  combo_seurat,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"),
  combine = FALSE, label.size = 2
)

# FastMNN
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn", verbose = TRUE
)

combo_seurat <- FindNeighbors(combo_seurat, reduction = "integrated.mnn", dims = 1:20)
combo_seurat <- FindClusters(combo_seurat, resolution = 2, cluster.name = "mnn_clusters")

combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.mnn", dims = 1:20, reduction.name = "umap.mnn")

p4 <- DimPlot(
  combo_seurat,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)

# scVI
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = scVIIntegration,
  new.reduction = "integrated.scvi", 
  conda_env = "/home/jcorn2/miniforge3/envs/scvi-env2",
  verbose = TRUE
)

combo_seurat <- FindNeighbors(combo_seurat, reduction = "integrated.scvi", dims = 1:20)
combo_seurat <- FindClusters(combo_seurat, resolution = 2, cluster.name = "scvi_clusters")

combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.scvi", dims = 1:20, reduction.name = "umap.scvi")

p5 <- DimPlot(
  combo_seurat,
  reduction = "umap.scvi",
  group.by = c("scvi_clusters"),
  combine = FALSE, label.size = 2
)


### Add metadata columns ###
meta_data <- combo_seurat@meta.data$orig.ident
meta_data <- strsplit(meta_data, split = "_")
meta_data <- as.data.frame(do.call(rbind, meta_data))
meta_data$V1 <- str_sub(meta_data$V1, end = -2)


combo_seurat <- AddMetaData(combo_seurat, metadata = meta_data[,1], col.name = "genotype")
combo_seurat <- AddMetaData(combo_seurat, metadata = meta_data[,3], col.name = "time_point")

#saveRDS(combo_seurat, file = "DMSO_all.RDS")

### Look at haq vs wt and 6hr vs 24 hr PCAs and UMAPs ###
DimPlot(combo_seurat, reduction = "pca", group.by = "genotype")
DimPlot(combo_seurat, reduction = "pca", group.by = "time_point")
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "genotype")
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "time_point")

### Plan is to evaluate all 5 integration methods and compare them. However, for our Jan 12th meeting, just going with RPCA for now.
### Doing cell type calls: Up first is Azimuth using reference guided annotation

combo_seurat <- RunAzimuth(combo_seurat, reference = "pbmcref", assay = "SCT") # This errored out, waiting to hear back on my github issue

# saveRDS(combo_seurat, file = "DMSO_all_azimuth.RDS")

# Look at predicted cell types

DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l1", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l2", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "orig.ident")

### ClusterMole is next

### Get average expression matrix

#avg_exp_mat <- AverageExpression(combo_seurat)
agg_exp_mat <- AggregateExpression(combo_seurat)
agg_exp_mat_RNA <- as.matrix(agg_exp_mat$RNA)
agg_exp_mat_SCT <- as.matrix(agg_exp_mat$SCT)

agg_exp_mat_RNA <- log1p(agg_exp_mat_RNA)
agg_exp_mat_SCT <- log1p(agg_exp_mat_SCT)

enrich_tbl <- clustermole_enrichment(expr_mat = agg_exp_mat_SCT, species = "hs")

write.csv(enrich_tbl, file = "cluster_mole_enrichments3.csv")

### scType up next

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = combo_seurat[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(combo_seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(combo_seurat@meta.data[combo_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(combo_seurat@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Generate umap
combo_seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  combo_seurat@meta.data$customclassif[combo_seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(combo_seurat, reduction = "umap.rpca", label = FALSE, repel = TRUE, group.by = 'customclassif', shuffle = TRUE) 







### DE Testing
genotype_markers <- FindMarkers(combo_seurat, ident.1 = "haq", ident.2 = "wt")
