setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_sub/Cornelius_analysis")

library(Seurat)
library(BPCells)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(glmGamPoi)
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

# specify that you would like to create a Seurat v5 assay
# note that we require setting this option to ensure that existing pipelines are not affected
options(Seurat.object.assay.version = "v5")

# Read in target file
target <- read.csv("./target_file.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Make vector of directories where all h5 files are stored
directories <- c("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/"
)




# Loop through h5 files and output BPCells matrices on-disk

# Read in the data and initialize the Seurat objects with the raw (non-normalized data)

data.list <- c()
cell_ids <- c()
cell_meta <- c()
sample_meta <- c()

for (j in directories){
  print(j)
  for (i in list.files(j)){
    print(i)
    seurat_data <- open_matrix_10x_hdf5(path = paste0(j,i,"/count/sample_raw_feature_bc_matrix.h5"))
    #print(summary(colSums(seurat_data)))
    # Write the matrix to a directory
    write_matrix_dir(
      mat = seurat_data,
      dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_sub/Cornelius_analysis/sample_files/",i)
    )
    
    # Load matrix from directory
    seurat.mat <- open_matrix_dir(dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_sub/Cornelius_analysis/sample_files/",i))
    seurat.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = seurat.mat, species = "human")
    
    target_test <- target[row.names(target)== i,]
    sample_meta <- rbind(sample_meta, target_test)
    
    seurat.meta <- rownames(t(seurat.mat))
    
    cell_temp <- cbind(target_test[rep(1, length(seurat.meta)),], data.frame(barcode = seurat.meta))
    
    rownames(cell_temp) <- cell_temp$barcode
    cell_temp$barcode <- NULL
    
    cell_meta[[i]] <- cell_temp
    
    
    
    
    data.list[[i]] <- seurat.mat
    
  }
}



new_suffixes <- c("-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15", "-16")
# Create 1 seurat object from list of matrices

data_mod <- Replace_Suffix(data = data.list, current_suffix = "-1", new_suffix = new_suffixes)

merged.obj <- CreateSeuratObject(counts = data_mod)


# Create metadata dataframe to add, and add it to the seurat object

cell_meta <- lapply(cell_meta, function(x){t(x)})

cell_mod <- Replace_Suffix(data = cell_meta, current_suffix = "-1", new_suffix = new_suffixes)

cell_mod <- lapply(cell_mod, function(x){t(x)})

cell_all <- do.call(rbind, cell_mod)

merged.obj <- AddMetaData(object = merged.obj, metadata = cell_all)

# Perform QC on merged object
merged.obj <- PercentageFeatureSet(merged.obj, "^MT-", col.name = "percent_mito")

merged.obj_filt <- subset(merged.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 5)

# saveRDS(merged.obj_filt, file = "merged.obj_filt.RDS")

merged.obj_filt <- readRDS("merged.obj_filt.RDS")


# Sketch based pipeline
merged.obj_filt <- NormalizeData(merged.obj_filt)
#saveRDS(merged.obj_filt, file = "merged.obj_filt_norm.RDS")
merged.obj_filt <- readRDS("merged.obj_filt_norm.RDS")

merged.obj_filt <- FindVariableFeatures(merged.obj_filt)
# saveRDS(merged.obj_filt, file = "merged.obj_filt_norm_var.RDS")
merged.obj_filt <- readRDS("merged.obj_filt_norm_var.RDS")

merged.obj_filt <- SketchData(object = merged.obj_filt, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch", verbose = TRUE)
# saveRDS(merged.obj_filt, file = "merged.obj_filt_norm_var_sketch.rds")
merged.obj_filt


DefaultAssay(merged.obj_filt) <- "sketch"
merged.obj_filt <- FindVariableFeatures(merged.obj_filt)
merged.obj_filt <- ScaleData(merged.obj_filt)
merged.obj_filt <- RunPCA(merged.obj_filt)
merged.obj_filt <- FindNeighbors(merged.obj_filt, dims = 1:50)
merged.obj_filt <- FindClusters(merged.obj_filt, resolution = 2)

merged.obj_filt <- RunUMAP(merged.obj_filt, dims = 1:50, return.model = T)
DimPlot(merged.obj_filt, label = T, label.size = 3, reduction = "umap") + NoLegend()

FeaturePlot(
  object = merged.obj_filt,
  features = c(
    "CD4"
  )
)

merged.obj_filt <- ProjectData(
  object = merged.obj_filt,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(merged.obj_filt) <- "RNA"

DimPlot(merged.obj_filt, label = T, label.size = 3, reduction = "full.umap", group.by = "cluster_full", alpha = 0.1) + NoLegend()
