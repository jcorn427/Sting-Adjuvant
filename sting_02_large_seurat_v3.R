setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat_v3")

library(Seurat)
library(BPCells)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(ggplot2)
library(tm)
library(stringr)

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

#Create Seurate Object
pbmc_seurat <- CreateSeuratObject(counts = pbmc.mat)

# Generate metadata data frame to add to seurat object
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
pbmc_seurat <- AddMetaData(object = pbmc_seurat, metadata = cell_data)


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

# Normalize data
pbmc_seurat <- NormalizeData(pbmc_seurat)

# Split into individual layers for each sample
pbmc_seurat[["RNA"]] <- split(pbmc_seurat[["RNA"]], f = pbmc_seurat$Sample)

# Find variable features for each sample
pbmc_seurat <- FindVariableFeatures(pbmc_seurat)

# Sketch your subsample for each individual sample
pbmc_seurat <- SketchData(object = pbmc_seurat, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")




# Run PCA
combo_seurat <- RunPCA(combo_seurat)



