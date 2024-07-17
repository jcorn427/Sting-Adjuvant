setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis")

library(Seurat)
library(BPCells)
library(scCustomize)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(glmGamPoi)
# needs to be set for large dataset analysis
options(future.globals.maxSize = 3e+09)

# specify that you would like to create a Seurat v5 assay
# note that we require setting this option to ensure that existing pipelines are not affected
options(Seurat.object.assay.version = "v5")

# Read in target file
target <- read.csv("./target_file.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Make vector of directories where all h5 files are stored
directories <- c("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/", 
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_2/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_3/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_A/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_B/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_C/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_D/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_E/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_F/outs/per_sample_outs/",
                 "/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_G/outs/per_sample_outs/"
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
    seurat_data <- open_matrix_10x_hdf5(path = paste0(j,i,"/count/sample_filtered_feature_bc_matrix.h5"))
    #print(summary(colSums(seurat_data)))
    # Write the matrix to a directory
    write_matrix_dir(
      mat = seurat_data,
      dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sample_files/",i)
    )
    
    # Load matrix from directory
    seurat.mat <- open_matrix_dir(dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sample_files/",i))
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



new_suffixes <- c("-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15", "-16", 
                  "-17", "-18", "-19", "-20", "-21", "-22", "-23", "-24", "-25", "-26", "-27", "-28", "-29", "-30", "-31", "-32", 
                  "-33", "-34", "-35", "-36", "-37", "-38", "-39", "-40", "-41", "-42", "-43", "-44", "-45", "-46", "-47", "-48", 
                  "-49", "-50", "-51", "-52", "-53", "-54", "-55", "-56", "-57", "-58", "-59", "-60", "-61", "-62", "-63", "-64", 
                  "-65", "-66", "-67", "-68", "-69", "-70", "-71", "-72", "-73", "-74", "-75", "-76", "-77", "-78", "-79", "-80",
                  "-81", "-82", "-83", "-84", "-85", "-86", "-87", "-88", "-89", "-90", "-91", "-92", "-93", "-94", "-95", "-96", 
                  "-97", "-98", "-99", "-100", "-101", "-102", "-103", "-104", "-105", "-106", "-107", "-108", "-109", "-110", "-111", "-112", 
                  "-113", "-114", "-115", "-116", "-117", "-118", "-119", "-120", "-121", "-122", "-123", "-124", "-125", "-126", "-127", "-128", 
                  "-129", "-130", "-131", "-132", "-133", "-134", "-135", "-136", "-137", "-138", "-139", "-140", "-141", "-142", "-143", "-144", 
                  "-145", "-146", "-147", "-148", "-149", "-150", "-151", "-152", "-153", "-154", "-155", "-156", "-157", "-158", "-159", "-160"
)
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

#saveRDS(merged.obj_filt, file = "merged.obj_filt.RDS")

merged.obj_filt <- readRDS("merged.obj_filt.RDS")



#Non-integrated analysis

merged.obj_filt <- SCTransform(merged.obj_filt, vars.to.regress = "percent_mito")
merged.obj_filt <- RunPCA(merged.obj_filt)

ElbowPlot(merged.obj_filt, ndims = 30, reduction = "pca") # based on the elbow plot, going with 20 dims

merged.obj_filt <- RunUMAP(merged.obj_filt, dims = 1:20)



# Tried normal pipeline, R kept crashing at the PCA step
merged.obj_filt <- NormalizeData(merged.obj_filt)
#saveRDS(merged.obj_filt, file = "merged.obj_filt_norm.RDS")
merged.obj_filt <- readRDS("merged.obj_filt_norm.RDS")

merged.obj_filt <- FindVariableFeatures(merged.obj_filt)
#saveRDS(merged.obj_filt, file = "merged.obj_filt_norm_var.RDS")
merged.obj_filt <- readRDS("merged.obj_filt_norm_var.RDS")

### Continue for sketch protocol

merged.obj_filt <- SketchData(object = merged.obj_filt, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch", verbose = TRUE)
merged.obj_filt

saveRDS(merged.obj_filt, file = "merged.obj_filt_norm_var_sketch.RDS")


merged.obj_filt_scale <- ScaleData(merged.obj_filt, features = rownames(merged.obj_filt))
merged.obj_filt_scale <- RunPCA(merged.obj_filt_scale, features = VariableFeatures(merged.obj_filt_scale))

merged.obj_filt <- FindNeighbors(merged.obj_filt, dims = 1:20, reduction = "pca")
merged.obj_filt <- FindClusters(merged.obj_filt, resolution = 2, cluster.name = "unintegrated_clusters")

### End of sketch protocol




