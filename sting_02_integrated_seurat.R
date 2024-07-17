# Analysis of single cell data which was prepared using the 10x flex protocol instead of standard 3' sequencing. This is the second batch, part 1
# Script by John Cornelius

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")
library(Seurat)
library(clustermole)
library(dplyr)
library(plotly)
library(ggplot2)
library(SingleR)
library(celldex)
library(UCell)
library(seurathelpeR)
library(dittoSeq)
library(tidyverse)
library(WebGestaltR)
library(ComplexHeatmap)
source("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/heatmap3LW_function.r")

####### Look into parallelizing this all later ######
###
# library(foreach)
# library(doParallel)
# 
# #Create cluster for parallel operation
# n.cores <- parallel::detectCores() - 1
# #create the cluster
# my.cluster <- parallel::makeCluster(
#   n.cores, 
#   type = "PSOCK"
# )
# 
# #check cluster definition (optional)
# print(my.cluster)
# 
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# 
# #check if it is registered (optional)
# foreach::getDoParRegistered()
# 
# #how many workers are available? (optional)
# foreach::getDoParWorkers()
# 
# # Read in the data and initialize the Seurat objects with the raw (non-normalized data) in parallel
# seurat_list <- foreach (i = list.files("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs"), .packages = "Seurat") %dopar% {
#   #print(i)
#   seurat_data <- Read10X(data.dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/",i,"/count/sample_raw_feature_bc_matrix"))
#   seurat_obj <- CreateSeuratObject(counts = seurat_data, 
#                                    min.features = 100, 
#                                    project = i)
#   assign(i, seurat_obj)
# }
###

#source stuff for sc-type
lapply(c("HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Function to count cells per gene
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}




# Read in the data and initialize the Seurat objects with the raw (non-normalized data)
for (i in list.files("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs")){
  print(i)
  seurat_data <- Read10X(data.dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_1/outs/per_sample_outs/",i,"/count/sample_raw_feature_bc_matrix"))
  print(summary(colSums(seurat_data)))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = i)
  assign(i, seurat_obj)
}

for (i in list.files("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_2/outs/per_sample_outs")){
  #print(i)
  seurat_data <- Read10X(data.dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_2/outs/per_sample_outs/",i,"/count/sample_raw_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = i)
  assign(i, seurat_obj)
}

for (i in list.files("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_3/outs/per_sample_outs")){
  #print(i)
  seurat_data <- Read10X(data.dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/sting_02_lib_3/outs/per_sample_outs/",i,"/count/sample_raw_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = i)
  assign(i, seurat_obj)
}

MK434_001 <- AddMetaData(MK434_001, metadata = "DMSO", col.name = "Treatment")
MK434_001 <- AddMetaData(MK434_001, metadata = "6hr", col.name = "Timepoint")
MK434_001 <- AddMetaData(MK434_001, metadata = "HAQ2", col.name = "Donor")

MK434_002 <- AddMetaData(MK434_002, metadata = "DiABZi 5nM", col.name = "Treatment")
MK434_002 <- AddMetaData(MK434_002, metadata = "6hr", col.name = "Timepoint")
MK434_002 <- AddMetaData(MK434_002, metadata = "HAQ2", col.name = "Donor")

MK434_003 <- AddMetaData(MK434_003, metadata = "DiABZi 25nM", col.name = "Treatment")
MK434_003 <- AddMetaData(MK434_003, metadata = "6hr", col.name = "Timepoint")
MK434_003 <- AddMetaData(MK434_003, metadata = "HAQ2", col.name = "Donor")

MK434_004 <- AddMetaData(MK434_004, metadata = "ST002 25uM", col.name = "Treatment")
MK434_004 <- AddMetaData(MK434_004, metadata = "6hr", col.name = "Timepoint")
MK434_004 <- AddMetaData(MK434_004, metadata = "HAQ2", col.name = "Donor")

MK434_005 <- AddMetaData(MK434_005, metadata = "ST002 50uM", col.name = "Treatment")
MK434_005 <- AddMetaData(MK434_005, metadata = "6hr", col.name = "Timepoint")
MK434_005 <- AddMetaData(MK434_005, metadata = "HAQ2", col.name = "Donor")

MK434_006 <- AddMetaData(MK434_006, metadata = "ST012 25uM", col.name = "Treatment")
MK434_006 <- AddMetaData(MK434_006, metadata = "6hr", col.name = "Timepoint")
MK434_006 <- AddMetaData(MK434_006, metadata = "HAQ2", col.name = "Donor")

MK434_007 <- AddMetaData(MK434_007, metadata = "ST012 50uM", col.name = "Treatment")
MK434_007 <- AddMetaData(MK434_007, metadata = "6hr", col.name = "Timepoint")
MK434_007 <- AddMetaData(MK434_007, metadata = "HAQ2", col.name = "Donor")

MK434_008 <- AddMetaData(MK434_008, metadata = "ST020 25uM", col.name = "Treatment")
MK434_008 <- AddMetaData(MK434_008, metadata = "6hr", col.name = "Timepoint")
MK434_008 <- AddMetaData(MK434_008, metadata = "HAQ2", col.name = "Donor")

MK434_009 <- AddMetaData(MK434_009, metadata = "ST020 50uM", col.name = "Treatment")
MK434_009 <- AddMetaData(MK434_009, metadata = "6hr", col.name = "Timepoint")
MK434_009 <- AddMetaData(MK434_009, metadata = "HAQ2", col.name = "Donor")

MK434_010 <- AddMetaData(MK434_010, metadata = "INI3069 25uM", col.name = "Treatment")
MK434_010 <- AddMetaData(MK434_010, metadata = "6hr", col.name = "Timepoint")
MK434_010 <- AddMetaData(MK434_010, metadata = "HAQ2", col.name = "Donor")

MK434_011 <- AddMetaData(MK434_011, metadata = "INI3069 50uM", col.name = "Treatment")
MK434_011 <- AddMetaData(MK434_011, metadata = "6hr", col.name = "Timepoint")
MK434_011 <- AddMetaData(MK434_011, metadata = "HAQ2", col.name = "Donor")

MK434_012 <- AddMetaData(MK434_012, metadata = "DMSO", col.name = "Treatment")
MK434_012 <- AddMetaData(MK434_012, metadata = "6hr", col.name = "Timepoint")
MK434_012 <- AddMetaData(MK434_012, metadata = "WT5", col.name = "Donor")

MK434_013 <- AddMetaData(MK434_013, metadata = "DiABZi 5nM", col.name = "Treatment")
MK434_013 <- AddMetaData(MK434_013, metadata = "6hr", col.name = "Timepoint")
MK434_013 <- AddMetaData(MK434_013, metadata = "WT5", col.name = "Donor")

MK434_014 <- AddMetaData(MK434_014, metadata = "DiABZi 25nM", col.name = "Treatment")
MK434_014 <- AddMetaData(MK434_014, metadata = "6hr", col.name = "Timepoint")
MK434_014 <- AddMetaData(MK434_014, metadata = "WT5", col.name = "Donor")

MK434_015 <- AddMetaData(MK434_015, metadata = "ST002 25uM", col.name = "Treatment")
MK434_015 <- AddMetaData(MK434_015, metadata = "6hr", col.name = "Timepoint")
MK434_015 <- AddMetaData(MK434_015, metadata = "WT5", col.name = "Donor")

MK434_016 <- AddMetaData(MK434_016, metadata = "ST002 50uM", col.name = "Treatment")
MK434_016 <- AddMetaData(MK434_016, metadata = "6hr", col.name = "Timepoint")
MK434_016 <- AddMetaData(MK434_016, metadata = "WT5", col.name = "Donor")

MK434_017 <- AddMetaData(MK434_017, metadata = "ST012 25uM", col.name = "Treatment")
MK434_017 <- AddMetaData(MK434_017, metadata = "6hr", col.name = "Timepoint")
MK434_017 <- AddMetaData(MK434_017, metadata = "WT5", col.name = "Donor")

MK434_018 <- AddMetaData(MK434_018, metadata = "ST012 50uM", col.name = "Treatment")
MK434_018 <- AddMetaData(MK434_018, metadata = "6hr", col.name = "Timepoint")
MK434_018 <- AddMetaData(MK434_018, metadata = "WT5", col.name = "Donor")

MK434_019 <- AddMetaData(MK434_019, metadata = "ST020 25uM", col.name = "Treatment")
MK434_019 <- AddMetaData(MK434_019, metadata = "6hr", col.name = "Timepoint")
MK434_019 <- AddMetaData(MK434_019, metadata = "WT5", col.name = "Donor")

MK434_020 <- AddMetaData(MK434_020, metadata = "ST020 50uM", col.name = "Treatment")
MK434_020 <- AddMetaData(MK434_020, metadata = "6hr", col.name = "Timepoint")
MK434_020 <- AddMetaData(MK434_020, metadata = "WT5", col.name = "Donor")

MK434_021 <- AddMetaData(MK434_021, metadata = "INI3069 25uM", col.name = "Treatment")
MK434_021 <- AddMetaData(MK434_021, metadata = "6hr", col.name = "Timepoint")
MK434_021 <- AddMetaData(MK434_021, metadata = "WT5", col.name = "Donor")

MK434_022 <- AddMetaData(MK434_022, metadata = "INI3069 50uM", col.name = "Treatment")
MK434_022 <- AddMetaData(MK434_022, metadata = "6hr", col.name = "Timepoint")
MK434_022 <- AddMetaData(MK434_022, metadata = "WT5", col.name = "Donor")

MK434_023 <- AddMetaData(MK434_023, metadata = "1148 10uM", col.name = "Treatment")
MK434_023 <- AddMetaData(MK434_023, metadata = "6hr", col.name = "Timepoint")
MK434_023 <- AddMetaData(MK434_023, metadata = "WT5", col.name = "Donor")

MK434_024 <- AddMetaData(MK434_024, metadata = "DMSO", col.name = "Treatment")
MK434_024 <- AddMetaData(MK434_024, metadata = "24hr", col.name = "Timepoint")
MK434_024 <- AddMetaData(MK434_024, metadata = "HAQ2", col.name = "Donor")

MK434_025 <- AddMetaData(MK434_025, metadata = "DiABZi 5nM", col.name = "Treatment")
MK434_025 <- AddMetaData(MK434_025, metadata = "24hr", col.name = "Timepoint")
MK434_025 <- AddMetaData(MK434_025, metadata = "HAQ2", col.name = "Donor")

MK434_026 <- AddMetaData(MK434_026, metadata = "DiABZi 25nM", col.name = "Treatment")
MK434_026 <- AddMetaData(MK434_026, metadata = "24hr", col.name = "Timepoint")
MK434_026 <- AddMetaData(MK434_026, metadata = "HAQ2", col.name = "Donor")

MK434_027 <- AddMetaData(MK434_027, metadata = "ST002 25uM", col.name = "Treatment")
MK434_027 <- AddMetaData(MK434_027, metadata = "24hr", col.name = "Timepoint")
MK434_027 <- AddMetaData(MK434_027, metadata = "HAQ2", col.name = "Donor")

MK434_028 <- AddMetaData(MK434_028, metadata = "ST002 50uM", col.name = "Treatment")
MK434_028 <- AddMetaData(MK434_028, metadata = "24hr", col.name = "Timepoint")
MK434_028 <- AddMetaData(MK434_028, metadata = "HAQ2", col.name = "Donor")

MK434_029 <- AddMetaData(MK434_029, metadata = "ST012 25uM", col.name = "Treatment")
MK434_029 <- AddMetaData(MK434_029, metadata = "24hr", col.name = "Timepoint")
MK434_029 <- AddMetaData(MK434_029, metadata = "HAQ2", col.name = "Donor")

MK434_030 <- AddMetaData(MK434_030, metadata = "ST012 50uM", col.name = "Treatment")
MK434_030 <- AddMetaData(MK434_030, metadata = "24hr", col.name = "Timepoint")
MK434_030 <- AddMetaData(MK434_030, metadata = "HAQ2", col.name = "Donor")

MK434_031 <- AddMetaData(MK434_031, metadata = "ST020 25uM", col.name = "Treatment")
MK434_031 <- AddMetaData(MK434_031, metadata = "24hr", col.name = "Timepoint")
MK434_031 <- AddMetaData(MK434_031, metadata = "HAQ2", col.name = "Donor")

MK434_032 <- AddMetaData(MK434_032, metadata = "ST020 50uM", col.name = "Treatment")
MK434_032 <- AddMetaData(MK434_032, metadata = "24hr", col.name = "Timepoint")
MK434_032 <- AddMetaData(MK434_032, metadata = "HAQ2", col.name = "Donor")

MK434_033 <- AddMetaData(MK434_033, metadata = "INI3069 25uM", col.name = "Treatment")
MK434_033 <- AddMetaData(MK434_033, metadata = "24hr", col.name = "Timepoint")
MK434_033 <- AddMetaData(MK434_033, metadata = "HAQ2", col.name = "Donor")

MK434_034 <- AddMetaData(MK434_034, metadata = "INI3069 50uM", col.name = "Treatment")
MK434_034 <- AddMetaData(MK434_034, metadata = "24hr", col.name = "Timepoint")
MK434_034 <- AddMetaData(MK434_034, metadata = "HAQ2", col.name = "Donor")

MK434_035 <- AddMetaData(MK434_035, metadata = "DMSO", col.name = "Treatment")
MK434_035 <- AddMetaData(MK434_035, metadata = "24hr", col.name = "Timepoint")
MK434_035 <- AddMetaData(MK434_035, metadata = "WT5", col.name = "Donor")

MK434_036 <- AddMetaData(MK434_036, metadata = "DiABZi 5nM", col.name = "Treatment")
MK434_036 <- AddMetaData(MK434_036, metadata = "24hr", col.name = "Timepoint")
MK434_036 <- AddMetaData(MK434_036, metadata = "WT5", col.name = "Donor")

MK434_037 <- AddMetaData(MK434_037, metadata = "DiABZi 25nM", col.name = "Treatment")
MK434_037 <- AddMetaData(MK434_037, metadata = "24hr", col.name = "Timepoint")
MK434_037 <- AddMetaData(MK434_037, metadata = "WT5", col.name = "Donor")

MK434_038 <- AddMetaData(MK434_038, metadata = "ST002 25uM", col.name = "Treatment")
MK434_038 <- AddMetaData(MK434_038, metadata = "24hr", col.name = "Timepoint")
MK434_038 <- AddMetaData(MK434_038, metadata = "WT5", col.name = "Donor")

MK434_039 <- AddMetaData(MK434_039, metadata = "ST002 50uM", col.name = "Treatment")
MK434_039 <- AddMetaData(MK434_039, metadata = "24hr", col.name = "Timepoint")
MK434_039 <- AddMetaData(MK434_039, metadata = "WT5", col.name = "Donor")

MK434_040 <- AddMetaData(MK434_040, metadata = "ST012 25uM", col.name = "Treatment")
MK434_040 <- AddMetaData(MK434_040, metadata = "24hr", col.name = "Timepoint")
MK434_040 <- AddMetaData(MK434_040, metadata = "WT5", col.name = "Donor")

MK434_041 <- AddMetaData(MK434_041, metadata = "ST012 50uM", col.name = "Treatment")
MK434_041 <- AddMetaData(MK434_041, metadata = "24hr", col.name = "Timepoint")
MK434_041 <- AddMetaData(MK434_041, metadata = "WT5", col.name = "Donor")

MK434_042 <- AddMetaData(MK434_042, metadata = "ST020 25uM", col.name = "Treatment")
MK434_042 <- AddMetaData(MK434_042, metadata = "24hr", col.name = "Timepoint")
MK434_042 <- AddMetaData(MK434_042, metadata = "WT5", col.name = "Donor")

MK434_043 <- AddMetaData(MK434_043, metadata = "ST020 50uM", col.name = "Treatment")
MK434_043 <- AddMetaData(MK434_043, metadata = "24hr", col.name = "Timepoint")
MK434_043 <- AddMetaData(MK434_043, metadata = "WT5", col.name = "Donor")

MK434_044 <- AddMetaData(MK434_044, metadata = "INI3069 25uM", col.name = "Treatment")
MK434_044 <- AddMetaData(MK434_044, metadata = "24hr", col.name = "Timepoint")
MK434_044 <- AddMetaData(MK434_044, metadata = "WT5", col.name = "Donor")

MK434_045 <- AddMetaData(MK434_045, metadata = "INI3069 50uM", col.name = "Treatment")
MK434_045 <- AddMetaData(MK434_045, metadata = "24hr", col.name = "Timepoint")
MK434_045 <- AddMetaData(MK434_045, metadata = "WT5", col.name = "Donor")

MK434_046 <- AddMetaData(MK434_046, metadata = "1148 10uM", col.name = "Treatment")
MK434_046 <- AddMetaData(MK434_046, metadata = "24hr", col.name = "Timepoint")
MK434_046 <- AddMetaData(MK434_046, metadata = "WT5", col.name = "Donor")

MK434_047 <- AddMetaData(MK434_047, metadata = "Mock", col.name = "Treatment")
MK434_047 <- AddMetaData(MK434_047, metadata = "24hr", col.name = "Timepoint")
MK434_047 <- AddMetaData(MK434_047, metadata = "WT5", col.name = "Donor")

MK434_048 <- AddMetaData(MK434_048, metadata = "BRIDGING", col.name = "Treatment")
MK434_048 <- AddMetaData(MK434_048, metadata = "BRIDGING", col.name = "Timepoint")
MK434_048 <- AddMetaData(MK434_048, metadata = "BRIDGING", col.name = "Donor")

# Create list of seurat objects
data.list <- list(MK434_001, MK434_002, MK434_003, MK434_004, MK434_005, MK434_006, MK434_007, MK434_008, MK434_009, MK434_010, MK434_011, MK434_012, MK434_013, MK434_014, MK434_015, MK434_016,
                  MK434_017, MK434_018, MK434_019, MK434_020, MK434_021, MK434_022, MK434_023, MK434_024, MK434_025, MK434_026, MK434_027, MK434_028, MK434_029, MK434_030, MK434_031, MK434_032,
                  MK434_033, MK434_034, MK434_035, MK434_036, MK434_037, MK434_038, MK434_039, MK434_040, MK434_041, MK434_042, MK434_043, MK434_044, MK434_045, MK434_046, MK434_047, MK434_048)

#Calculate % mt
data.list2 <- lapply(X = data.list, FUN = function(x){
  x$percent_MT <- PercentageFeatureSet(x, pattern = "^MT-")
  return(x)
}
  )

# VlnPlot(data.list2[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[11]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[12]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[13]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[14]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[15]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[16]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[17]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[18]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[19]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[20]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[21]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[22]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[23]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[24]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[25]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[26]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[27]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[28]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[29]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[30]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[31]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[32]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[33]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[34]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[35]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[36]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[37]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[38]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[39]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[40]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[41]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[42]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[43]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[44]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[45]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[46]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[47]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)
# VlnPlot(data.list2[[48]], features = c("nFeature_RNA", "nCount_RNA", "percent_MT"), ncol = 3)

#Filter based on nFeatures and MT%
data.list2[[1]] <- subset(data.list2[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[2]] <- subset(data.list2[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[3]] <- subset(data.list2[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[4]] <- subset(data.list2[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[5]] <- subset(data.list2[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[6]] <- subset(data.list2[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[7]] <- subset(data.list2[[7]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[8]] <- subset(data.list2[[8]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[9]] <- subset(data.list2[[9]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[10]] <- subset(data.list2[[10]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[11]] <- subset(data.list2[[11]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[12]] <- subset(data.list2[[12]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[13]] <- subset(data.list2[[13]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[14]] <- subset(data.list2[[14]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[15]] <- subset(data.list2[[15]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[16]] <- subset(data.list2[[16]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[17]] <- subset(data.list2[[17]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[18]] <- subset(data.list2[[18]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[19]] <- subset(data.list2[[19]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[20]] <- subset(data.list2[[20]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[21]] <- subset(data.list2[[21]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[22]] <- subset(data.list2[[22]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[23]] <- subset(data.list2[[23]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[24]] <- subset(data.list2[[24]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[25]] <- subset(data.list2[[25]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[26]] <- subset(data.list2[[26]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[27]] <- subset(data.list2[[27]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[28]] <- subset(data.list2[[28]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[29]] <- subset(data.list2[[29]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[30]] <- subset(data.list2[[30]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[31]] <- subset(data.list2[[31]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[32]] <- subset(data.list2[[32]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[33]] <- subset(data.list2[[33]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[34]] <- subset(data.list2[[34]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[35]] <- subset(data.list2[[35]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[36]] <- subset(data.list2[[36]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[37]] <- subset(data.list2[[37]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[38]] <- subset(data.list2[[38]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[39]] <- subset(data.list2[[39]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[40]] <- subset(data.list2[[40]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[41]] <- subset(data.list2[[41]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[42]] <- subset(data.list2[[42]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[43]] <- subset(data.list2[[43]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[44]] <- subset(data.list2[[44]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[45]] <- subset(data.list2[[45]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[46]] <- subset(data.list2[[46]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[47]] <- subset(data.list2[[47]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)
data.list2[[48]] <- subset(data.list2[[48]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_MT < 5)


################# Trying a merge instead of integrate. Find integration pipeline below this one ################


long_list <- ls(pattern = "^MK434")


big_seurat_obj <- merge(MK434_001, y = data.list2)

big_seurat_obj <- NormalizeData(big_seurat_obj)
big_seurat_obj <- FindVariableFeatures(big_seurat_obj, selection.method = "vst", nfeatures = 2000)


# Run the standard workflow for visualization and clustering
big_seurat_obj <- ScaleData(big_seurat_obj, verbose = FALSE)
big_seurat_obj <- RunPCA(big_seurat_obj, verbose = FALSE)

ElbowPlot(big_seurat_obj, ndims = 50, reduction = "pca") # based on this plot, going with 20 dims

big_seurat_obj <- RunUMAP(big_seurat_obj, reduction = "pca", dims = 1:20)
big_seurat_obj <- FindNeighbors(big_seurat_obj, reduction = "pca", dims = 1:20)
big_seurat_obj <- FindClusters(big_seurat_obj, resolution = 0.5)

DimPlot(big_seurat_obj, group.by = "orig.ident", raster = FALSE)
DimPlot(big_seurat_obj, reduction = "umap", label = TRUE, raster = FALSE, repel = FALSE)
DimPlot(big_seurat_obj, group.by = "Treatment", raster = FALSE)


# Plot UMAP in 3D
yourseuratobject <- big_seurat_obj

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~seurat_clusters, 
               colors = c("lightseagreen",
                          "gray50",
                          "darkgreen",
                          "red4",
                          "red",
                          "turquoise4",
                          "black",
                          "yellow4",
                          "royalblue1",
                          "lightcyan3",
                          "peachpuff3",
                          "khaki3",
                          "gray20",
                          "orange2",
                          "royalblue4",
                          "yellow3",
                          "gray80",
                          "darkorchid1",
                          "blue3",
                          "green3",
                          "magenta3",
                          "salmon",
                          "orange"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
fig

################################################################################################################






# normalize and identify variable features for each dataset independently
data.list2 <- lapply(X = data.list2, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list2)

#Workflow for large dataset integration#

#Scale data and run PCA
data.list2 <- lapply(X= data.list2, FUN= function(x){
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = TRUE)
})


#Perform the integration
integration.anchors <- FindIntegrationAnchors(object.list = data.list2, reference = c(1,13,24,35), reduction = "rpca")
# this command creates an 'integrated' data assay
integration.combined <- IntegrateData(anchorset = integration.anchors)


integration.combined <- ScaleData(integration.combined, verbose = FALSE)
integration.combined <- RunPCA(integration.combined, verbose = FALSE)

ElbowPlot(integration.combined, ndims = 50, reduction = "pca") # based on this plot, going with 10 dims

integration.combined <- RunUMAP(integration.combined, dims = 1:10)
integration.combined <- FindNeighbors(integration.combined, reduction = "pca", dims = 1:10)
integration.combined <- FindClusters(integration.combined, resolution = 0.5) # adjusted from 0.5 to 0.1 due to number of small clusters

#Reorder timepoints
integration.combined$Timepoint <- factor(x = integration.combined$Timepoint, levels = c("6hr", "24hr"))

saveRDS(integration.combined, file = "./integration.combined_no_cell_type.rds")
#integration.combined <- readRDS("./integration.combined_no_cell_type.rds")

## Data Exploration ##
integration.combined.sub <- subset(x = integration.combined, subset = Timepoint != "BRIDGING")

png("dim_plot_all_orig_ident.png",width = 10, height = 8, units = 'in', res = 900)
DimPlot(integration.combined.sub, group.by = "orig.ident", raster = FALSE)
dev.off()

png("dim_plot_all_treatment.png",width = 50, height = 10, units = 'cm', res = 900)
DimPlot(integration.combined.sub, split.by = "Treatment", raster = FALSE)
dev.off()

png("dim_plot_all_donor.png",width = 10, height = 8, units = 'in', res = 900)
DimPlot(integration.combined.sub, split.by = "Donor", raster = FALSE)
dev.off()

png("dim_plot_all_timepoint.png",width = 10, height = 8, units = 'in', res = 900)
DimPlot(integration.combined.sub, split.by = "Timepoint", raster = FALSE)
dev.off()

p1 <- DimPlot(integration.combined.sub, reduction = "umap", split.by = "Timepoint", group.by = "Treatment", raster = FALSE)
p1

p1 <- DimPlot(integration.combined.sub, reduction = "umap", group.by = "Treatment", raster = FALSE)
p2 <- DimPlot(integration.combined.sub, reduction = "umap", label = TRUE, raster = FALSE, repel = FALSE)
p1 + p2

#Find markers for individual treatments
DefaultAssay(integration.combined) <- "RNA"
mock.markers <- FindMarkers(integration.combined, ident.1 = "Mock", group.by = "Treatment")
dmso.markers <- FindMarkers(integration.combined, ident.1 = "DMSO", group.by = "Treatment")
DA_5nm.markers <- FindMarkers(integration.combined, ident.1 = "DiABZi 5nM", group.by = "Treatment")
DA_25nm.markers <- FindMarkers(integration.combined, ident.1 = "DiABZi 25nM", group.by = "Treatment")
ST002_25nm.markers <- FindMarkers(integration.combined, ident.1 = "ST002 25uM", group.by = "Treatment")
ST002_50nm.markers <- FindMarkers(integration.combined, ident.1 = "ST002 50uM", group.by = "Treatment")
ST012_25nm.markers <- FindMarkers(integration.combined, ident.1 = "ST012 25uM", group.by = "Treatment")
ST012_50nm.markers <- FindMarkers(integration.combined, ident.1 = "ST012 50uM", group.by = "Treatment")
ST020_25nm.markers <- FindMarkers(integration.combined, ident.1 = "ST020 25uM", group.by = "Treatment")
ST020_50nm.markers <- FindMarkers(integration.combined, ident.1 = "ST020 50uM", group.by = "Treatment")
INI3069_25uM.markers <- FindMarkers(integration.combined, ident.1 = "INI3069 25uM", group.by = "Treatment")
INI3069_50uM.markers <- FindMarkers(integration.combined, ident.1 = "INI3069 50uM", group.by = "Treatment")
treat_1148_10uM.markers <- FindMarkers(integration.combined, ident.1 = "1148 10uM", group.by = "Treatment")
BRIDGING.markers <- FindMarkers(integration.combined, ident.1 = "BRIDGING", group.by = "Treatment")

write.csv(mock.markers, "mock.markers.csv")
write.csv(dmso.markers, "dmso.markers.csv")
write.csv(DA_5nm.markers, "DA_5nm.markers.csv")
write.csv(DA_25nm.markers, "DA_25nm.markers.csv")
write.csv(ST002_25nm.markers, "ST002_25nm.markers.csv")
write.csv(ST002_50nm.markers, "ST002_50nm.markers.csv")
write.csv(ST012_25nm.markers, "ST012_25nm.markers.csv")
write.csv(ST012_50nm.markers, "ST012_50nm.markers.csv")
write.csv(ST020_25nm.markers, "ST020_25nm.markers.csv")
write.csv(ST020_50nm.markers, "ST020_50nm.markers.csv")
write.csv(INI3069_25uM.markers, "INI3069_25uM.markers.csv")
write.csv(INI3069_50uM.markers, "INI3069_50uM.markers.csv")
write.csv(treat_1148_10uM.markers, "treat_1148_10uM.markers.csv")
write.csv(BRIDGING.markers, "BRIDGING.markers.csv")

#Find markers for individual time points
DefaultAssay(integration.combined) <- "RNA"
six_hr.markers <- FindMarkers(integration.combined, ident.1 = "6hr", group.by = "Timepoint")
twentyfour_hr.markers <- FindMarkers(integration.combined, ident.1 = "24hr", group.by = "Timepoint")

write.csv(six_hr.markers, "6hr.markers.csv")
write.csv(twentyfour_hr.markers, "24hr.markers.csv")

#Find markers for individual donors
DefaultAssay(integration.combined) <- "RNA"
HAQ2.markers <- FindMarkers(integration.combined, ident.1 = "HAQ2", group.by = "Donor")
WT5.markers <- FindMarkers(integration.combined, ident.1 = "WT5", group.by = "Donor")

write.csv(HAQ2.markers, "HAQ2.markers.csv")
write.csv(WT5.markers, "WT5.markers.csv")

#Perform DE between groups
dmso.mock.markers <- FindMarkers(integration.combined, ident.1 = "DMSO", ident.2 = "Mock", group.by = "Treatment")
DA5nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "DiABZi 5nM", ident.2 = "DMSO", group.by = "Treatment")
DA25nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")

write.csv(dmso.mock.markers, "./dmso_vs_mock_markers.csv")
write.csv(DA5nm.dmso.markers, "./DA5nm_vs_dmso_markers.csv")
write.csv(DA25nm.dmso.markers, "./DA25nm_vs_dmso_markers.csv")

ST002_25nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST002 25uM", ident.2 = "DMSO", group.by = "Treatment")
ST002_50nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST002 50uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_25nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST012 25uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_50nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST012 50uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_25nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST020 25uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_50nm.dmso.markers <- FindMarkers(integration.combined, ident.1 = "ST020 50uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_25uM.dmso.markers <- FindMarkers(integration.combined, ident.1 = "INI3069 25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50uM.dmso.markers <- FindMarkers(integration.combined, ident.1 = "INI3069 50uM", ident.2 = "DMSO", group.by = "Treatment")
treat_1148_10uM.dmso.markers <- FindMarkers(integration.combined, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")
BRIDGING.dmso.markers <- FindMarkers(integration.combined, ident.1 = "BRIDGING", ident.2 = "DMSO", group.by = "Treatment")


write.csv(ST002_25nm.dmso.markers, "ST002_25nm.dmso_markers.csv")
write.csv(ST002_50nm.dmso.markers, "ST002_50nm.dmso_markers.csv")
write.csv(ST012_25nm.dmso.markers, "ST012_25nm.dmso_markers.csv")
write.csv(ST012_25nm.dmso.markers, "ST012_25nm.dmso_markers.csv")
write.csv(ST020_25nm.dmso.markers, "ST020_25nm.dmso_markers.csv")
write.csv(ST020_25nm.dmso.markers, "ST020_25nm.dmso_markers.csv")
write.csv(INI3069_25uM.dmso.markers, "INI3069_25uM.dmso_markers.csv")
write.csv(INI3069_50uM.dmso.markers, "INI3069_50uM.dmso_markers.csv")
write.csv(treat_1148_10uM.dmso.markers, "treat_1148_10uM.dmso_markers.csv")
write.csv(BRIDGING.dmso.markers, "BRIDGING.dmso_markers.csv")

# Generate more figures
integration.combined.no.umap <- integration.combined.sub
integration.combined.no.umap[["umap"]] <- NULL
DimPlot(integration.combined.no.umap, raster = FALSE) + RotatedAxis()

# Get counts of genes 
mean(integration.combined@meta.data$nFeature_RNA)

# Plot UMAP in 3D
yourseuratobject <- integration.combined.sub

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~seurat_clusters, 
               colors = c("lightseagreen",
                          "gray50",
                          "darkgreen",
                          "red4",
                          "red",
                          "turquoise4",
                          "black",
                          "yellow4",
                          "royalblue1",
                          "lightcyan3",
                          "peachpuff3",
                          "khaki3",
                          "gray20",
                          "orange2",
                          "royalblue4",
                          "yellow3",
                          "gray80",
                          "darkorchid1"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# Lets do some cell type calls! WOOT


# First need to subset the integrated object to get just the baseline
integration.combined.baseline <- subset(integration.combined.sub, subset = Treatment %in% c("DMSO", "Mock"))
levels(integration.combined.baseline)

# Find all markers for different clusters
cluster0.markers <- FindMarkers(integration.combined.baseline, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(integration.combined.baseline, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(integration.combined.baseline, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(integration.combined.baseline, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(integration.combined.baseline, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(integration.combined.baseline, ident.1 = 5, min.pct = 0.25)
cluster6.markers <- FindMarkers(integration.combined.baseline, ident.1 = 6, min.pct = 0.25)
cluster7.markers <- FindMarkers(integration.combined.baseline, ident.1 = 7, min.pct = 0.25)
cluster8.markers <- FindMarkers(integration.combined.baseline, ident.1 = 8, min.pct = 0.25)
cluster9.markers <- FindMarkers(integration.combined.baseline, ident.1 = 9, min.pct = 0.25)
cluster10.markers <- FindMarkers(integration.combined.baseline, ident.1 = 10, min.pct = 0.25)
cluster11.markers <- FindMarkers(integration.combined.baseline, ident.1 = 11, min.pct = 0.25)
cluster12.markers <- FindMarkers(integration.combined.baseline, ident.1 = 12, min.pct = 0.25)
cluster13.markers <- FindMarkers(integration.combined.baseline, ident.1 = 13, min.pct = 0.25)
cluster14.markers <- FindMarkers(integration.combined.baseline, ident.1 = 14, min.pct = 0.25)
cluster15.markers <- FindMarkers(integration.combined.baseline, ident.1 = 15, min.pct = 0.25)
cluster16.markers <- FindMarkers(integration.combined.baseline, ident.1 = 16, min.pct = 0.25)
cluster17.markers <- FindMarkers(integration.combined.baseline, ident.1 = 17, min.pct = 0.25)
cluster18.markers <- FindMarkers(integration.combined.baseline, ident.1 = 18, min.pct = 0.25)


write.csv(cluster0.markers, "./cluster0.markers.csv", row.names = TRUE)
write.csv(cluster1.markers, "./cluster1.markers.csv", row.names = TRUE)
write.csv(cluster2.markers, "./cluster2.markers.csv", row.names = TRUE)
write.csv(cluster3.markers, "./cluster3.markers.csv", row.names = TRUE)
write.csv(cluster4.markers, "./cluster4.markers.csv", row.names = TRUE)
write.csv(cluster5.markers, "./cluster5.markers.csv", row.names = TRUE)
write.csv(cluster6.markers, "./cluster6.markers.csv", row.names = TRUE)
write.csv(cluster7.markers, "./cluster7.markers.csv", row.names = TRUE)
write.csv(cluster8.markers, "./cluster8.markers.csv", row.names = TRUE)
write.csv(cluster9.markers, "./cluster9.markers.csv", row.names = TRUE)
write.csv(cluster10.markers, "./cluster10.markers.csv", row.names = TRUE)
write.csv(cluster11.markers, "./cluster11.markers.csv", row.names = TRUE)
write.csv(cluster12.markers, "./cluster12.markers.csv", row.names = TRUE)
write.csv(cluster13.markers, "./cluster13.markers.csv", row.names = TRUE)
write.csv(cluster14.markers, "./cluster14.markers.csv", row.names = TRUE)
write.csv(cluster15.markers, "./cluster15.markers.csv", row.names = TRUE)
write.csv(cluster16.markers, "./cluster16.markers.csv", row.names = TRUE)
write.csv(cluster17.markers, "./cluster17.markers.csv", row.names = TRUE)
write.csv(cluster18.markers, "./cluster18.markers.csv", row.names = TRUE)

# Going to try clustermole for automated cell type assignments
# it can be found here: https://igordot.github.io/clustermole/index.html




avg_exp_mat <- AverageExpression(integration.combined.baseline)
avg_exp_mat <- as.matrix(avg_exp_mat$RNA)
avg_exp_mat <- log1p(avg_exp_mat)

enrich_tbl <- clustermole_enrichment(expr_mat = avg_exp_mat, species = "hs")
cluster0 <- enrich_tbl %>% filter(cluster == "0") %>% head(50)
cluster1 <- enrich_tbl %>% filter(cluster == "1") %>% head(50)
cluster2 <- enrich_tbl %>% filter(cluster == "2") %>% head(50)
cluster3 <- enrich_tbl %>% filter(cluster == "3") %>% head(50)
cluster4 <- enrich_tbl %>% filter(cluster == "4") %>% head(50)
cluster5 <- enrich_tbl %>% filter(cluster == "5") %>% head(50)
cluster6 <- enrich_tbl %>% filter(cluster == "6") %>% head(50)
cluster7 <- enrich_tbl %>% filter(cluster == "7") %>% head(50)
cluster8 <- enrich_tbl %>% filter(cluster == "8") %>% head(50)
cluster9 <- enrich_tbl %>% filter(cluster == "9") %>% head(50)
cluster10 <- enrich_tbl %>% filter(cluster == "10") %>% head(50)
cluster11 <- enrich_tbl %>% filter(cluster == "11") %>% head(50)
cluster12 <- enrich_tbl %>% filter(cluster == "12") %>% head(50)
cluster13 <- enrich_tbl %>% filter(cluster == "13") %>% head(50)
cluster14 <- enrich_tbl %>% filter(cluster == "14") %>% head(50)
cluster15 <- enrich_tbl %>% filter(cluster == "15") %>% head(50)
cluster16 <- enrich_tbl %>% filter(cluster == "16") %>% head(50)
cluster17 <- enrich_tbl %>% filter(cluster == "17") %>% head(50)
cluster18 <- enrich_tbl %>% filter(cluster == "18") %>% head(50)

write.csv(cluster0, "./cluster0.csv", row.names = TRUE)
write.csv(cluster1, "./cluster1.csv", row.names = TRUE)
write.csv(cluster2, "./cluster2.csv", row.names = TRUE)
write.csv(cluster3, "./cluster3.csv", row.names = TRUE)
write.csv(cluster4, "./cluster4.csv", row.names = TRUE)
write.csv(cluster5, "./cluster5.csv", row.names = TRUE)
write.csv(cluster6, "./cluster6.csv", row.names = TRUE)
write.csv(cluster7, "./cluster7.csv", row.names = TRUE)
write.csv(cluster8, "./cluster8.csv", row.names = TRUE)
write.csv(cluster9, "./cluster9.csv", row.names = TRUE)
write.csv(cluster10, "./cluster10.csv", row.names = TRUE)
write.csv(cluster11, "./cluster11.csv", row.names = TRUE)
write.csv(cluster12, "./cluster12.csv", row.names = TRUE)
write.csv(cluster13, "./cluster13.csv", row.names = TRUE)
write.csv(cluster14, "./cluster14.csv", row.names = TRUE)
write.csv(cluster15, "./cluster15.csv", row.names = TRUE)
write.csv(cluster16, "./cluster16.csv", row.names = TRUE)
write.csv(cluster17, "./cluster17.csv", row.names = TRUE)
write.csv(cluster18, "./cluster18.csv", row.names = TRUE)

# trying another method: sc-type
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = integration.combined.baseline[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(integration.combined.baseline@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(integration.combined.baseline@meta.data[integration.combined.baseline@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integration.combined.baseline@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

integration.combined.baseline@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  integration.combined.baseline@meta.data$customclassif[integration.combined.baseline@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(integration.combined.baseline, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', raster = FALSE) 



#### Next up, SingleR #### Didn't work, I'll have to come back to this

hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

pred.hesc <- SingleR(test = GetAssayData(integration.combined), ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)

#### Verify predicted cell type assignments using canonical marker genes ####
DefaultAssay(integration.combined.baseline) <- "RNA"
tcells <- FeaturePlot(integration.combined.baseline, features = "CD3D")
tcells

cd8_t <- FeaturePlot(integration.combined.baseline, features = c("CD8A", "CD8B"))
cd8_t

cd4_t <- FeaturePlot(integration.combined.baseline, features = c("CD4"))
cd4_t

cd8_effector <- FeaturePlot(integration.combined.baseline, features = c("GZMK", "GZMH", "PRF1", "CCL5"))
cd8_effector

cd8_memory <- FeaturePlot(integration.combined.baseline, features = c("ITGB1"))
cd8_memory

cd8_naive <- FeaturePlot(integration.combined.baseline, features = c("CCR7"))
cd8_naive

cd4_naive <- FeaturePlot(integration.combined.baseline, features = c("IL7R", "CCR7"))
cd4_naive

cd4_memory <- FeaturePlot(integration.combined.baseline, features = c("IL7R", "S100A4"))
cd4_memory

treg <- FeaturePlot(integration.combined.baseline, features = c("FOXP3", "IL2RA"))
treg

MAIT <- FeaturePlot(integration.combined.baseline, features = c("SLC4A10", "TRAV1-2"))
integration.combined.baseline$RNA@data["TRAV1-2",1]=1
MAIT

gamma_delta <- FeaturePlot(integration.combined.baseline, features = c("TRGV9", "TRDV2"))
integration.combined.baseline$RNA@data["TRGV9",1]=1
integration.combined.baseline$RNA@data["TRDV2",1]=1
gamma_delta

b_cells <- FeaturePlot(integration.combined.baseline, features = c("CD79A", "CD79B"))
b_cells

b_naive <- FeaturePlot(integration.combined.baseline, features = c("CD27"))
b_naive

b_plasma <- FeaturePlot(integration.combined.baseline, features = c("SDC1", "MZB1", "XBP1"))
b_plasma

mono_classic <- FeaturePlot(integration.combined.baseline, features = c("CD14", "LYZ"))
mono_classic

mono_nonclassic <- FeaturePlot(integration.combined.baseline, features = c("FCGR3A", "MS4A7"))
mono_nonclassic

dcs <- FeaturePlot(integration.combined.baseline, features = c("IL3RA", "CLEC4C", "CST3", "NRP1"))
dcs

dc_myeloid <- FeaturePlot(integration.combined.baseline, features = c("FCER1A", "CD1C"))
dc_myeloid

dc_plasmacytoid <- FeaturePlot(integration.combined.baseline, features = c("FCER1A", "LILRA4"))
dc_plasmacytoid

nk_cells <- FeaturePlot(integration.combined.baseline, features = c("KLRB1", "KLRC1", "KLRD1", "GNLY", "NKG7", "NCAM1"))
nk_cells

platelets <- FeaturePlot(integration.combined.baseline, features = c("PPBP"))
platelets


FeaturePlot(integration.combined.baseline, features = c("IL3RA", "CCR3"), raster = FALSE)
FeaturePlot(integration.combined.baseline, features = c("MMP8"), raster = FALSE)


#### Cell Type Assignments ####

# New cluster assignments
new.cluster.ids <- c("CD4+ T", "Naive T", "Monocytes", "CD8+ T", "Monocytes", "Monocytes", "CD4+ T", "CD4+ T", "CD4+ T", "B", "NK", "Naive B", "CD8+ T", "NK", "CD8+ T", "B", "Stress Induced Monocytes", "Myeloid dendritic", "Basophils")
names(new.cluster.ids) <- levels(integration.combined)
integration.combined <- RenameIdents(integration.combined, new.cluster.ids)
integration.combined@meta.data$celltype_clusters <- Idents(integration.combined)

plot1 <- DimPlot(integration.combined, reduction = "umap", label = TRUE, pt.size = 0.5, raster = FALSE, repel = TRUE)
plot1

plot1 <- FeaturePlot(integration.combined, features = c("CD4"), order = TRUE, reduction = "umap", label = TRUE, raster = FALSE)
plot1

plot1 <- FeaturePlot(integration.combined, features = c("CD8A"), order = TRUE, reduction = "umap", label = TRUE, raster = FALSE)
plot1

#### Figures for presentation ####
plot1 <- DimPlot(integration.combined.sub, reduction = "umap", label = TRUE, pt.size = 0.5, raster = FALSE, repel = TRUE, group.by = "Treatment")
plot1

plot1 <- DimPlot(integration.combined.sub, reduction = "umap", label = TRUE, pt.size = 0.5, raster = FALSE, group.by = "ident", split.by = "Treatment")
plot1

integration.combined.wt_only <- subset(x = integration.combined, subset = Donor == "WT5")

integration.combined.wt_only.6hr <- subset(x = integration.combined.wt_only, subset = Timepoint == "6hr")
integration.combined.wt_only.24hr <- subset(x = integration.combined.wt_only, subset = Timepoint == "24hr")

plot1 <- DimPlot(integration.combined.wt_only.6hr, reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, group.by = "ident", split.by = "Treatment", ncol = 3)
plot1

plot1 <- DimPlot(integration.combined.wt_only.24hr, reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, group.by = "ident", split.by = "Treatment", ncol = 3)
plot1


plot1 <- DimPlot(integration.combined.wt_only.6hr, reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, group.by = "Treatment", split.by = "ident", ncol = 3)
plot1


plot1 <- DimPlot(subset(integration.combined.wt_only.6hr, subset = Treatment %in% c("DMSO", "DiABZi 5nM", "1148 10uM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1$data$Treatment <- factor(x = plot1$data$Treatment, levels = c("DMSO", "DiABZi 5nM", "1148 10uM"))
plot1

plot1 <- DimPlot(subset(integration.combined.wt_only.24hr, subset = Treatment %in% c("DMSO", "DiABZi 5nM", "1148 10uM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1$data$Treatment <- factor(x = plot1$data$Treatment, levels = c("DMSO", "DiABZi 5nM", "1148 10uM"))
plot1



plot1 <- DimPlot(subset(integration.combined.wt_only.6hr, subset = Treatment %in% c("1148 10uM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1

plot1 <- DimPlot(subset(integration.combined.wt_only.24hr, subset = Treatment %in% c("1148 10uM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1



plot1 <- DimPlot(subset(integration.combined.wt_only.6hr, subset = Treatment %in% c("DiABZi 5nM", "DiABZi 25nM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1$data$Treatment <- factor(x = plot1$data$Treatment, levels = c("DiABZi 5nM", "DiABZi 25nM"))
plot1

plot1 <- DimPlot(subset(integration.combined.wt_only.24hr, subset = Treatment %in% c("DiABZi 5nM", "DiABZi 25nM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1$data$Treatment <- factor(x = plot1$data$Treatment, levels = c("DiABZi 5nM", "DiABZi 25nM"))
plot1



plot1 <- DimPlot(subset(integration.combined.wt_only.6hr, subset = Treatment %in% c("DMSO")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1

plot1 <- DimPlot(subset(integration.combined.wt_only.24hr, subset = Treatment %in% c("DMSO", "Mock")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, split.by = "Treatment")
plot1



integration.combined.haq2_only <- subset(x = integration.combined, subset = Donor == "HAQ2")
integration.combined.haq2_only.6hr <- subset(x = integration.combined.haq2_only, subset = Timepoint == "6hr")
integration.combined.haq2_only.24hr <- subset(x = integration.combined.haq2_only, subset = Timepoint == "24hr")

integration.combined.6hr <- subset(x = integration.combined, subset = Timepoint == "6hr")
plot1 <- DimPlot(subset(integration.combined.6hr, subset = Treatment %in% c("DMSO", "DiABZi 25nM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, group.by = "Treatment", split.by =  "Donor")
plot1

integration.combined.24hr <- subset(x = integration.combined, subset = Timepoint == "24hr")
plot1 <- DimPlot(subset(integration.combined.24hr, subset = Treatment %in% c("DMSO", "DiABZi 25nM")), reduction = "umap", label = FALSE, pt.size = 0.5, raster = FALSE, group.by = "Treatment", split.by =  "Donor")
plot1

plot1 <- VlnPlot(subset(integration.combined.wt_only, downsample = 500), features = "CXCL10", idents = c("CD4+ T"), split.by = "Treatment", raster = TRUE, pt.size = 0, add.noise = T)
plot1


### Time for UCell ###
markers <- list()
str <- "IFIT2
HERC5
OASL
RSAD2
IFIT1
PMAIP1
IFIT3
ISG15
CMPK2
IFIH1
ZC3HAV1
MX1
EPSTI1
OAS2
OAS3
OAS1
AFAP1
CXCL10
ISG20
HERC6
IFI6
TNFSF10
DDX58
IFI44
IFI44L
CD38
DHX58
DDX60
PARP14
RGS1
ZBP1
GBP5
HELB
GBP4
GBP1
SRPK2
TNIK
USP18
IRF7
FAS
CCL3
NFKBIZ
MX2
SOCS1
EIF2AK2
RNF149
SAT1
RIPK1
RIPOR2
SLAMF7"

markers$DA5nm_up <- as.list(strsplit(str, "\n")[[1]])

str <-"IL32
KLF2
LYZ
AKNA
S100A9
THBS1
VCAN
GLUL
GNLY
MYO1F
TRBC1
AHNAK
CD6
CST7
PLXND1
SH2D2A
ITGAL
SUN2
ZFP36L2
S100A8
SEPTIN9
SELPLG
TCF7
S1PR4
CORO1A
CD5
RASSF5
ARHGAP45
EFHD2
ATP2B4
FAM102A
CSF1R
ARHGDIB
FMNL1
NLRC3
SESN3
ARHGEF1
LIMD2
SASH3
ZAP70
BCL9L
TNFRSF1B
SH3BGRL3
TRABD
ITGB2
SBNO2
ENO1
GRK2
S1PR1
ARRB2"

markers$DA5nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "IFIT2
HERC5
OASL
PMAIP1
RSAD2
IFIT3
IFIT1
ISG15
IFIH1
ZC3HAV1
CMPK2
OAS3
OAS1
AFAP1
MX1
OAS2
EPSTI1
CXCL10
IFI44L
DDX58
IFI44
IFNB1
ISG20
DDX60
TNFSF10
PARP14
IFI6
DHX58
GBP4
HERC6
EIF2AK2
IRF7
HELB
CD38
GCH1
SAT1
USP18
SAMD9L
SRPK2
NFKBIZ
GBP5
HSPA1B
MX2
GBP1
RGS1
PPP1R15A
TNIK
IRF8
RIPOR2
SLAMF7"

markers$DA25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "LYZ
IL32
S100A9
S100A8
KLF2
MYO1F
THBS1
AKNA
TCF7
VCAN
TRBC1
FAM102A
SUN2
CD5
CST7
GLUL
SELPLG
ZFP36L2
RASSF5
CD6
GNLY
SESN3
ARHGEF1
ITGAL
TRAC
AHNAK
S1PR1
ARHGDIB
BCL9L
CD4
EFHD2
SEPTIN9
PLXND1
YPEL3
ATP2B4
RBL2
CSF1R
ARHGAP45
NLRC3
ZAP70
ITGB2
SH2D2A
S1PR4
NLRP1
ARRB2
FXYD5
FGFBP2
TRABD
SORL1
SBNO2
"

markers$DA25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "ISG15
IFI44L
CXCL8
RSAD2
CMPK2
MX1
OAS1
IDO1
IFI44
XAF1
IFIT1
IFIT3
CTSL
OAS2
IFI6
TRIM22
EPSTI1
MX2
OAS3
LY6E
MT-ND3
IRF7
PLAUR
GBP5
OASL
SLAMF7
HERC5
GBP1
SAMD9L
BASP1
C15orf48
HERC6
UBE2L6
MT2A
EIF2AK2
DDX60
MOV10
HELZ2
PARP9
DEPRECATED-ENSG00000205413
SAMD9
SNX10
IL1B
ZBP1
PARP12
KIAA0319L
CREG1
DTX3L
DHX58
APOL6"

markers$KIN1148_up <- as.list(strsplit(str, "\n")[[1]])

str <- "THBS1
LYZ
VCAN
CD14
S100A9
S100A8
CXCL5
CYBB
MPEG1
CSF1R
EEF1G
JAML
KLRG1
TOMM7
MT-ND4L
KLRB1
AIF1
CD52
SH3BGRL3
MS4A7
ITGB2
ANPEP
CD48
CXCL16
ENG
SPI1
PECAM1
SEPTIN9
EVI2A
AHNAK
ALOX5
KLRK1
LY9
LTA4H
CSF2RA
IFNGR2
IDS
ENO1
LAPTM5
EPB41
HK2
EEF1B2
ACTG1
TOB1
KCNAB2
VNN2
S100A4
PTPRC
UHMK1
USP53"

markers$KIN1148_down <- as.list(strsplit(str, "\n")[[1]])

integration.combined <- AddModuleScore_UCell(integration.combined, features = markers)
signature.names <- paste0(names(markers), "_UCell")


VlnPlot(integration.combined, features = signature.names, group.by = "Treatment", raster = F, pt.size = 0)
FeaturePlot(integration.combined, order = TRUE, reduction = "umap", features = signature.names, label = F, raster = F) &
  scale_color_viridis_c()

markers_heatmap <- unlist(unique(markers))
DoHeatmap(integration.combined, features = markers_heatmap, cells = 1:250, size = 3) + theme(axis.text.y = element_text(size = 5))


### Cell Type Distribution ###
ggplot(integration.combined@meta.data, aes(celltype_clusters))+geom_bar(stat="count") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggplot(integration.combined.wt_only@meta.data, aes(celltype_clusters))+geom_bar(stat="count") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggplot(integration.combined.haq2_only@meta.data, aes(celltype_clusters))+geom_bar(stat="count") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggplot(integration.combined.baseline@meta.data, aes(celltype_clusters))+geom_bar(stat="count") +
  scale_x_discrete(guide = guide_axis(angle = 45))


table(integration.combined@meta.data$celltype_clusters)
asdf <- count_cells(seurat_obj = integration.combined, group_by_var = "celltype_clusters")
write.csv(asdf, "./cell_type_percents.csv")

asdf <- count_cells(seurat_obj = integration.combined.wt_only, group_by_var = "celltype_clusters")

asdf <- count_cells(seurat_obj = integration.combined.haq2_only, group_by_var = "celltype_clusters")

asdf <- count_cells(seurat_obj = integration.combined.baseline, group_by_var = "celltype_clusters")

asdf <- count_cells(seurat_obj = integration.combined, group_by_var = "Condition", subgroup_var = "celltype_clusters")
write.csv(asdf, "./cell_types_per_condition.csv")







### Heatmaps for DE genes for each cell type ###
integration.combined.6hr.cd4t.wt <- subset(integration.combined.6hr, idents = "CD4+ T", subset = Donor == "WT5")
DA25nm.dmso.markers.6.cd4.wt <- FindMarkers(integration.combined.6hr.cd4t.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.cd4.wt <- FindMarkers(integration.combined.6hr.cd4t.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.6hr.mono.wt <- subset(integration.combined.6hr, idents = "Monocytes", subset = Donor == "WT5")
DA25nm.dmso.markers.6.mono.wt <- FindMarkers(integration.combined.6hr.mono.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.mono.wt <- FindMarkers(integration.combined.6hr.mono.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.6hr.cd8t.wt <- subset(integration.combined.6hr, idents = "CD8+ T", subset = Donor == "WT5")
DA25nm.dmso.markers.6.cd8t.wt <- FindMarkers(integration.combined.6hr.cd8t.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.cd8t.wt <- FindMarkers(integration.combined.6hr.cd8t.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.6hr.b.wt <- subset(integration.combined.6hr, idents = "B", subset = Donor == "WT5")
DA25nm.dmso.markers.6.b.wt <- FindMarkers(integration.combined.6hr.b.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.b.wt <- FindMarkers(integration.combined.6hr.b.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.6hr.nk.wt <- subset(integration.combined.6hr, idents = "NK", subset = Donor == "WT5")
DA25nm.dmso.markers.6.nk.wt <- FindMarkers(integration.combined.6hr.nk.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.nk.wt <- FindMarkers(integration.combined.6hr.nk.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.6hr.md.wt <- subset(integration.combined.6hr, idents = "Myeloid dendritic", subset = Donor == "WT5")
DA25nm.dmso.markers.6.md.wt <- FindMarkers(integration.combined.6hr.md.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.6.md.wt <- FindMarkers(integration.combined.6hr.md.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)


# 6 hr all cell types
list_of_data <- list(DA25nm.dmso.markers.6.cd4.wt, DA25nm.dmso.markers.6.cd8t.wt, DA25nm.dmso.markers.6.nk.wt, DA25nm.dmso.markers.6.mono.wt, DA25nm.dmso.markers.6.b.wt)

DA25_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

DA25_combined[is.na(DA25_combined)] <- 0

DA25_combined2 <- DA25_combined[,-1]

rownames(DA25_combined2) <- DA25_combined[,1]

head(DA25_combined2)

hmap <- heatmap.L.4(as.matrix(DA25_combined2),
            colcolorlist = colcolormatrix,
            distmethod = "euclidean",
            clustermethod='ward.D2',
            cutoffmethod = 'depth',
            cutoff=0
            )

###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/6hr_25nm_diabzi")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(DA25_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")




integration.combined.24hr.cd4t.wt <- subset(integration.combined.24hr, idents = "CD4+ T", subset = Donor == "WT5")
DA25nm.dmso.markers.24.cd4.wt <- FindMarkers(integration.combined.24hr.cd4t.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.cd4.wt <- FindMarkers(integration.combined.24hr.cd4t.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.24hr.mono.wt <- subset(integration.combined.24hr, idents = "Monocytes", subset = Donor == "WT5")
DA25nm.dmso.markers.24.mono.wt <- FindMarkers(integration.combined.24hr.mono.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.mono.wt <- FindMarkers(integration.combined.24hr.mono.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.24hr.cd8t.wt <- subset(integration.combined.24hr, idents = "CD8+ T", subset = Donor == "WT5")
DA25nm.dmso.markers.24.cd8t.wt <- FindMarkers(integration.combined.24hr.cd8t.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.cd8t.wt <- FindMarkers(integration.combined.24hr.cd8t.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.24hr.b.wt <- subset(integration.combined.24hr, idents = "B", subset = Donor == "WT5")
DA25nm.dmso.markers.24.b.wt <- FindMarkers(integration.combined.24hr.b.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.b.wt <- FindMarkers(integration.combined.24hr.b.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.24hr.nk.wt <- subset(integration.combined.24hr, idents = "NK", subset = Donor == "WT5")
DA25nm.dmso.markers.24.nk.wt <- FindMarkers(integration.combined.24hr.nk.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.nk.wt <- FindMarkers(integration.combined.24hr.nk.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")

integration.combined.24hr.md.wt <- subset(integration.combined.24hr, idents = "Myeloid dendritic", subset = Donor == "WT5")
DA25nm.dmso.markers.24.md.wt <- FindMarkers(integration.combined.24hr.md.wt, ident.1 = "DiABZi 25nM", ident.2 = "DMSO", group.by = "Treatment")
kin1148.dmso.markers.24.md.wt <- FindMarkers(integration.combined.24hr.md.wt, ident.1 = "1148 10uM", ident.2 = "DMSO", group.by = "Treatment")


colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)

# CD4 T cells
list_of_data <- list(DA25nm.dmso.markers.6.cd4.wt, DA25nm.dmso.markers.24.cd4.wt, kin1148.dmso.markers.6.cd4.wt, kin1148.dmso.markers.24.cd4.wt)

cd4_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

cd4_combined[is.na(cd4_combined)] <- 0

cd4_combined2 <- cd4_combined[,-1]

rownames(cd4_combined2) <- cd4_combined[,1]

hmap <- heatmap.L.4(as.matrix(cd4_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "spearman",
                    clustermethod='ward.D2',
                    cutoffmethod = 'depth'
                    #cutoff = 1
                    )


###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/cd4t")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cd4_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


# CD8 T cells
list_of_data <- list(DA25nm.dmso.markers.6.cd8t.wt, DA25nm.dmso.markers.24.cd8t.wt, kin1148.dmso.markers.6.cd8t.wt, kin1148.dmso.markers.24.cd8t.wt)

cd8_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

cd8_combined[is.na(cd8_combined)] <- 0

cd8_combined2 <- cd8_combined[,-1]

rownames(cd8_combined2) <- cd8_combined[,1]

hmap <- heatmap.L.4(as.matrix(cd8_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "spearman",   #choices are 'euclidean', 'spearman', 'pearson'
                    clustermethod='ward.D2',   #choices are 'ward.D2' (default), 'ward.D', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
                    cutoffmethod = 'depth',    #choices are 'depth' (default), 'height', 'number'
                    cutoff = 2                 #default is '3'
)

###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/cd8t")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cd8_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


# NK cells
list_of_data <- list(DA25nm.dmso.markers.6.nk.wt, DA25nm.dmso.markers.24.nk.wt, kin1148.dmso.markers.6.nk.wt, kin1148.dmso.markers.24.nk.wt)

nk_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

nk_combined[is.na(nk_combined)] <- 0

nk_combined2 <- nk_combined[,-1]

rownames(nk_combined2) <- nk_combined[,1]

hmap <- heatmap.L.4(as.matrix(nk_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "euclidean",   #choices are 'euclidean', 'spearman', 'pearson'
                    clustermethod='ward.D2',   #choices are 'ward.D2' (default), 'ward.D', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
                    cutoffmethod = 'depth',    #choices are 'depth' (default), 'height', 'number'
                    cutoff = 0                 #default is '3'
)

###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/nk")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(nk_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")

# Monocytes
list_of_data <- list(DA25nm.dmso.markers.6.mono.wt, DA25nm.dmso.markers.24.mono.wt, kin1148.dmso.markers.6.mono.wt, kin1148.dmso.markers.24.mono.wt)

mono_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

mono_combined[is.na(mono_combined)] <- 0

mono_combined2 <- mono_combined[,-1]

rownames(mono_combined2) <- mono_combined[,1]

hmap <- heatmap.L.4(as.matrix(mono_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "euclidean",   #choices are 'euclidean', 'spearman', 'pearson'
                    clustermethod='ward.D',   #choices are 'ward.D2' (default), 'ward.D', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
                    cutoffmethod = 'number',    #choices are 'depth' (default), 'height', 'number'
                    cutoff = 5                 #default is '3'
)


###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/mono")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(mono_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


# B cells
list_of_data <- list(DA25nm.dmso.markers.6.b.wt, DA25nm.dmso.markers.24.b.wt, kin1148.dmso.markers.6.b.wt, kin1148.dmso.markers.24.b.wt)

b_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

b_combined[is.na(b_combined)] <- 0

b_combined2 <- b_combined[,-1]

rownames(b_combined2) <- b_combined[,1]

hmap <- heatmap.L.4(as.matrix(b_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "pearson",   #choices are 'euclidean', 'spearman', 'pearson'
                    clustermethod='ward.D',   #choices are 'ward.D2' (default), 'ward.D', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
                    cutoffmethod = 'depth',    #choices are 'depth' (default), 'height', 'number'
                    cutoff = 3                 #default is '3'
)

###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/b_cells")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(b_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)


# 24 hr all cell types DiABZi
list_of_data <- list(DA25nm.dmso.markers.24.cd4.wt, DA25nm.dmso.markers.24.cd8t.wt, DA25nm.dmso.markers.24.nk.wt, DA25nm.dmso.markers.24.mono.wt, DA25nm.dmso.markers.24.b.wt)

DA25_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

DA25_combined[is.na(DA25_combined)] <- 0

DA25_combined2 <- DA25_combined[,-1]

rownames(DA25_combined2) <- DA25_combined[,1]

head(DA25_combined2)

hmap <- heatmap.L.4(as.matrix(DA25_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "euclidean",
                    clustermethod='ward.D2',
                    cutoffmethod = 'depth',
                    cutoff=0)

###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/24hr_25nm_diabzi")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(DA25_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


# 6 hr all cell types KIN 1148
list_of_data <- list(kin1148.dmso.markers.6.cd4.wt, kin1148.dmso.markers.6.cd8t.wt, kin1148.dmso.markers.6.nk.wt, kin1148.dmso.markers.6.mono.wt, kin1148.dmso.markers.6.b.wt)

kin1148_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

kin1148_combined[is.na(kin1148_combined)] <- 0

kin1148_combined2 <- kin1148_combined[,-1]

rownames(kin1148_combined2) <- kin1148_combined[,1]

head(kin1148_combined2)

hmap <- heatmap.L.4(as.matrix(kin1148_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "euclidean",
                    clustermethod='ward.D2',
                    cutoffmethod = 'depth',
                    cutoff=0)


###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/6hr_kin1148")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(kin1148_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")


# 24 hr all cell types KIN 1148
list_of_data <- list(kin1148.dmso.markers.24.cd4.wt, kin1148.dmso.markers.24.cd8t.wt, kin1148.dmso.markers.24.nk.wt, kin1148.dmso.markers.24.mono.wt, kin1148.dmso.markers.24.b.wt)

kin1148_combined <- list_of_data %>% 
  map(~select(.x, avg_log2FC)) %>% 
  map(rownames_to_column) %>% 
  reduce(full_join, by = "rowname")

kin1148_combined[is.na(kin1148_combined)] <- 0

kin1148_combined2 <- kin1148_combined[,-1]

rownames(kin1148_combined2) <- kin1148_combined[,1]

head(kin1148_combined2)

hmap <- heatmap.L.4(as.matrix(kin1148_combined2),
                    colcolorlist = colcolormatrix,
                    distmethod = "euclidean",
                    clustermethod='ward.D2',
                    cutoffmethod = 'number',
                    cutoff=8)


###WebGestaltR analysis####
setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis/24hr_kin1148")
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism = "hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(kin1148_combined2),
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  temp <- hmap$clustermatrix[row.names(hmap$clustermatrix)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_pvals.csv"))
}

setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02/Cornelius_analysis")










#### Some other stuff I was trying with making heatmaps ####


CD4t.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "CD4+ T", min.pct = 0.25)
mono.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "Monocytes", min.pct = 0.25)
CD8t.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "CD8+ T", min.pct = 0.25)
b.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "B", min.pct = 0.25)
nk.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "NK", min.pct = 0.25)
md.6.markers <- FindMarkers(integration.combined.6hr, ident.1 = "Myeloid dendritic", min.pct = 0.25)

CD4t.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "CD4+ T", min.pct = 0.25)
mono.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "Monocytes", min.pct = 0.25)
CD8t.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "CD8+ T", min.pct = 0.25)
b.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "B", min.pct = 0.25)
nk.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "NK", min.pct = 0.25)
md.24.markers <- FindMarkers(integration.combined.24hr, ident.1 = "Myeloid dendritic", min.pct = 0.25)

cd4.6.genes <- rownames(CD4t.6.markers)

#dittoHeatmap(integration.combined, cd4.6.genes, show_colnames = FALSE, show_rownames = FALSE)
DefaultAssay(integration.combined.6hr.cd4t) <- "RNA"
integration.combined.6hr.cd4t <- subset(integration.combined.6hr, idents = "CD4+ T")

DoHeatmap(integration.combined.6hr)







heatmap.L.4(as.matrix(DA25nm.dmso.markers.6.b))
table(is.nan(DA5nm.dmso.markers[,2]))


avgexp = AverageExpression(integration.combined, return.seurat = T, group.by = "celltype_clusters")
