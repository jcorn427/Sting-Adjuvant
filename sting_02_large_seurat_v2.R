setwd("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat")

library(Seurat)
library(BPCells)
library(scCustomize)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(glmGamPoi)
library(tm)
library(svglite)
library(Azimuth)
library(clustifyr)
library(ExperimentHub)
library(clustifyrdata)
library(clustermole)
library(HGNChelper)
library(SingleR)
library(celldex)
library(decontX)
library(SingleCellExperiment)
library(UCell)
library(scater)
library(ChromSCape)
library(dittoSeq)
library(WebGestaltR)
library(singleCellTK)
library(viridis)
library(CellMembrane)

# needs to be set for large dataset analysis
options(future.globals.maxSize = 3e+09)

# specify that you would like to create a Seurat v5 assay
# note that we require setting this option to ensure that existing pipelines are not affected
options(Seurat.object.assay.version = "v5")

# Read in target file
target <- read.csv("../target_file.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# in case of crashing R, run this:
combo_seurat <- readRDS("final_seurat.RDS")

Idents(combo_seurat) <- "Genotype"
combo_seurat_wt <- subset(combo_seurat, idents = "WT")
combo_seurat_haq <- subset(combo_seurat, idents = "HAQ")













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

# Read in the data and initialize the Seurat objects with the filtered

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
      dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat/sample_files/",i)
    )
    
    # Load matrix from directory
    seurat.mat <- open_matrix_dir(dir = paste0("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/Cornelius_analysis/sting_02_large_seurat/sample_files/",i))
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

combo_seurat <- CreateSeuratObject(counts = data_mod)

# Create metadata dataframe to add, and add it to the seurat object

cell_meta <- lapply(cell_meta, function(x){t(x)})

cell_mod <- Replace_Suffix(data = cell_meta, current_suffix = "-1", new_suffix = new_suffixes)

cell_mod <- lapply(cell_mod, function(x){t(x)})

cell_all <- do.call(rbind, cell_mod)

cell_all <- as.data.frame(cell_all)

cell_all$Genotype <- removeNumbers(cell_all$Donor)

combo_seurat <- AddMetaData(object = combo_seurat, metadata = cell_all)

# Add metadata column for samples based on #
meta_data <- rownames(combo_seurat@meta.data)
meta_data <- strsplit(meta_data, split = "-")
meta_data <- as.data.frame(do.call(rbind, meta_data))
combo_seurat <- AddMetaData(combo_seurat, metadata = meta_data[,2], col.name = "Sample")

# Get MT percentages
combo_seurat <- PercentageFeatureSet(combo_seurat, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(combo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Subset based on nFeature_RNA and MT%
combo_seurat <- subset(combo_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Add cell #s here for pre and post filtering

VlnPlot(combo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Remove IL15 data
combo_seurat <- subset(combo_seurat, subset = Treatment == "IL15", invert = TRUE)


# Use SCTransform to normalize - got this to work with integration, this is the preferred method
combo_seurat <- SCTransform(combo_seurat, vst.flavor = "v2", conserve.memory = TRUE, return.only.var.genes = F)

# Run PCA
combo_seurat <- RunPCA(combo_seurat)

# Take a look at PCA
# print(combo_seurat[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(combo_seurat, dims = 1:2, reduction = "pca")
DimPlot(combo_seurat, reduction = "pca", group.by = "Treatment", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Timepoint", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Treatment", split.by = "Treatment", shuffle = TRUE, ncol = 3)
DimPlot(combo_seurat, reduction = "pca", group.by = "Timepoint", split.by = "Timepoint", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Genotype", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Genotype", split.by = "Genotype", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Donor", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Donor", split.by = "Donor", shuffle = TRUE, ncol = 2)
p1 <- DimPlot(combo_seurat, reduction = "pca", group.by = "Sample", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "pca", group.by = "Sample", split.by = "Sample", shuffle = TRUE, ncol = 4) + NoLegend()
ggsave(filename = "full_pca_sample.png", plot = p1)


ElbowPlot(combo_seurat, ndims = 50)

# Run UMAP
combo_seurat <- RunUMAP(combo_seurat, dims = 1:20)

# Find neigbors and clusters
combo_seurat <- FindNeighbors(combo_seurat, dims = 1:20)
combo_seurat <- FindClusters(combo_seurat)

saveRDS(combo_seurat, file = "seurat_full_SCT.RDS")

# Take a look at UMAPs
DimPlot(combo_seurat, group.by = "Sample", shuffle = TRUE)
DimPlot(combo_seurat, shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Treatment", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Timepoint", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Treatment", split.by = "Treatment", shuffle = TRUE, ncol = 3)
DimPlot(combo_seurat, reduction = "umap", group.by = "Timepoint", split.by = "Timepoint", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Genotype", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Genotype", split.by = "Genotype", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Donor", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap", group.by = "Donor", split.by = "Donor", shuffle = TRUE, ncol = 2)
DimPlot(combo_seurat, group.by = "seurat_clusters", shuffle = TRUE)
DimPlot(combo_seurat, group.by = "seurat_clusters", split.by = "seurat_clusters", shuffle = TRUE, ncol = 4) + NoLegend()

# RPCA
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT", assay = "SCT", k.weight = 55
)

saveRDS(combo_seurat, file = "seurat_full_SCT_int.RDS")

combo_seurat <- FindNeighbors(combo_seurat, reduction = "integrated.rpca", dims = 1:20)


# Harmony
combo_seurat <- IntegrateLayers(
  object = combo_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony", normalization.method = "SCT", assay = "SCT", k.weight = 55
)














## Run with resolution at 0.4
combo_seurat <- FindClusters(combo_seurat, resolution = 0.4, cluster.name = "clusters_res_0.4")
combo_seurat <- FindClusters(combo_seurat, resolution = 0.8, cluster.name = "clusters_res_0.8")

combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")
combo_seurat <- RunUMAP(combo_seurat, reduction = "integrated.harmony", dims = 1:20, reduction.name = "umap.harmony")


DimPlot(combo_seurat, reduction = "integrated.rpca", group.by = "Treatment", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "Treatment", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "clusters_res_0.8", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "clusters_res_0.4", shuffle = TRUE, label = TRUE)

DimPlot(combo_seurat, reduction = "umap.harmony", group.by = "Treatment", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.harmony", group.by = "clusters_res_0.8", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.harmony", group.by = "clusters_res_0.4", shuffle = TRUE, label = TRUE)





saveRDS(combo_seurat, file = "seurat_full_SCT_int_umap.RDS")

##### Cell type annotations #####

# Run azimuth
combo_seurat <- RunAzimuth(combo_seurat, reference = "pbmcref", assay = "SCT")

# Check azimuth assignments
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l2", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l1", shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l2", split.by = "predicted.celltype.l2", ncol = 5) + NoLegend()
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "predicted.celltype.l1", split.by = "predicted.celltype.l1", ncol = 3) + NoLegend()

# Run Clustifyr

# build reference data set

pbmc_ref <- average_clusters(mat = pbmc_matrix, metadata = pbmc_meta$classified, if_log = TRUE)

# get normalized matrix from seurat object

norm_mat <- LayerData(combo_seurat, layer = "data", assay = "SCT")


# run it

res <- clustify(input = norm_mat, ref_mat = pbmc_ref, metadata = combo_seurat@meta.data, cluster_col = "rpca_clusters_res_0.8")

# get metadata to add to seurat object, and add it
cl <- cor_to_call(res, cluster_col = "rpca_clusters_res_0.8")

new_metadata <- call_to_metadata(cl, 
                                 metadata = combo_seurat@meta.data, 
                                 cluster_col = "rpca_clusters_res_0.8")

combo_seurat@meta.data <- new_metadata

# Check out the results
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "type", shuffle = TRUE)


### ClusterMole is next

### Get average expression matrix

#avg_exp_mat <- AverageExpression(combo_seurat)
agg_exp_mat <- AggregateExpression(combo_seurat)
agg_exp_mat_RNA <- as.matrix(agg_exp_mat$RNA)
agg_exp_mat_SCT <- as.matrix(agg_exp_mat$SCT)

agg_exp_mat_RNA <- log1p(agg_exp_mat_RNA)
agg_exp_mat_SCT <- log1p(agg_exp_mat_SCT)

enrich_tbl <- clustermole_enrichment(expr_mat = agg_exp_mat_SCT, species = "hs")

write.csv(enrich_tbl, file = "cluster_mole_enrichments.csv")

DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "rpca_clusters_res_0.8", shuffle = TRUE, label = TRUE)

# ScType
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = combo_seurat[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(combo_seurat@meta.data$rpca_clusters_res_0.8), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(combo_seurat@meta.data[combo_seurat@meta.data$rpca_clusters_res_0.8==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(combo_seurat@meta.data$rpca_clusters_res_0.8==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

combo_seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  combo_seurat@meta.data$customclassif[combo_seurat@meta.data$rpca_clusters_res_0.8 == j] = as.character(cl_type$type[1])
}

DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'customclassif')        




# SingleR
# convert seurat object to single cell experiment
combo_seurat[["RNA"]] <- as(combo_seurat[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(combo_seurat)

# get our reference data
ref <- MonacoImmuneData()

# run SingleR
monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.main)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.fine)

# Check cell type assignments using heatmap
# df is too big to make heatmap, trying subsampling
set.seed(42)
monaco.main.sub <- monaco.main[sample(nrow(monaco.main), 100000), ]
plotScoreHeatmap(monaco.main.sub)

plotScoreHeatmap(monaco.main.sub, 
                 annotation_col=as.data.frame(colData(sce)[,"Donor",drop=FALSE]))

plotDeltaDistribution(monaco.main.sub)


monaco.fine.sub <- monaco.fine[sample(nrow(monaco.fine), 100000), ]
plotScoreHeatmap(monaco.fine.sub)

plotScoreHeatmap(monaco.fine.sub, 
                 annotation_col=as.data.frame(colData(sce)[,"Donor",drop=FALSE]))

plotDeltaDistribution(monaco.fine.sub)






# add metadata to seurat object
combo_seurat@meta.data$monaco.main <- monaco.main$pruned.labels
combo_seurat@meta.data$monaco.fine <- monaco.fine$pruned.labels

DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'monaco.main')
DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'monaco.fine')

DimPlot(combo_seurat, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'monaco.main')
DimPlot(combo_seurat, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'monaco.fine')


# Use singleR to look at other databases
# Human primary cell atlas data
hpca.ref <- celldex::HumanPrimaryCellAtlasData()

# DICE
dice.ref <- celldex::DatabaseImmuneCellExpressionData()

hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)

# add metadata to seurat object
combo_seurat@meta.data$hpca.main <- hpca.main$pruned.labels
combo_seurat@meta.data$hpca.fine <- hpca.fine$pruned.labels
combo_seurat@meta.data$dice.main <- dice.main$pruned.labels
combo_seurat@meta.data$dice.fine <- dice.fine$pruned.labels

DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'hpca.main')


plot <- DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'hpca.fine')
AugmentPlot(plot = plot)

DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'dice.main')
DimPlot(combo_seurat, reduction = "umap.rpca", label = TRUE, repel = TRUE, group.by = 'dice.fine')

#saveRDS(combo_seurat, file = "seurat_full_SCT_int_umap_pred.RDS")
combo_seurat <- readRDS("seurat_full_SCT_int_umap_pred.RDS")


DimPlot(combo_seurat, reduction = "umap.rpca",  label = TRUE, shuffle = TRUE)
DimPlot(combo_seurat, reduction = "umap.rpca", group.by = "monaco.main", shuffle = TRUE, split.by = "monaco.main", ncol = 2) + NoLegend()

## Reorder clusters to make more sense ##
new_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
combo_seurat@active.ident <- factor(combo_seurat@active.ident, levels = new_order)

## Fix an error (I left out the 0 on the new_order)
#combo_seurat$rpca_clusters_res_0.8[is.na(combo_seurat$rpca_clusters_res_0.8)] <- 0

new.cluster.ids <- c("T cell", "Naive CD4 T", "Naive CD4 T", "NK", "CD8 T", "Memory CD4 T", "Naive B", "CD4 T", "Follicular Helper T", "T cell", "NK", "Mono", "T cell", "Naive B", "Unknown", "Tregs", "Memory B", "T cell", "Unknown", "DC", "Tregs", "Mono", "Megakaryocyte", "NK", "DC", "Memory B", "Unknown")
names(new.cluster.ids) <- levels(combo_seurat)
combo_seurat <- RenameIdents(combo_seurat, new.cluster.ids)





#### Verify predicted cell type assignments using canonical marker genes ####

tcells <- FeaturePlot(combo_seurat, features = "CD3D", reduction = "umap.rpca")
tcells

cd8_t <- FeaturePlot(combo_seurat, features = c("CD8A", "CD8B"), reduction = "umap.rpca")
cd8_t

cd4_t <- FeaturePlot(combo_seurat, features = c("CD4"), reduction = "umap.rpca")
cd4_t

cd8_effector <- FeaturePlot(combo_seurat, features = c("GZMK", "GZMH", "PRF1", "CCL5"), reduction = "umap.rpca")
cd8_effector

cd8_memory <- FeaturePlot(combo_seurat, features = c("ITGB1"), reduction = "umap.rpca")
cd8_memory

cd8_naive <- FeaturePlot(combo_seurat, features = c("CCR7"), reduction = "umap.rpca")
cd8_naive

cd4_naive <- FeaturePlot(combo_seurat, features = c("IL7R", "CCR7"), reduction = "umap.rpca")
cd4_naive

cd4_memory <- FeaturePlot(combo_seurat, features = c("IL7R", "S100A4"), reduction = "umap.rpca")
cd4_memory

treg <- FeaturePlot(combo_seurat, features = c("FOXP3", "IL2RA"), reduction = "umap.rpca")
treg

MAIT <- FeaturePlot(combo_seurat, features = c("SLC4A10", "TRAV1-2"), reduction = "umap.rpca")
MAIT

gamma_delta <- FeaturePlot(combo_seurat, features = c("TRGV9", "TRDV2"), reduction = "umap.rpca")
gamma_delta

b_cells <- FeaturePlot(combo_seurat, features = c("CD79A", "CD79B"), reduction = "umap.rpca")
b_cells

b_naive <- FeaturePlot(combo_seurat, features = c("CD27"), reduction = "umap.rpca")
b_naive

b_plasma <- FeaturePlot(combo_seurat, features = c("SDC1", "MZB1", "XBP1"), reduction = "umap.rpca")
b_plasma

mono_classic <- FeaturePlot(combo_seurat, features = c("CD14", "LYZ"), reduction = "umap.rpca")
mono_classic

mono_nonclassic <- FeaturePlot(combo_seurat, features = c("FCGR3A", "MS4A7"), reduction = "umap.rpca")
mono_nonclassic

dcs <- FeaturePlot(combo_seurat, features = c("IL3RA", "CLEC4C", "CST3", "NRP1"), reduction = "umap.rpca")
dcs

dc_myeloid <- FeaturePlot(combo_seurat, features = c("FCER1A", "CD1C"), reduction = "umap.rpca")
dc_myeloid

dc_plasmacytoid <- FeaturePlot(combo_seurat, features = c("FCER1A", "LILRA4"), reduction = "umap.rpca")
dc_plasmacytoid

nk_cells <- FeaturePlot(combo_seurat, features = c("KLRB1", "KLRC1", "KLRD1", "GNLY", "NKG7", "NCAM1"), reduction = "umap.rpca")
nk_cells

platelets <- FeaturePlot(combo_seurat, features = c("PPBP"), reduction = "umap.rpca")
platelets

## Work on resolving "unknown" clusters ##

Idents(combo_seurat)
combo_seurat <- PrepSCTFindMarkers(combo_seurat, assay = "SCT", verbose = TRUE)
cluster_18_markers <- FindMarkers(combo_seurat, group.by = "rpca_clusters_res_0.8", ident.1 = 18)

cluster_14_markers <- FindMarkers(combo_seurat, group.by = "rpca_clusters_res_0.8", ident.1 = 14)

cluster_26_markers <- FindMarkers(combo_seurat, group.by = "rpca_clusters_res_0.8", ident.1 = 26)
# Errored out saying there's fewer than 3 cells

table(combo_seurat@meta.data$rpca_clusters_res_0.8)
# sure enough, there's only 2 cells in cluster 26. Wat.

# Prep for meet with Leanne
FeaturePlot(combo_seurat, features = c("CD19"), reduction = "umap.rpca")
FeaturePlot(combo_seurat, features = c("CD3D"), reduction = "umap.rpca")
FeaturePlot(combo_seurat, features = c("CD14"), reduction = "umap.rpca")
FeaturePlot(combo_seurat, features = c("CD3D", "CD3E", "CD3G"), reduction = "umap.rpca")

FeaturePlot(combo_seurat_orig, features = "percent.mt", reduction = "umap.rpca")
FeaturePlot(combo_seurat, features = "percent.mt", reduction = "umap.rpca")
FeaturePlot(combo_seurat_orig, features = "nFeature_RNA", reduction = "umap.rpca")
FeaturePlot(combo_seurat_orig, features = "nCount_RNA", reduction = "umap.rpca")

### SoupX work ### Turns out it won't work because there's no support for cell ranger multi

# filt.matrix <- Read10X("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_A/outs/per_sample_outs/MK434_049/count/sample_filtered_feature_bc_matrix/")
# raw.matrix  <- Read10X("/vol08/ngs/Sting-Adjuvant/Sting-Adjuvant_02_full/sting_02_part_2_set2_A/outs/per_sample_outs/MK434_049/count/sample_raw_feature_bc_matrix/")
# srat  <- CreateSeuratObject(counts = filt.matrix)
# soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

### DecontX work ### Skipping this for now, may revisit it later. Was causing broadway to crash when running SCT again
# combo_seurat <- readRDS("seurat_full_SCT.RDS")
# counts <- GetAssayData(object = combo_seurat, layer = "counts")
# sce <- SingleCellExperiment(list(counts = counts))
# sce <- decontX(sce)
# combo_seurat[["decontXcounts"]] <- CreateAssay5Object(counts = decontXcounts(sce))
# 
# DefaultAssay(object = combo_seurat) <- "decontXcounts"
# 
# combo_seurat <- SCTransform(combo_seurat, assay = "decontXcounts", new.assay.name = "decontSCT")


####Perform DE between groups in wt only####
combo_seurat_wt <- PrepSCTFindMarkers(combo_seurat, assay = "SCT")

dmso.mock.markers <- FindMarkers(combo_seurat_wt, ident.1 = "DMSO", ident.2 = "Mock", group.by = "Treatment")
DA5nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "DiABZi_5nM", ident.2 = "DMSO", group.by = "Treatment")
DA25nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "DiABZi_25nM", ident.2 = "DMSO", group.by = "Treatment")
ST002_25nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST002_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST002_50nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST002_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_25nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST012_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_50nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST012_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_25nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST020_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_50nm.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "ST020_50uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_25uM.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "INI3069_25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50uM.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "INI3069_50uM", ident.2 = "DMSO", group.by = "Treatment")
treat_1148_10uM.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "1148_10uM", ident.2 = "DMSO", group.by = "Treatment")
BRIDGING.dmso.markers <- FindMarkers(combo_seurat_wt, ident.1 = "BRIDGING", ident.2 = "DMSO", group.by = "Treatment")

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
write.csv(BRIDGING.dmso.markers, "BRIDGING.dmso_markers.csv")


##### Run UCell #####
markers <- list()
str <- "CH25H
FZD4
SHISA2
OASL
IGFBP6
IFIT2
CCL3L3
HERC5
LRRC4
RSAD2
PLAAT2
PMAIP1
AFAP1
IFIT3
ENPP2
IFNG
IFIT1
CMPK2
CCL3
RNF152
ZC3HAV1
HES4
IFIH1
IL6
PTGS2
NFE2L3
C2orf66
USP18
SOX8
HRH1
SMPD3
RIGI
CCL19
ISG15
GPR161
GBP5
HERC6
CCL4L2
RTP4
IFI44
RHEBL1
TRIM69
OAS3
TNFSF10
EPSTI1
MYO1A
GRAMD1C
LAMP5
GBP4
OR52K2
"

markers$DA5nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "WNT7A
STING1
BMERB1
DSC1
SNED1
CDHR3
SARDH
ITPRIPL1
RGS14
RAB3A
GPA33
FBXL16
HPCAL4
TMEM204
TSPAN32
IGSF9B
AQP3
SDK2
NMUR1
PLEKHG3
FGFBP2
ZNF365
LIME1
KRT73
NHSL2
COLGALT2
TTC16
CCR4
KLF2
NT5DC2
NRG1
FCRL6
GZMH
ADGRG1
OLFM2
SPON1
TMEM45B
ZNF683
TSPAN18
CCDC65
SEMA4C
KRT72
ADAM23
MPV17L
RGL4
THBD
VSIG1
CACNA1I
THBS1
TAFA1
"
markers$DA5nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "IFNB1
IFNW1
IFNA1
CH25H
IGFBP6
FZD4
TNFSF15
SHISA2
IFNG
OASL
IFIT2
PMAIP1
HERC5
RSAD2
PLAAT2
CCL3L3
AFAP1
GBP6
IFIT3
ENPP2
CXCL11
SEMA3A
IFIH1
IFIT1
RNF152
CMPK2
ZBTB32
CHRNA6
CXCL9
USP18
CCL3
HES4
NFE2L3
NEXN
ISG15
CXCL10
ZC3HAV1
CD38
RIGI
RGS1
TRIM69
TNFSF10
HERC6
GBP5
GBP4
IL6
EPSTI1
SMPD3
OAS3
RTP4
"

markers$DA25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "NRG1
NHSL2
GZMH
PLXDC2
RASAL1
SEMA4C
SEMA6B
COLGALT2
IL9R
CLEC7A
TGFBI
PTGDR
KCNE3
TTC16
OLFML2B
RCN3
CYP27A1
ADAM23
FGFBP2
SDK2
CD14
PLPP3
PLEKHG3
SIGLEC9
HS3ST1
NLRC4
HOMER3
ZNF683
FCRL6
MPV17L
GPR162
MCEMP1
C5AR2
SDC2
MERTK
CD300LB
ADGRG1
EEPD1
NMUR1
PID1
OLR1
KRT72
FCAR
SLCO2B1
KRT73
HTRA1
PDK4
PTGFRN
TAFA1
THBS1
"
markers$DA25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "CSF3
SERPINB7
MUCL1
IL6
LINC02898
MMP1
TNFAIP6
GSTM1
TBC1D3B
CCL3L3
INHBA
SMAD1
PTGES
PTGS2
FBXO39
CEMIP
MT1G
CCL19
C11orf96
C1QTNF1
GREM2
CCL3
MT1M
CXCL13
CA12
DRAXIN
OR7D2
IFNG
ADGRG3
IGHG3
SNED1
SLC27A2
G0S2
GNLY
ZNRF1
ACOD1
PLEKHG2
ASTL
B3GAT1
DNAH17
IL1B
SLC7A5
PARD6G
CEL
ARHGAP33
PCBP3
PDGFRB
PTPN13
B3GNTL1
MT1H
"
markers$KIN1148_up <- as.list(strsplit(str, "\n")[[1]])

str <- "DPYSL4
C3
NRGN
FOS
SAXO2
CYP27A1
CCL24
ST14
KRT73
TGFBI
CD300LB
LCN10
CXCL9
LRRN3
OLR1
MERTK
PID1
UTS2
CD109
COX20
FPR3
CD36
HS3ST1
RAB7B
CD209
DSC1
THBS1
PTGFRN
HTRA1
MMP9
CXCL10
GPNMB
SEMA6B
PPBP
THBD
RNASE1
MRC1
COL23A1
MS4A6A
MDGA1
SLCO2B1
NUPR1
NRG1
CXCL11
STAB1
CCL7
PDK4
CCL2
CLEC5A
CLEC10A
"

markers$KIN1148_down <- as.list(strsplit(str, "\n")[[1]])

str <- "IFNB1
IFNG
FZD4
SHISA2
OASL
IFIT2
CSF3
RSAD2
HERC5
PMAIP1
GBP6
TNFAIP6
IL6
PLAAT2
IFIT3
AFAP1
NFE2L3
RGS1
RNF152
CMPK2
IFIT1
CCL3L3
ENPP2
MT1G
USP18
MT1M
ZC3HAV1
ISG15
NEXN
IFIH1
CHI3L1
EPSTI1
CXCL11
GBP5
RIGI
DDX60
RTP4
GRAMD1C
LAMP3
CD38
OR52K2
CES1
GBP4
RUFY4
PTGS2
HES4
RHEBL1
FAM124A
MX1
CCL3
"
markers$INI3069_25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "ITGAM
SPON1
SERPINB2
ZNF683
VSIG1
ADAMTS10
THBS1
SPRED1
CCR2
LCN10
KLHL14
RGL4
EPHA1
MYLIP
MRC1
CLEC7A
MAP2K6
CAMK1
CACNA1I
MGAT5B
CCL7
SARDH
NT5DC2
THBD
ADAM23
MS4A6A
C5AR2
CD163
TGFBI
TAFA1
OLR1
FXYD2
SDC3
C16orf74
CLEC5A
CLEC10A
ST14
NRG1
HTRA1
MERTK
RASAL1
CCL2
STAB1
CD300LB
PTGFRN
SLCO2B1
EEPD1
HS3ST1
ABCG1
PDK4
"
markers$INI3069_25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "IFNA4
IFNA1
IFNB1
IFNA8
IFNA16
IFNW1
IFNL1
IFNA7
IFNG
TNFSF15
FZD4
RAD21L1
IFIT2
OASL
SHISA2
HERC5
PLEKHA4
RORB
RSAD2
NGFR
NRTN
SEMA3A
PMAIP1
AFAP1
NFE2L3
IFIT3
ENPP2
IL6
PLAAT2
GBP6
IFIH1
IFIT1
RNF152
ZC3HAV1
PTGS2
CCL3L3
CCL3
CMPK2
RGS1
RIGI
DDX60
TNFAIP6
USP18
RHEBL1
RTP4
GRAMD1C
ISG15
EPSTI1
LTA
MAMLD1
"
markers$INI3069_50nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "TLR4
SPRED1
KLHL14
CCR2
MS4A6A
CACNA1I
CLEC10A
NHSL2
GPR68
CD101
RCN3
WNT10A
CFAP96
ST14
GPER1
IL9R
C16orf74
PLXDC2
FCAR
CLEC5A
CXCR3
RNASE6
NLRC4
GPR162
COL23A1
DIPK1B
MAP2K6
TGFBI
HTRA1
OLR1
VPREB3
C5AR2
STAB1
CLEC7A
OLFML2B
EPHA1
RASAL1
FCRL6
ABCG1
THBS1
CCL2
MERTK
EEPD1
CCL7
NRG1
CD300LB
HS3ST1
PTGFRN
SLCO2B1
PDK4
"
markers$INI3069_50nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "RANBP3L
STPG3
PPP1R27
PHYHIP
CIST1
GARIN4
DAND5
GAST
ACTRT3
PERM1
HSPA1B
IDI2
TRIM72
PLEKHA4
PTGER1
TMPRSS9
TRIM10
PRPH
VASN
PLK2
HSPB9
FBXO47
TMEM145
HSPA6
PDE4C
THSD8
DHH
GARIN5B
FAP
TSGA10IP
SMIM6
LY6G6E
C16orf92
GNG14
RNASE10
PPP1R36
HSPA1A
LY6G6C
PLVAP
PGF
GCG
TMEM262
RLN3
TMEM151B
CETP
HIPK4
FOSB
LDHAL6B
RXRG
IFNB1
"
markers$ST002_25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "B3GNTL1
RGS14
MSS51
TMA7B
CCDC180
PNOC
VASH1
DDHD2
IGIP
RGL4
GPI
GGACT
RASAL1
MAP2K6
CIDEB
NCKIPSD
SORBS3
RCAN2
SGSM1
SLC46A1
PSD4
H1-3
ACKR3
MPI
TAS2R30
TMEM63A
FAM162A
GZMH
GAS6
GZMB
S1PR5
KRT72
ASB13
PPFIA4
SDK2
F2R
OLFM2
RIMKLA
KIAA1671
EPHA1
ANKZF1
NOG
IL9R
KLHL14
ADAM23
IGFBP2
TAFA1
HS3ST1
PFKFB4
KRT73
"
markers$ST002_25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "PTGER1
TRIM72
DAND5
PPP1R27
C16orf92
LY6G6E
FBXO47
TMEM145
CIST1
TRIM10
TMPRSS9
GNAT1
LY6G6C
STPG3
HSPB9
TMEM262
NR0B2
CDSN
PHYHIP
GARIN5B
HIPK4
TSGA10IP
SPEM2
PDE4C
CPB2
IDI2
CRLF1
AMPD1
LSMEM2
TMEM210
AHSG
CETP
BSX
PLVAP
LRRC73
GCG
GAST
SMIM6
TARM1
PRPH
NME5
PNCK
PNLDC1
APC2
FAP
CRX
FAM83E
OR2Y1
ACTL6B
TMEM151B
"
markers$ST002_50nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "EHMT1
RALGAPA1
WDR1
LONP2
PABIR1
EP400
BBS2
FANCF
SKIC3
TMEM86B
ZNF770
CSNK2A1
ECHDC2
NPEPL1
KLHL14
MAP10
PRKDC
TAFA1
NSD2
TMEM177
USP9Y
FCRL5
PRKN
MUSTN1
COL18A1
DDX28
MTLN
ZNF133
PNOC
CATSPER2
FCRL1
FAM30A
ZC3H12B
PET117
XAF1
ZNF287
TENM1
IGHA2
TRAP1
UTY
TLR3
MT-CO3
IGIP
TAS2R14
RIC3
KDM5D
CTC1
EBLN2
CXCL11
TAS2R30
"
markers$ST002_50nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "FZD4
SHISA2
OASL
IFIT2
PLAAT2
HERC5
RSAD2
PMAIP1
IFIT3
OR52K2
AFAP1
CMPK2
USP18
RTP4
IFIT1
RHEBL1
DDX60
MAMLD1
ZC3HAV1
RUFY4
RIGI
NEXN
MX1
NFE2L3
PIK3R3
IFI44
HERC6
IFIH1
OAS3
ZBP1
EPSTI1
ISG15
IFI44L
CYP2J2
OAS1
GBP5
SPATS2L
LAMP3
ENPP2
TNFSF10
CD38
SOX8
IFI6
IGFBP3
GRAMD1C
TMEM229B
APOL4
SMAD1
ETV7
FRMD3
"
markers$ST012_25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "ARAP3
PTGDR
PDK4
NCS1
CD209
MMP19
ARMC9
CD93
MRC1
NMUR1
PLAU
HMOX1
PHLDA1
SLAMF9
TM4SF19
FAM20C
CD300LB
SDC2
RASAL1
STAB1
CD14
COL23A1
ANPEP
C3
ADGRG1
TGFBI
ST14
PTGFRN
SEMA6B
NT5DC2
CCL24
ABCG1
HOMER3
C5AR2
PPBP
CXCL3
HTRA1
FCAR
MMP9
MERTK
MCEMP1
PID1
PLPP3
EEPD1
CLEC5A
HS3ST1
THBD
THBS1
SERPINB2
CXCL5
"
markers$ST012_25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "FZD4
OASL
SHISA2
IFIT2
HERC5
RSAD2
OR52K1
IFIT3
PMAIP1
IFIT1
OR52K2
CMPK2
AFAP1
USP18
RUFY4
NEXN
DDX60
RTP4
IFI44
OAS1
NFE2L3
EPSTI1
ZC3HAV1
RIGI
MX1
ISG15
IFNG
IFI44L
TNFSF10
LAMP3
CYP2J2
GBP5
OAS3
IFIH1
HERC6
SPATS2L
ZBP1
ENPP2
GRAMD1C
TNFAIP6
CD38
GBP6
ETV7
RHEBL1
OAS2
SOX8
SAMD9L
EIF2AK2
CXCL11
PLSCR1
"
markers$ST012_50nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "TBXAS1
COL23A1
MMP19
NMUR1
SEMA4C
ADAM23
WNT7A
PPBP
AHRR
FAM20C
ARAP3
VASH1
C16orf74
B3GAT1
SPON1
TGFBI
CYP27A1
C3
CD300LB
STAB1
TSPAN18
C5AR2
SLCO2B1
FXYD2
TAFA1
RASAL1
FCAR
ADGRG1
CD14
ST14
PTGFRN
CCL24
ABCG1
CLEC5A
NT5DC2
MERTK
SEMA6B
PLPP3
HTRA1
PID1
HOMER3
MCEMP1
MMP9
PDK4
EEPD1
THBD
HS3ST1
SERPINB2
CXCL5
THBS1
"
markers$ST012_50nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "OR52K1
OR52K2
CYP2J2
USP18
NEXN
IFIT3
IFIT1
ETV7
MX1
OAS1
LAMP3
OAS3
SPATS2L
CARD17P
IFI6
CMPK2
RSAD2
HERC5
IFIT2
IFI44
SAMD9L
RTP4
TNFSF10
IFI44L
MDK
EIF2AK2
RIGI
PTGS2
PLSCR1
PLK2
KLF5
HELZ2
DDX60
ZBP1
DDX60L
SAMD9
ISG15
IFNG
IFITM1
OAS2
IFIT5
MX2
NMI
PARP9
HERC6
XAF1
PGAP1
NT5C3A
STAT1
LGALS9
"

markers$ST020_25nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "NEK6
SNED1
SARDH
WNT7A
NLRC4
SHTN1
HS3ST1
COLGALT2
NEIL1
ARMC9
C14orf132
KDM5D
SEMA4C
OSBPL10
PID1
GNAO1
NIBAN3
DTX4
STAB1
NRCAM
FADS3
HOMER3
MGLL
OLFML2B
MGAT5B
ZNF704
TNS3
UTS2
CYREN
NMUR1
NCS1
FXYD2
TAFA1
B3GAT1
MMP9
RASAL1
ST14
HMOX1
VASH1
EEPD1
LCN10
DIPK1B
ERBB3
PLPP3
COL19A1
SYT17
KLHL14
C16orf74
CYP27A1
C3
"

markers$ST020_25nm_down <- as.list(strsplit(str, "\n")[[1]])

str <- "OR52K1
OR52K2
IFNG
CYP2J2
USP18
IFIT1
ETV7
NEXN
IFIT3
MX1
OAS1
LAMP3
OAS3
HERC5
SPATS2L
RSAD2
IFI6
CMPK2
CARD17P
IFIT2
EIF2AK2
RTP4
RIGI
TNFSF10
KLF5
SAMD9L
IFI44
PTGS2
HELZ2
PLSCR1
IL1B
DDX60L
DDX60
IFI44L
SAMD9
IFIT5
IFITM1
OAS2
NMI
MX2
STAT1
ZBP1
MDK
HERC6
EPSTI1
PARP9
RGS1
GBP4
PI4K2B
ISG15
"

markers$ST020_50nm_up <- as.list(strsplit(str, "\n")[[1]])

str <- "NRCAM
HOMER3
ADAM23
VASH1
DTX4
RBM47
NT5DC2
ZNF860
CEP295NL
OSBPL10
N4BP3
MGST1
CLEC17A
SH3TC1
CNTLN
WNT7A
MITF
MYOF
EPHA1
NEIL1
GAS6
NSUN7
TNS3
CCL7
PDK4
MERTK
SHTN1
PLD4
NRP2
NIBAN3
HS3ST1
RASAL1
ERBB3
MGLL
STAB1
NCS1
MMP9
MGAT5B
FXYD2
CYP27A1
ST14
EEPD1
SYT17
DIPK1B
C16orf74
LCN10
COL19A1
C3
PLPP3
KLHL14
"
markers$ST020_50nm_down <- as.list(strsplit(str, "\n")[[1]])

combo_seurat <- AddModuleScore_UCell(combo_seurat, features = markers)
signature.names <- paste0(names(markers), "_UCell")

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[7], signature.names[8]), label = F) &
  scale_color_viridis_c()


FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[9], signature.names[10]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[5], signature.names[6]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[3], signature.names[4]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[1], signature.names[2]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[11], signature.names[12]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[13], signature.names[14]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[15], signature.names[16]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[17], signature.names[18]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[19], signature.names[20]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[21], signature.names[22]), label = F) &
  scale_color_viridis_c()

# saveRDS(combo_seurat, file = "seurat_full_SCT_int_harmony_pred_ucell.RDS")
combo_seurat <- readRDS("seurat_full_SCT_int_harmony_pred_ucell.RDS")

#### Post meeting work ####

DimPlot(combo_seurat, reduction = "umap.harmony")
DimPlot(combo_seurat, reduction = "umap.harmony", group.by = "monaco.main", shuffle = TRUE, split.by = "monaco.main", ncol = 3) + NoLegend()

# INI3069 DE work
combo_seurat <- PrepSCTFindMarkers(combo_seurat, assay = "SCT")

INI3069_25uM.dmso.markers <- FindMarkers(combo_seurat, ident.1 = "INI3069_25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50uM.dmso.markers <- FindMarkers(combo_seurat, ident.1 = "INI3069_50uM", ident.2 = "DMSO", group.by = "Treatment")


foldchange = 0.26
pval = 0.05
pctthreshdefault=0.1

INI3069_25uM.dmso.markers[["gene_name"]] <- rownames(INI3069_25uM.dmso.markers)
DE_sig_final <- INI3069_25uM.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
              organism="hsapiens",
              interestGene = genes,
              interestGeneType="genesymbol",
              referenceSet = "genome",
              
              projectName = "INI3069_25uM_ORA",
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



INI3069_50uM.dmso.markers[["gene_name"]] <- rownames(INI3069_50uM.dmso.markers)
DE_sig_final <- INI3069_50uM.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_50uM_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Pseudobulk work

pseudo_combo <- AggregateExpression(combo_seurat, assays = "SCT", return.seurat = T, group.by = c("Treatment", "Timepoint", "Genotype", "monaco.main"))
pseudo_combo <- NormalizeData(pseudo_combo)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)

DimPlot(pseudo_combo, group.by = "monaco.main", reduction = "pca")

## Prep for next meeting ##

DoHeatmap(combo_seurat, features = VariableFeatures(combo_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main") + NoLegend()

plot <- DimPlot(combo_seurat, reduction = "pca") + NoLegend()
LabelClusters(plot = plot, id = "ident")

dittoDimPlot(combo_seurat, "monaco.main", reduction.use = "umap.harmony")
dittoBarPlot(combo_seurat, "monaco.main", group.by = "Treatment")

## Subset, no bridging ##
combo_seurat_sub <- subset(x = combo_seurat, subset = Treatment != "BRIDGING")
dittoBarPlot(combo_seurat_sub, "monaco.main", group.by = "Treatment")







# SingleR
# convert seurat object to single cell experiment
combo_seurat[["RNA"]] <- as(combo_seurat[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(combo_seurat_sub)

# get our reference data
ref <- MonacoImmuneData()

# Subset Single R ref data to remove basophils and neutrophils, then rerunning SingleR ##
pbmc_ref <- ref[, ref$label.main != "Basophils"]
pbmc_ref <- pbmc_ref[, pbmc_ref$label.main != "Neutrophils"]

# run SingleR
monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = pbmc_ref, labels = pbmc_ref$label.main)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = pbmc_ref, labels = pbmc_ref$label.fine)

# Check cell type assignments using heatmap
# df is too big to make heatmap, trying subsampling
set.seed(42)
monaco.main.sub <- monaco.main[sample(nrow(monaco.main), 100000), ]
plotScoreHeatmap(monaco.main.sub)

plotScoreHeatmap(monaco.main.sub, 
                 annotation_col=as.data.frame(colData(sce)[,"Donor",drop=FALSE]))

plotDeltaDistribution(monaco.main.sub)


monaco.fine.sub <- monaco.fine[sample(nrow(monaco.fine), 100000), ]
plotScoreHeatmap(monaco.fine.sub)

plotScoreHeatmap(monaco.fine.sub, 
                 annotation_col=as.data.frame(colData(sce)[,"Donor",drop=FALSE]))

plotDeltaDistribution(monaco.fine.sub)






# add metadata to seurat object
combo_seurat_sub@meta.data$monaco.main <- monaco.main$pruned.labels
combo_seurat_sub@meta.data$monaco.fine <- monaco.fine$pruned.labels


DimPlot(combo_seurat_sub, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'monaco.main')
DimPlot(combo_seurat_sub, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'monaco.fine')

combo_seurat <- combo_seurat_sub

saveRDS(combo_seurat_sub, file = "final_seurat.RDS")
combo_seurat <- readRDS("final_seurat.RDS")

# More work for upcoming meeting

dittoBarPlot(combo_seurat, "monaco.main", group.by = "Timepoint")

INI3069_25uM_seurat <- subset(x = combo_seurat, idents = "INI3069_25uM")
dittoBarPlot(INI3069_25uM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "INI3069 25uM Cell Composition")
dittoBarPlot(INI3069_25uM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "INI3069 25uM Cell Composition", scale = "count")

INI3069_50uM_seurat <- subset(x = combo_seurat, idents = "INI3069_50uM")
dittoBarPlot(INI3069_50uM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "INI3069 50uM Cell Composition")
dittoBarPlot(INI3069_50uM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "INI3069 50uM Cell Composition", scale = "count")

DiABZI_5nM_seurat <- subset(x = combo_seurat, idents = "DiABZi_5nM")
dittoBarPlot(DiABZI_5nM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "DiABZi 5nM Cell Composition")

DiABZI_25nM_seurat <- subset(x = combo_seurat, idents = "DiABZi_25nM")
dittoBarPlot(DiABZI_25nM_seurat, "monaco.main", group.by = "Timepoint", x.reorder = c(2,1), main = "DiABZi 25nM Cell Composition")




DiABZi_5nM.dmso.markers <- FindMarkers(combo_seurat, ident.1 = "DiABZi_5nM", ident.2 = "DMSO", group.by = "Treatment")
DiABZi_25nM.dmso.markers <- FindMarkers(combo_seurat, ident.1 = "DiABZi_25nM", ident.2 = "DMSO", group.by = "Treatment")


DiABZi_5nM.dmso.markers[["gene_name"]] <- rownames(DiABZi_5nM.dmso.markers)
DE_sig_final <- DiABZi_5nM.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_5nM_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


DiABZi_25nM.dmso.markers[["gene_name"]] <- rownames(DiABZi_25nM.dmso.markers)
DE_sig_final <- DiABZi_25nM.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_25nM_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

DoHeatmap(INI3069_25uM_seurat, features = VariableFeatures(INI3069_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F, split.by = "Timepoint")
DoHeatmap(INI3069_50uM_seurat, features = VariableFeatures(INI3069_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DoHeatmap(DiABZI_5nM_seurat, features = VariableFeatures(DiABZI_5nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_25nM_seurat, features = VariableFeatures(DiABZI_25nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)


dittoHeatmap(combo_seurat, genes = VariableFeatures(combo_seurat)[1:50], cells.use = 1:500)


# More UCell work, see line 1769 for more code

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[7], signature.names[8]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[1], signature.names[2]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[3], signature.names[4]), label = F) &
  scale_color_viridis_c()

# Get cell numbers
Idents(DiABZI_5nM_seurat) <- "Timepoint"
table(Idents(DiABZI_5nM_seurat))

Idents(DiABZI_25nM_seurat) <- "Timepoint"
table(Idents(DiABZI_25nM_seurat))

Idents(INI3069_25uM_seurat) <- "Timepoint"
table(Idents(INI3069_25uM_seurat))

Idents(INI3069_50uM_seurat) <- "Timepoint"
table(Idents(INI3069_50uM_seurat))

# Generate more heatmaps
str <- "ZFP36L2
SASH3
PELI1
PIK3CD
CD38
BST2
TBC1D10C
LEF1
TNFAIP3
IKZF3
MIF
ZAP70
ZBTB7A
CD28
ITGA4
RNF168
ATM
MFNG
ITM2A
PRKCB
ZBTB1
IL7R
CD300A
PTPN6
NBN
XBP1
CD79B
PTPRJ
PTPRC
JAK3
SAMSN1
CEBPG
ITFG2
HDAC5
CASP8
IRF2BP2
STAT6
CD79A
BCL3
BCL6
POU2F2
PTK2B
PKN1
LRRC8A
TPD52
RBPJ
NCKAP1L
EP300
PIK3R1
PTPN2
KLF6
POLM
HSPD1
AHR
CD27
CARD11
PLCG2
LAX1
PLCL2
FNIP1
AKAP17A
LGALS1
ID2
CD74
SWAP70
SLA2
PRKDC
BAX
RIF1
RNF8
DNAJB9
PPP2R3C
"
bcell_markers <- strsplit(str, "\n")
bcell_markers <- unlist(bcell_markers)

set.seed(42)
combo_seurat.subsampled <- combo_seurat[, sample(colnames(combo_seurat), size =500, replace=F)]

combo_seurat <- FindMarkers(combo_seurat, features = bcell_markers)

DefaultAssay(combo_seurat)

DoHeatmap(combo_seurat, features = bcell_markers, cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DoHeatmap(INI3069_50uM_seurat, features = rownames(DE_sig_final))
dittoHeatmap(combo_seurat.subsampled, genes = bcell_markers, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE)

NormalizeData(combo_seurat)
ScaleData(combo_seurat)



dittoBarPlot(combo_seurat,  "monaco.main", group.by = "Timepoint", split.by = "Treatment", x.reorder = c(2,1), main = "Cell Composition, All Treatments")
dittoDimHex(combo_seurat, "monaco.main", reduction.use = "umap.harmony", main = "Cell Density UMAP")


dittoHeatmap(INI3069_50uM_seurat, genes = bcell_markers, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE, scaled.to.max = TRUE, heatmap.colors.max.scaled = viridis(100), cluster_rows = FALSE)

dittoHeatmap(INI3069_25uM_seurat, genes = bcell_markers, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE, scaled.to.max = TRUE, heatmap.colors.max.scaled = viridis(100), cluster_rows = FALSE)

dittoHeatmap(DiABZI_5nM_seurat, genes = bcell_markers, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE, scaled.to.max = TRUE, heatmap.colors.max.scaled = viridis(100), cluster_rows = FALSE)

dittoHeatmap(DiABZI_25nM_seurat, genes = bcell_markers, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE, scaled.to.max = TRUE, heatmap.colors.max.scaled = viridis(100), cluster_rows = FALSE)

# Post-meeting work
dittoDimPlot(combo_seurat, "IP10")
DoHeatmap(combo_seurat, features = "CXCL10")
FetchData(combo_seurat, vars = c("ident", "CXCL10"))
RidgePlot(combo_seurat, features = "CXCL10", layer = "counts")
FeaturePlot(DiABZI_25nM_seurat, features = "CXCL10", reduction = "umap.harmony") 
dittoBarPlot(combo_seurat, "monaco.main", group.by = "Timepoint", split.by = "Treatment", x.reorder = c(2,1) )

# Subsetting to make seurat object for May 10th data
combo_seurat_sub <- subset(x = combo_seurat, subset = Treatment != "1140_10uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST002_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST002_50uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST012_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST012_50uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST020_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST020_50uM")

# This is still too big. Subsetting further
set.seed(42)
combo_seurat_sub_sub <- combo_seurat_sub[, sample(colnames(combo_seurat_sub), size =50000, replace=F)]


dittoBarPlot(combo_seurat_sub, "monaco.main", group.by = "Timepoint", split.by = "Treatment", x.reorder = c(2,1), main = "Cell Composition: Counts", scale = "count")

# Requested genes of interest
gene_list <- c("IP10", "TNF-a", "IFN-g", "MIP1a", "MIP1b", "MIP3a", "MIP3b", "IL-1b", "IL-6", "IL-8")
# Used gene cards to get the gene names we have in our dataset
gene_list <- c("CXCL10", "TNF", "IFNG", "CCL3", "CCL4", "CCL20", "CCL19", "IL1B", "IL6", "CXCL8")
gene_list <- unlist(gene_list)

FeaturePlot(combo_seurat_sub, features = "CXCL8", reduction = "umap.harmony")


DoHeatmap(combo_seurat_sub_sub, features = gene_list, slot = "counts")
dittoHeatmap(combo_seurat_sub, genes = gene_list, annot.by = c("monaco.main", "Treatment"), complex = TRUE, use_raster = TRUE, scaled.to.max = TRUE, heatmap.colors.max.scaled = viridis(100), cluster_rows = FALSE)
VlnPlot(combo_seurat_sub, features = gene_list)

VlnPlot(combo_seurat_sub_sub, features = "CXCL9", layer = "data", assay = "SCT", split.by = "monaco.main")


# Monocyte work
# First off, need to subset out just the monocyte data using the May 10th seurat object
# Set annotations as active identities
Idents(object = combo_seurat_sub) <- "monaco.main"
combo_seurat_sub <- subset(x = combo_seurat_sub, idents = "Monocytes")

# Set idents to timepoint
Idents(object = combo_seurat_sub) <- "Timepoint"

# Explore a bit
DimPlot(combo_seurat_sub, reduction = "umap.harmony", split.by = "Timepoint", group.by = "Treatment")

# Lets make some new metadata columns to aid splitting
combo_seurat_sub$Genotype_Timepoint <- paste(combo_seurat_sub$Genotype, combo_seurat_sub$Timepoint, sep = "_")
combo_seurat_sub$Genotype_Timepoint_Treatment <- paste(combo_seurat_sub$Genotype, combo_seurat_sub$Timepoint, combo_seurat_sub$Treatment, sep = "_")

DimPlot(combo_seurat_sub, reduction = "umap.harmony", split.by = "Genotype_Timepoint_Treatment", group.by = "Treatment", ncol = 3)

# Split into wt and haq
Idents(combo_seurat_sub) <- "Genotype"
combo_seurat_sub_wt <- subset(combo_seurat_sub, idents = "WT")
combo_seurat_sub_haq <- subset(combo_seurat_sub, idents = "HAQ")

# Rerun SCT pipeline on just monocyte cells

combo_seurat_sub_sct <- SCTransform(combo_seurat_sub)

combo_seurat_sub_sct <- RunPCA(combo_seurat_sub_sct)
ElbowPlot(combo_seurat_sub_sct, ndims = 50)
combo_seurat_sub_sct <- RunUMAP(combo_seurat_sub_sct, dims = 1:40)
combo_seurat_sub_sct <- FindNeighbors(combo_seurat_sub_sct, dims = 1:40)
combo_seurat_sub_sct <- FindClusters(combo_seurat_sub_sct)
DimPlot(combo_seurat_sub_sct, label = TRUE, group.by = "seurat_clusters")

Idents(combo_seurat_sub_sct) <- "Treatment"

# Generate new heatmaps similar to slide 11 but for monocytes only
DoHeatmap(combo_seurat_sub_sct, features = VariableFeatures(combo_seurat_sub_sct)[1:100], size = 4,
          angle = 90, label = F)
# Need to add MOAR metadata column
combo_seurat_sub_sct$Treatment_Timepoint <- paste(combo_seurat_sub$Treatment, combo_seurat_sub$Timepoint, sep = "_")

Idents(combo_seurat_sub_sct) <- "Treatment_Timepoint"

DoHeatmap(combo_seurat_sub_sct, features = VariableFeatures(combo_seurat_sub_sct)[1:100], size = 4,
          angle = 45, label = T)


# Got official reply, going with wt only for now
DimPlot(combo_seurat_sub_wt, group.by = "seurat_clusters", reduction = "umap.harmony")

Idents(combo_seurat) <- "Genotype"
combo_seurat_wt <- subset(combo_seurat, idents = "WT")
DimPlot(combo_seurat_wt, group.by = "seurat_clusters", reduction = "umap.harmony")

# Now we're making heatmaps like slide 11 *again*
combo_seurat_wt$Cell_Type_Timepoint <- paste(combo_seurat_wt$monaco.main, combo_seurat_wt$Timepoint, sep = "_")
Idents(combo_seurat_wt) <- "Treatment"
INI3069_25uM_seurat <- subset(x = combo_seurat_wt, idents = "INI3069_25uM")
INI3069_50uM_seurat <- subset(x = combo_seurat_wt, idents = "INI3069_50uM")
DiABZI_5nM_seurat <- subset(x = combo_seurat_wt, idents = "DiABZi_5nM")
DiABZI_25nM_seurat <- subset(x = combo_seurat_wt, idents = "DiABZi_25nM")

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
          angle = 90, group.by = "monaco.main", label = F)
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

# Making venn diagrams for monocytes with overlapping genes between DiABZi vs INI3069
Idents(combo_seurat_wt) <- "monaco.main"
combo_seurat_mono <- subset(combo_seurat_wt, idents = "Monocytes")
combo_seurat_mono$Treatment_Timepoint <- paste(combo_seurat_mono$Treatment, combo_seurat_mono$Timepoint, sep = "_")
Idents(combo_seurat_mono) <- "Treatment_Timepoint"

combo_seurat_mono <- PrepSCTFindMarkers(combo_seurat_mono)

INI3069_25uM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_25uM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
INI3069_50uM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_50uM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
DiABZi_5nM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_5nM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
DiABZi_25nM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_25nM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)

write.csv(INI3069_50uM_6hr.markers, file = "INI3069_50uM_6hr_wt_mono.markers.csv")
write.csv(INI3069_25uM_6hr.markers, file = "INI3069_25uM_6hr_wt_mono.markers.csv")
write.csv(DiABZi_5nM_6hr.markers, file = "DiABZi_5nM_6hr_wt_mono.markers.csv")
write.csv(DiABZi_25nM_6hr.markers, file = "DiABZi_25nM_6hr_wt_mono.markers.csv")


INI3069_25uM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_25uM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
INI3069_50uM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_50uM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
DiABZi_5nM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_5nM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
DiABZi_25nM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_25nM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)

write.csv(INI3069_50uM_24hr.markers, file = "INI3069_50uM_24hr_wt_mono.markers.csv")
write.csv(INI3069_25uM_24hr.markers, file = "INI3069_25uM_24hr_wt_mono.markers.csv")
write.csv(DiABZi_5nM_24hr.markers, file = "DiABZi_5nM_24hr_wt_mono.markers.csv")
write.csv(DiABZi_25nM_24hr.markers, file = "DiABZi_25nM_24hr_wt_mono.markers.csv")
# Used venny to generate venn diagrams from lists of genes

# pathway analysis for Monocytes
foldchange = 0.26
pval = 0.05
pctthreshdefault=0.1

INI3069_25uM_6hr.markers[["gene_name"]] <- rownames(INI3069_25uM_6hr.markers)
DE_sig_final <- INI3069_25uM_6hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_25uM_6hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

INI3069_50uM_6hr.markers[["gene_name"]] <- rownames(INI3069_50uM_6hr.markers)
DE_sig_final <- INI3069_50uM_6hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_50uM_6hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

DiABZi_5nM_6hr.markers[["gene_name"]] <- rownames(DiABZi_5nM_6hr.markers)
DE_sig_final <- DiABZi_5nM_6hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_5nM_6hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

DiABZi_25nM_6hr.markers[["gene_name"]] <- rownames(DiABZi_25nM_6hr.markers)
DE_sig_final <- DiABZi_25nM_6hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_25nM_6hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


INI3069_25uM_24hr.markers[["gene_name"]] <- rownames(INI3069_25uM_24hr.markers)
DE_sig_final <- INI3069_25uM_24hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_25uM_24hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

INI3069_50uM_24hr.markers[["gene_name"]] <- rownames(INI3069_50uM_24hr.markers)
DE_sig_final <- INI3069_50uM_24hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_50uM_24hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

DiABZi_5nM_24hr.markers[["gene_name"]] <- rownames(DiABZi_5nM_24hr.markers)
DE_sig_final <- DiABZi_5nM_24hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_5nM_24hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


DiABZi_25nM_24hr.markers[["gene_name"]] <- rownames(DiABZi_25nM_24hr.markers)
DE_sig_final <- DiABZi_25nM_24hr.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DiABZi_25nM_24hr_mono_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Vln plot of CXCL10 for WT and HAQ

VlnPlot(combo_seurat_wt, "CXCL10", pt.size = 0)

Idents(combo_seurat) <- "Genotype"

combo_seurat_haq <- subset(combo_seurat, idents = "HAQ")

Idents(combo_seurat_haq) <- "monaco.main"

VlnPlot(combo_seurat_haq, "CXCL10")

RidgePlot(combo_seurat_wt, "CXCL10")
RidgePlot(combo_seurat_haq, "CXCL10")

# reordering cell types on violin plots
combo_seurat_wt_copy <- combo_seurat_wt
combo_seurat_haq_copy <- combo_seurat_haq

my_levels <- c("T cells", "CD4+ T cells", "NK cells", "Monocytes", "Dendritic cells", "B cells", "Progenitors", "NA")

Idents(combo_seurat_wt_copy) <- factor(Idents(combo_seurat_wt_copy), levels = my_levels)

VlnPlot(combo_seurat_wt_copy, "CXCL10", pt.size = 0)
VlnPlot(combo_seurat_wt_copy, "CXCL10")

Idents(combo_seurat_haq_copy) <- factor(Idents(combo_seurat_haq_copy), levels = my_levels)

VlnPlot(combo_seurat_haq_copy, "CXCL10", pt.size = 0)
VlnPlot(combo_seurat_haq_copy, "CXCL10")

RidgePlot(combo_seurat_wt_copy, "CXCL10")
RidgePlot(combo_seurat_haq_copy, "CXCL10")


# Repeating analysis performed on INI3069 but in ST002
# Subsetting to make seurat object for ST002
combo_seurat_sub <- subset(x = combo_seurat_wt, subset = Treatment != "1148_10uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "INI3069_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "INI3069_50uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST012_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST012_50uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST020_25uM")
combo_seurat_sub <- subset(x = combo_seurat_sub, subset = Treatment != "ST020_50uM")

dittoBarPlot(combo_seurat_sub, "monaco.main", group.by = "Timepoint", split.by = "Treatment", x.reorder = c(2,1), main = "Cell Composition: Counts", scale = "count")
dittoBarPlot(combo_seurat_sub, "monaco.main", group.by = "Timepoint", split.by = "Treatment", x.reorder = c(2,1), main = "Cell Composition")

# Post-crash recovery code
Idents(combo_seurat) <- "Genotype"
combo_seurat_wt <- subset(combo_seurat, idents = "WT")
combo_seurat_haq <- subset(combo_seurat, idents = "HAQ")

# Re-generate ST002 only seurat object

# Now we're making heatmaps like slide 11 *again*
combo_seurat_wt$Cell_Type_Timepoint <- paste(combo_seurat_wt$monaco.main, combo_seurat_wt$Timepoint, sep = "_")
Idents(combo_seurat_wt) <- "Treatment"
ST002_25uM_seurat <- subset(x = combo_seurat_wt, idents = "ST002_25uM")
ST002_50uM_seurat <- subset(x = combo_seurat_wt, idents = "ST002_50uM")
DiABZI_5nM_seurat <- subset(x = combo_seurat_wt, idents = "DiABZi_5nM")
DiABZI_25nM_seurat <- subset(x = combo_seurat_wt, idents = "DiABZi_25nM")

Idents(ST002_25uM_seurat) <- "Timepoint"
Idents(ST002_50uM_seurat) <- "Timepoint"
Idents(DiABZI_25nM_seurat) <- "Timepoint"
Idents(DiABZI_5nM_seurat) <- "Timepoint"


ST002_25uM_seurat_6 <- subset(x = ST002_25uM_seurat, idents = "6hr")
ST002_50uM_seurat_6 <- subset(x = ST002_50uM_seurat, idents = "6hr")
DiABZI_5nM_seurat_6 <- subset(x = DiABZI_5nM_seurat, idents = "6hr")
DiABZI_25nM_seurat_6 <- subset(x = DiABZI_25nM_seurat, idents = "6hr")

ST002_25uM_seurat_24 <- subset(x = ST002_25uM_seurat, idents = "24hr")
ST002_50uM_seurat_24 <- subset(x = ST002_50uM_seurat, idents = "24hr")
DiABZI_5nM_seurat_24 <- subset(x = DiABZI_5nM_seurat, idents = "24hr")
DiABZI_25nM_seurat_24 <- subset(x = DiABZI_25nM_seurat, idents = "24hr")



DoHeatmap(ST002_25uM_seurat_6, features = VariableFeatures(ST002_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST002_50uM_seurat_6, features = VariableFeatures(ST002_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_5nM_seurat_6, features = VariableFeatures(DiABZI_5nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_25nM_seurat_6, features = VariableFeatures(DiABZI_25nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DoHeatmap(ST002_25uM_seurat_24, features = VariableFeatures(ST002_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST002_50uM_seurat_24, features = VariableFeatures(ST002_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_5nM_seurat_24, features = VariableFeatures(DiABZI_5nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(DiABZI_25nM_seurat_24, features = VariableFeatures(DiABZI_25nM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

# Generate more figures for other compounds as well

INI3069_25uM_seurat <- subset(x = combo_seurat_wt, idents = "INI3069_25uM")
INI3069_50uM_seurat <- subset(x = combo_seurat_wt, idents = "INI3069_50uM")
ST012_25uM_seurat <- subset(x = combo_seurat_wt, idents = "ST012_25uM")
ST012_50uM_seurat <- subset(x = combo_seurat_wt, idents = "ST012_50uM")
ST020_25uM_seurat <- subset(x = combo_seurat_wt, idents = "ST020_25uM")
ST020_50uM_seurat <- subset(x = combo_seurat_wt, idents = "ST020_50uM")

Idents(INI3069_25uM_seurat) <- "Timepoint"
Idents(INI3069_50uM_seurat) <- "Timepoint"
Idents(ST012_25uM_seurat) <- "Timepoint"
Idents(ST012_50uM_seurat) <- "Timepoint"
Idents(ST020_25uM_seurat) <- "Timepoint"
Idents(ST020_50uM_seurat) <- "Timepoint"

INI3069_25uM_seurat_6 <- subset(x = INI3069_25uM_seurat, idents = "6hr")
INI3069_50uM_seurat_6 <- subset(x = INI3069_50uM_seurat, idents = "6hr")
ST012_25uM_seurat_6 <- subset(x = ST012_25uM_seurat, idents = "6hr")
ST012_50uM_seurat_6 <- subset(x = ST012_50uM_seurat, idents = "6hr")
ST020_25uM_seurat_6 <- subset(x = ST020_25uM_seurat, idents = "6hr")
ST020_50uM_seurat_6 <- subset(x = ST020_50uM_seurat, idents = "6hr")

INI3069_25uM_seurat_24 <- subset(x = INI3069_25uM_seurat, idents = "24hr")
INI3069_50uM_seurat_24 <- subset(x = INI3069_50uM_seurat, idents = "24hr")
ST012_25uM_seurat_24 <- subset(x = ST012_25uM_seurat, idents = "24hr")
ST012_50uM_seurat_24 <- subset(x = ST012_50uM_seurat, idents = "24hr")
ST020_25uM_seurat_24 <- subset(x = ST020_25uM_seurat, idents = "24hr")
ST020_50uM_seurat_24 <- subset(x = ST020_50uM_seurat, idents = "24hr")

DoHeatmap(INI3069_25uM_seurat_6, features = VariableFeatures(INI3069_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(INI3069_50uM_seurat_6, features = VariableFeatures(INI3069_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST012_25uM_seurat_6, features = VariableFeatures(ST012_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST012_50uM_seurat_6, features = VariableFeatures(ST012_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST020_25uM_seurat_6, features = VariableFeatures(ST020_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST020_50uM_seurat_6, features = VariableFeatures(ST020_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

DoHeatmap(INI3069_25uM_seurat_24, features = VariableFeatures(INI3069_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(INI3069_50uM_seurat_24, features = VariableFeatures(INI3069_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST012_25uM_seurat_24, features = VariableFeatures(ST012_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST012_50uM_seurat_24, features = VariableFeatures(ST012_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST020_25uM_seurat_24, features = VariableFeatures(ST020_25uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
DoHeatmap(ST020_50uM_seurat_24, features = VariableFeatures(ST020_50uM_seurat)[1:100], cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)

# Generate figure using top 15 DE genes from each cell type and generate a heatmap from that

# get top genes per cell type
Idents(combo_seurat_wt) <- "monaco.main"

combo_seurat_wt <- PrepSCTFindMarkers(combo_seurat_wt, assay = "SCT", verbose = TRUE)

T_cells.markers <- FindMarkers(combo_seurat_wt, ident.1 = "T cells", group.by = "monaco.main")
CD4_T_cells.markers <- FindMarkers(combo_seurat_wt, ident.1 = "CD4+ T cells", group.by = "monaco.main")
B_cells.markers <- FindMarkers(combo_seurat_wt, ident.1 = "B cells", group.by = "monaco.main")
CD8_T_cells.markers <- FindMarkers(combo_seurat_wt, ident.1 = "CD8+ T cells", group.by = "monaco.main")
NK_cells.markers <- FindMarkers(combo_seurat_wt, ident.1 = "NK cells", group.by = "monaco.main")
Monocytes.markers <- FindMarkers(combo_seurat_wt, ident.1 = "Monocytes", group.by = "monaco.main")
Dendritic.markers <- FindMarkers(combo_seurat_wt, ident.1 = "Dendritic cells", group.by = "monaco.main")
Progenitors.markers <- FindMarkers(combo_seurat_wt, ident.1 = "Progenitors", group.by = "monaco.main")

markers <- list(T_cells.markers, CD4_T_cells.markers, B_cells.markers, CD8_T_cells.markers, NK_cells.markers,
             Monocytes.markers, Dendritic.markers, Progenitors.markers)

for (i in markers) {
  top_genes <- head(i, n = 15)
  print(top_genes)
}

all_markers <- FindAllMarkers(combo_seurat_wt)

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(combo_seurat_wt, features = top10$gene, cells = 1:500, size = 4,
          angle = 90, group.by = "monaco.main", label = F)
# We'll come back to this in the future, maybe.



#Run UCell
Idents(combo_seurat_wt) <- "Timepoint"

combo_seurat_wt_6hr <- subset(x = combo_seurat_wt, idents = "6hr")
combo_seurat_wt_24hr <- subset(x = combo_seurat_wt, idents = "24hr")

combo_seurat_wt_6hr <- PrepSCTFindMarkers(combo_seurat_wt_6hr, assay = "SCT", verbose = TRUE)
combo_seurat_wt_24hr <- PrepSCTFindMarkers(combo_seurat_wt_24hr, assay = "SCT", verbose = TRUE)

DA5nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "DiABZi_5nM", ident.2 = "DMSO", group.by = "Treatment")
DA25nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "DiABZi_25nM", ident.2 = "DMSO", group.by = "Treatment")
ST002_25nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST002_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST002_50nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST002_50uM", ident.2 = "DMSO", group.by = "Treatment")
DA5nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "DiABZi_5nM", ident.2 = "DMSO", group.by = "Treatment")
DA25nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "DiABZi_25nM", ident.2 = "DMSO", group.by = "Treatment")
ST002_25nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST002_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST002_50nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST002_50uM", ident.2 = "DMSO", group.by = "Treatment")





# Using previous UCell marker genes
# Need to regenerate the "signature.names" list above.

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[11], signature.names[12]), label = F) &
  scale_color_viridis_c()

FeaturePlot(combo_seurat, order = TRUE, reduction = "umap.harmony", features = c(signature.names[13], signature.names[14]), label = F) &
  scale_color_viridis_c()




#ST002 ORA
foldchange = 0.26
pval = 0.05
pctthreshdefault=0.1


DA5nm.6hr.dmso.markers[["gene_name"]] <- rownames(DA5nm.6hr.dmso.markers)
DE_sig_final <- DA5nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DA5nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



DA5nm.24hr.dmso.markers[["gene_name"]] <- rownames(DA5nm.24hr.dmso.markers)
DE_sig_final <- DA5nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DA5nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



DA25nm.6hr.dmso.markers[["gene_name"]] <- rownames(DA25nm.6hr.dmso.markers)
DE_sig_final <- DA25nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DA25nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



DA25nm.24hr.dmso.markers[["gene_name"]] <- rownames(DA25nm.24hr.dmso.markers)
DE_sig_final <- DA25nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "DA25nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



ST002_25nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST002_25nm.6hr.dmso.markers)
DE_sig_final <- ST002_25nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST002_25nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



ST002_25nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST002_25nm.24hr.dmso.markers)
DE_sig_final <- ST002_25nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST002_25nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)




ST002_50nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST002_50nm.6hr.dmso.markers)
DE_sig_final <- ST002_50nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST002_50nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)



ST002_50nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST002_50nm.24hr.dmso.markers)
DE_sig_final <- ST002_50nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST002_50nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Regenerate INI3069 ORA, this time for 6 and 24 hrs separately in wt only

INI3069_25nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "INI3069_25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "INI3069_50uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_25nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "INI3069_25uM", ident.2 = "DMSO", group.by = "Treatment")
INI3069_50nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "INI3069_50uM", ident.2 = "DMSO", group.by = "Treatment")

# Find final list of DE markers

INI3069_25nm.6hr.dmso.markers[["gene_name"]] <- rownames(INI3069_25nm.6hr.dmso.markers)
DE_sig_final <- INI3069_25nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_25nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

INI3069_50nm.6hr.dmso.markers[["gene_name"]] <- rownames(INI3069_50nm.6hr.dmso.markers)
DE_sig_final <- INI3069_50nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_50nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

INI3069_25nm.24hr.dmso.markers[["gene_name"]] <- rownames(INI3069_25nm.24hr.dmso.markers)
DE_sig_final <- INI3069_25nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_25nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

INI3069_50nm.24hr.dmso.markers[["gene_name"]] <- rownames(INI3069_50nm.24hr.dmso.markers)
DE_sig_final <- INI3069_50nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "INI3069_50nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Get markers for ST012

ST012_25nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST012_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_50nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST012_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_25nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST012_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST012_50nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST012_50uM", ident.2 = "DMSO", group.by = "Treatment")

# Get markers for ST020

ST020_25nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST020_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_50nm.6hr.dmso.markers <- FindMarkers(combo_seurat_wt_6hr, ident.1 = "ST020_50uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_25nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST020_25uM", ident.2 = "DMSO", group.by = "Treatment")
ST020_50nm.24hr.dmso.markers <- FindMarkers(combo_seurat_wt_24hr, ident.1 = "ST020_50uM", ident.2 = "DMSO", group.by = "Treatment")

#ORA for ST012

# Find final list of DE markers

ST012_25nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST012_25nm.6hr.dmso.markers)
DE_sig_final <- ST012_25nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST012_25nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


# Find final list of DE markers

ST012_50nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST012_50nm.6hr.dmso.markers)
DE_sig_final <- ST012_50nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST012_50nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

ST012_25nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST012_25nm.24hr.dmso.markers)
DE_sig_final <- ST012_25nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST012_25nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

ST012_50nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST012_50nm.24hr.dmso.markers)
DE_sig_final <- ST012_50nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST012_50nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# ST020 ORA

# Find final list of DE markers

ST020_25nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST020_25nm.6hr.dmso.markers)
DE_sig_final <- ST020_25nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST020_25nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


# Find final list of DE markers

ST020_50nm.6hr.dmso.markers[["gene_name"]] <- rownames(ST020_50nm.6hr.dmso.markers)
DE_sig_final <- ST020_50nm.6hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST020_50nm.6hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

ST020_25nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST020_25nm.24hr.dmso.markers)
DE_sig_final <- ST020_25nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST020_25nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)

# Find final list of DE markers

ST020_50nm.24hr.dmso.markers[["gene_name"]] <- rownames(ST020_50nm.24hr.dmso.markers)
DE_sig_final <- ST020_50nm.24hr.dmso.markers %>%
  filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
DE_sig_final <- DE_sig_final %>%
  filter(p_val_adj <= pval) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
DE_sig_final <- DE_sig_final %>%
  filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
  dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

## Web Gestalt ##

###WebGestaltR analysis####

genes <- rownames(DE_sig_final)
print(genes)
WebGestaltR(enrichMethod="ORA",
            organism="hsapiens",
            interestGene = genes,
            interestGeneType="genesymbol",
            referenceSet = "genome",
            
            projectName = "ST020_50nm.24hr.dmso_ORA",
            minNum=5,
            enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                               "geneontology_Molecular_Function_noRedundant",
                               "geneontology_Cellular_Component_noRedundant",
                               "pathway_KEGG"
            ),
            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)


# Monocyte population analysis
# DiABZi and INI3069 first
Idents(combo_seurat_wt) <- "monaco.main"
combo_seurat_mono <- subset(combo_seurat_wt, idents = "Monocytes")
combo_seurat_mono$Treatment_Timepoint <- paste(combo_seurat_mono$Treatment, combo_seurat_mono$Timepoint, sep = "_")
Idents(combo_seurat_mono) <- "Treatment_Timepoint"

combo_seurat_mono <- PrepSCTFindMarkers(combo_seurat_mono)

INI3069_25uM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_25uM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
INI3069_50uM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_50uM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
DiABZi_5nM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_5nM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)
DiABZi_25nM_6hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_25nM_6hr", ident.2 = "DMSO_6hr", recorrect_umi = FALSE)

INI3069_25uM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_25uM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
INI3069_50uM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "INI3069_50uM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
DiABZi_5nM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_5nM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)
DiABZi_25nM_24hr.markers <- FindMarkers(combo_seurat_mono, ident.1 = "DiABZi_25nM_24hr", ident.2 = "DMSO_24hr", recorrect_umi = FALSE)

write.csv(INI3069_50uM_6hr.markers, file = "INI3069_50uM_6hr_wt_mono.markers.csv")
write.csv(INI3069_25uM_6hr.markers, file = "INI3069_25uM_6hr_wt_mono.markers.csv")
write.csv(DiABZi_5nM_6hr.markers, file = "DiABZi_5nM_6hr_wt_mono.markers.csv")
write.csv(DiABZi_25nM_6hr.markers, file = "DiABZi_25nM_6hr_wt_mono.markers.csv")

write.csv(INI3069_50uM_24hr.markers, file = "INI3069_50uM_24hr_wt_mono.markers.csv")
write.csv(INI3069_25uM_24hr.markers, file = "INI3069_25uM_24hr_wt_mono.markers.csv")
write.csv(DiABZi_5nM_24hr.markers, file = "DiABZi_5nM_24hr_wt_mono.markers.csv")
write.csv(DiABZi_25nM_24hr.markers, file = "DiABZi_25nM_24hr_wt_mono.markers.csv")

# Run SCTransform on monocytes
combo_seurat_mono <- SCTransform(combo_seurat_mono, vst.flavor = "v2", conserve.memory = TRUE, return.only.var.genes = F)

# Run PCA
combo_seurat_mono <- RunPCA(combo_seurat_mono)

ElbowPlot(combo_seurat_mono, ndims = 50)

# Run UMAP
combo_seurat_mono <- RunUMAP(combo_seurat_mono, dims = 1:20)

# Find neigbors and clusters
combo_seurat_mono <- FindNeighbors(combo_seurat_mono, dims = 1:20)
combo_seurat_mono <- FindClusters(combo_seurat_mono)

# Harmony
combo_seurat_mono <- IntegrateLayers(
  object = combo_seurat_mono, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony", normalization.method = "SCT", assay = "SCT", k.weight = 55
)


DimPlot(combo_seurat, reduction = "umap.harmony", group.by = "seurat_clusters", shuffle = TRUE)


# Running cellmembrane to help clean up cell calls

combo_seurat <- CellMembrane::RunSingleR(seuratObj = combo_seurat, datasets = c("hpca", "blueprint", "dice", "monaco", "immgen"), assay = "SCT", nThreads = 50 )


pheatmap::pheatmap(chisq.test(table(combo_seurat$SCT_snn_res.0.8, combo_seurat$immgen.label.fine))$residuals, scale = "row")

DimPlot(combo_seurat, reduction = "umap.harmony", split.by = "SCT_snn_res.0.8", shuffle = TRUE, ncol = 3)

Idents(combo_seurat) <- "monaco.main"
combo_seurat_sub <- subset(combo_seurat, idents = "28")
DimPlot(combo_seurat_sub, reduction = "umap.harmony", split.by = "SCT_snn_res.0.8", shuffle = TRUE)

combo_seurat_sub <- subset(combo_seurat, idents = "6")
DimPlot(combo_seurat, reduction = "umap.harmony", split.by = "clusters_res_0.4", shuffle = TRUE, ncol = 3)
