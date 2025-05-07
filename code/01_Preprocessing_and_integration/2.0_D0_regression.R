title: "Regression of the heterogeneity in control cells"
author: "Wouter Saelens / Wangjie Liu"
date: "2024-10-23"



setwd("./")

library(dplyr)
library(Seurat)
library(readr)
library(ggplot2)
library(purrr)
library(Matrix)
library(irlba)


OutputDir <- "./results"

seu <- read_rds("results/C3H10_10X_all_exps_merged_genefiltered_integrated.rds")
seu$Dose <- seu$Log_Vector_UMI
seu$Dose[seu$TF %in% c("D0","D0_confluent")] <- 0

seu$batch_overall <- "batch1"
seu$batch_overall[seu$batch %in% c("exp07","exp08")] <- "batch2"
seu$batch_overall[seu$batch %in% c("exp09")] <- "batch3"
seu$batch_overall[seu$batch %in% c("exp10","exp11")] <- "batch4"
seu$batch_overall[seu$batch %in% c("exp12-13")] <- "batch5"
seu$batch_overall[seu$batch %in% c("exp14")] <- "batch6"

select_control_dataset <- function(seu, batch) {
    filtered <- seu@meta.data |> filter((batch_overall == !!batch) & (TF %in% c("D0", "D0_confluent")) & (Phase_corrected == "G1"))
    # note: Wouter used Phase instead of Phase_corrected
    seu <- subset(seu, cells = rownames(filtered))
    return(seu)
}

select_batch_dataset <- function(seu, batch) {
    filtered <- seu@meta.data |> filter((batch_overall == !!batch))
    seu <- subset(seu, cells = rownames(filtered))
    return(seu)
}

select_tf_dataset <- function(seu, batch, tf) {
    filtered_control <- seu@meta.data |> filter((batch_overall == !!batch) & (TF %in% c("D0", "D0_confluent")))
    filtered_tf <- seu@meta.data |> filter((batch_overall == !!batch) & (TF == !!tf))
    filtered <- rbind(filtered_control, filtered_tf)
    seu <- subset(seu, cells = rownames(filtered))
    return(seu)
}

batches <- c("batch1", "batch2", "batch3", "batch4", "batch5", "batch6")

nPCs_to_regress_D0 <- 10

for (batch in batches) {
    # calculate PCA on control dataset
    seu2 <- select_control_dataset(seu, batch)

    normData <- as.matrix(seu2@assays$RNA$data)
    rownames(normData) <- rownames(seu2)
    featuresUsed <- rownames(seu2) # all genes
    normData <- normData[featuresUsed,]
    normData <- t(normData)
    dim(normData) %>% print()
 
    means <- colMeans(normData) # gene expression means
    scaleData <- (t(normData) - means)
    npcs <- nPCs_to_regress_D0
    pca.results <- prcomp_irlba(t(scaleData), n = npcs)

    cell.embeddings <- predict(pca.results, t(scaleData))

    # Transfer and regress out PCA components on all cells
    seu_batch <- select_batch_dataset(seu, batch)

    # get normData
    normData2 <- as.matrix(seu_batch@assays$RNA$data)
    rownames(normData2) <- rownames(seu_batch)
    normData2 <- normData2[featuresUsed,]
    scaleData2 <- t(normData2 - means)

    cell.embeddings2 <- predict(pca.results, scaleData2)

    correction <- cell.embeddings2 %*% t(pca.results$rotation) 
    # pca.results$rotation:
      ## contains the principal component loadings, also known as the eigenvectors of the covariance matrix of the original data. 
      ## It tells how much each original feature contributes to each PC
      ## these loadings are used to project the original data onto the principal component space. Eg, original data X, the projection is done by X_centered %*% pca.results$rotation
    # pca.results$x:
      ## contains principal component scores, which are the coordinates of the data in the new PC space.
    colnames(correction) <- rownames(seu_batch)
    corrected <- seu_batch@assays$RNA@data - t(correction)

    write_rds(corrected, paste0(OutputDir, "corrected_", batch, "_", nPCs_to_regress_D0,"pc.rds"))
}

# combine
correcteds <- map(c("batch1", "batch2", "batch3", "batch4", "batch5", "batch6"), function(batch) {
    corrected <- read_rds(paste0(OutputDir,"corrected_", batch, "_", nPCs_to_regress_D0,"pc.rds"))
    corrected
})

corrected <- do.call(cbind, correcteds)
corrected <- as.matrix(corrected)
all(colnames(seu) %in% colnames(corrected))
seu@assays$corrected <- CreateAssayObject(data = corrected[, colnames(corrected)], key = "corrected_")
seu[["corrected"]]$counts <- seu[["RNA"]]$counts
seu@assays$RNA <- NULL
seu@assays$integrated <- NULL

## ---------------------------------------- run seurat integration to correct batch on the D0_regressed data
DefaultAssay(seu) <- "corrected" # normalized and D0_regressed
# find HVG on batch
seu.list <- SplitObject(seu, split.by = "batch")
exps <- names(seu.list)
seu.list <- lapply(exps, function(x){
  seu.list[[x]] <- FindVariableFeatures(seu.list[[x]], selection.method = "vst", nfeatures = 2000, verbose = F, assay = "corrected")
})
names(seu.list) <- exps
int_features <- SelectIntegrationFeatures(object.list = seu.list) # default 2000 features

seu.list <- lapply(seu.list, function(x) {
  x <- ScaleData(x, features = int_features, verbose = FALSE)
  x <- RunPCA(x, features = int_features, npcs = 50, verbose = FALSE)
})
names(seu.list) <- exps

options(future.globals.maxSize = 8000 * 1024^2)
npcs.integration <- 50
# Run the default Seurat pipeline on the integrated object
int_anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = int_features, reduction = "rpca", dims = 1:npcs.integration) 
data.seurat_integrated <- IntegrateData(anchorset = int_anchors, verbose = F, new.assay.name = "corrected_integrated")
DefaultAssay(data.seurat_integrated) <- "corrected_integrated"
data.seurat_integrated <- ScaleData(data.seurat_integrated, assay = "corrected_integrated", verbose = F)
data.seurat_integrated <- RunPCA(data.seurat_integrated, assay = "corrected_integrated", npcs = 200, verbose = F)
data.seurat_integrated <- RunTSNE(data.seurat_integrated, reduction = "pca", dims = 1:200)
data.seurat_integrated <- RunUMAP(data.seurat_integrated, reduction = "pca", dims = 1:200)
data.seurat_integrated <- SetIdent(data.seurat_integrated, value = "TF")
print(dim(data.seurat_integrated))

# Save Integrated Seurat object
saveRDS(data.seurat_integrated, "results/C3H10_10X_all_exps_D0regressed_integrated.rds")
