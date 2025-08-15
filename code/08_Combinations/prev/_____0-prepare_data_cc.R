### Author: Wouter
### Date: 28.04.2023
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Calculate differential expression

library(Seurat)
library(tidyverse)
library(edgeR)

setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-diffexp.R")


##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##

Seurat_int <- readRDS("output/3-integration/Seurat_integration_functional_TFcells.Rds")
Seurat_int <- readRDS("output/3-integration/Seurat_integration_allTFcells.Rds")

# seu <- Seurat_int[, Seurat_int@meta.data$batch %in% "exp12-13"]
seu <- Seurat_int[, (
  ((Seurat_int@meta.data$batch %in% c("exp7", "exp8", "exp9")) & (Seurat_int@meta.data$TF %in% c("D0"))) |
    ((Seurat_int@meta.data$batch %in% c("exp12-13")) & !(Seurat_int@meta.data$TF %in% c("D0")))
)]
# seu <- seu[, seu$Phase == "G1"]
seu@meta.data$cell <- rownames(seu@meta.data)

seu <- seu %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:20)
seu$controls <- ifelse(seu$TF_typeD0specified == "D0", "D0", ifelse(seu$TF_typeD0specified == "D0_confluent", "D0_confluent", "TF"))
plot <- patchwork::wrap_plots(
  DimPlot(seu, group.by = "Phase"),
  DimPlot(seu, group.by = "controls"),
  DimPlot(seu, group.by = "TF", label = TRUE),
  FeaturePlot(seu, "log_vector")
)
plot

data.annot <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Metadata/GRCm38.96_Vector_data.annot.txt", sep = "\t")
rownames(data.annot) <- data.annot$ens_id
data.annot$gene = data.annot$ens_id
seu@assays$RNA@meta.features <- data.annot

## Get TF/barcode counts

xs <- list()

x <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Analysis/Exp12_13/TFEnrichment/Results.Matrix.EXP11_13.FILTERED.txt", check.names = FALSE)
colnames(x) <- paste0(colnames(x), "-12-13")

TFs <- tibble(TF = unique(map_chr(str_split(rownames(x), "_"), 1)))
conditions <- tibble(condition = factor(unique(seu$TF)))

y <- tibble(cell = colnames(x))

for (TF in TFs$TF) {
  y[[TF]] <- colSums(x[startsWith(rownames(x), TF), ])
}

y <- y %>% 
  filter(cell %in% seu@meta.data$cell) %>% 
  complete(cell = rownames(seu@meta.data)) %>% 
  mutate_all(replace_na, 0) %>% 
  mutate(order = match(cell, seu@meta.data$cell)) %>% 
  arrange(order)

oe_assay <- CreateAssayObject(counts = y %>% select(-order) %>% column_to_rownames("cell") %>% as.matrix() %>% t)
seu[["oe"]] <- oe_assay

write_rds(seu, file.path("output/12-combinations/seu_cc.rds"))
