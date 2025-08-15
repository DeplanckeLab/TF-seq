### Author: Wouter Saelens
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Create the Seurat object with only relevant data for the combinations

library(Seurat)
library(tidyverse)
library(edgeR)

# TODO: should be the main directory containing the code, output, etc folders
setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-combinations.R")

### Load Data

# TODO: Should be replaced with loading in the main HDF5 file
Seurat_int <- readRDS("output/3-integration/Seurat_integration_allTFcells.Rds")

# we use the control cells from exp7/8/9
seu <- Seurat_int[, (
  ((Seurat_int@meta.data$batch %in% c("exp7", "exp8", "exp9")) & (Seurat_int@meta.data$TF %in% c("D0"))) |
    ((Seurat_int@meta.data$batch %in% c("exp12-13")) & !(Seurat_int@meta.data$TF %in% c("D0")))
)]

seu@meta.data$cell <- rownames(seu@meta.data)

# qc about the phase and TF
seu <- seu %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:20)
seu$controls <- ifelse(seu$TF_typeD0specified == "D0", "D0", ifelse(seu$TF_typeD0specified == "D0_confluent", "D0_confluent", "TF"))

plot <- patchwork::wrap_plots(
  DimPlot(seu, group.by = "Phase"),
  DimPlot(seu, group.by = "controls"),
  DimPlot(seu, group.by = "TF", label = TRUE),
  FeaturePlot(seu, "log_vector")
)
plot

# load the annotatons of genes and add it to the seurat object
# TODO: Should be replaced
data.annot <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Metadata/GRCm38.96_Vector_data.annot.txt", sep = "\t")
rownames(data.annot) <- data.annot$ens_id
data.annot$gene = data.annot$ens_id
seu@assays$RNA@meta.features <- data.annot

### Get TF/barcode counts
## Get enrichment data and put it in the seurat object
# raw counts
x <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Analysis/Exp12_13/TFEnrichment/Results.Matrix.EXP11_13.FILTERED.txt", check.names = FALSE)
colnames(x) <- paste0(colnames(x), "-12-13")

# calculate both "oe_reads" and "bc_reads"
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
seu[["oereads"]] <- oe_assay

y_bc <- tibble(cell = colnames(x))
for (bc in rownames(x)) {
  y_bc[[bc]] <- as.numeric(x[bc, ])
}
y_bc <- y_bc %>% 
  filter(cell %in% seu@meta.data$cell) %>% 
  complete(cell = rownames(seu@meta.data)) %>% 
  mutate_all(replace_na, 0) %>% 
  mutate(order = match(cell, seu@meta.data$cell)) %>% 
  arrange(order)

bc_assay <- CreateAssayObject(counts = y_bc %>% select(-order) %>% column_to_rownames("cell") %>% as.matrix() %>% t)
seu[["bcreads"]] <- bc_assay

## UMIs
library(data.table)
x1_df <- fread("~/SVRAW1/prainer/TF_scRNAseq_04.2019/TF_enrichment/exp12_13/Results.Matrix.UMI.Exp_12.txt")
x2_df <- fread("~/SVRAW1/prainer/TF_scRNAseq_04.2019/TF_enrichment/exp12_13/Results.Matrix.UMI.Exp_13.txt")

x1 <- as.matrix(x1_df[, 3:ncol(x1_df)])
rownames(x1) <- x1_df$TFId
colnames(x1) <- paste0(colnames(x1), "-1")

x2 <- as.matrix(x2_df[, 3:ncol(x2_df)])
rownames(x2) <- x2_df$TFId
colnames(x2) <- paste0(colnames(x2), "-2")

x <- cbind(x1, x2)
colnames(x) <- paste0(colnames(x), "-12-13")

sum(colnames(seu) %in% colnames(x))

# calculate both "oe_reads" and "bc_reads"
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
seu[["oeumi"]] <- oe_assay

y_bc <- tibble(cell = colnames(x))
for (bc in rownames(x)) {
  y_bc[[bc]] <- as.numeric(x[bc, ])
}
y_bc <- y_bc %>% 
  filter(cell %in% seu@meta.data$cell) %>% 
  complete(cell = rownames(seu@meta.data)) %>% 
  mutate_all(replace_na, 0) %>% 
  mutate(order = match(cell, seu@meta.data$cell)) %>% 
  arrange(order)

bc_assay <- CreateAssayObject(counts = y_bc %>% select(-order) %>% column_to_rownames("cell") %>% as.matrix() %>% t)
seu[["bcumi"]] <- bc_assay

## Save
write_rds(seu, file.path("output/12-combinations/seu.rds"))


