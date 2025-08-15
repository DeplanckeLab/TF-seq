### Author: Wouter
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Create separate seurat objects for each combination

library(Seurat)
library(tidyverse)
library(edgeR)

setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-combinations.R")

seu <- read_rds(file.path("output/12-combinations/seu.rds"))

## Prepare differential objects
prepare_data <- function(combination, seu) {
  TFsoi <- combination
  seu_diffexp <- seu[, seu@meta.data %>% filter(TF %in% c(TFsoi[[1]], TFsoi[[2]], "D0", paste0(TFsoi[[1]], "-", TFsoi[[2]]), paste0(TFsoi[[2]], "-", TFsoi[[1]]))) %>% pull(cell)]
  
  seu_diffexp$vector1 <- seu@assays[["oeumi"]]@counts[TFsoi[[1]],]
  seu_diffexp$vector2 <- seu@assays[["oeumi"]]@counts[TFsoi[[2]],]
  
  seu_diffexp$logvector1 <- log1p(seu_diffexp$vector1)
  seu_diffexp$logvector2 <- log1p(seu_diffexp$vector2)
  
  seu_diffexp$readvector1 <- seu@assays[["oereads"]]@counts[TFsoi[[1]],]
  seu_diffexp$readvector2 <- seu@assays[["oereads"]]@counts[TFsoi[[2]],]
  
  seu_diffexp$logreadvector1 <- log1p(seu_diffexp$readvector1)
  seu_diffexp$logreadvector2 <- log1p(seu_diffexp$readvector2)
  
  # adipo.markers <-  c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36",
  #                     "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb",
  #                     "Cidec", "Cidea")
  # adipo.markers <-  data.annot[data.annot$gene_short_name %in% adipo.markers ,"ens_id"]
  # DefaultAssay(seu_diffexp) <- "RNA"
  # seu_diffexp <- AddModuleScore(seu_diffexp, features = list(adipo = adipo.markers), name = "adiposcore")
  
  diffexp <- read_rds(file.path(diffexp_folder, paste0(combination[[1]], "_", combination[[2]], ".rds")))
  
  genes_oi <- diffexp %>% filter(logFC_2 > 0.1) %>% filter(FDR_2 < 0.05) %>% filter((logFC_1 < 0.1) & (logFC_1 > -0.1)) %>% pull(gene)
  seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = paste0(combination[[1]], "_unique"))
  genes_oi <- diffexp %>% filter(logFC_2 > 0.1) %>% filter(FDR_2 < 0.05) %>% pull(gene)
  seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = combination[[1]])
  
  genes_oi <- diffexp %>% filter(logFC_1 > 0.1) %>% filter(FDR_1 < 0.05) %>% filter((logFC_2 < 0.1) & (logFC_2 > -0.1)) %>% pull(gene)
  seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = paste0(combination[[2]], "_unique"))
  genes_oi <- diffexp %>% filter(logFC_1 > 0.1) %>% filter(FDR_1 < 0.05) %>% pull(gene)
  seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = combination[[2]])
  
  bins1 <- c(0, seq(0.001, max(seu_diffexp@meta.data$logvector1)+0.1, length.out = 5))
  bins2 <- c(0, seq(0.001, max(seu_diffexp@meta.data$logvector2)+0.1, length.out = 5))
  
  seu_diffexp@meta.data$bin1 <- cut(seu_diffexp$logvector1, bins1, right = FALSE)
  seu_diffexp@meta.data$bin2 <- cut(seu_diffexp$logvector2, bins2, right = FALSE)
  
  bin_info1 <- data.frame(
    name = levels(seu_diffexp$bin1),
    label = c(0, round(bins1[2:(length(bins1)-1)] + (bins1[3:length(bins1)] - bins1[2:(length(bins1)-1)])/2, 1))
  )
  bin_info2 <- data.frame(
    name = levels(seu_diffexp$bin2),
    label = c(0, round(bins2[2:(length(bins2)-1)] + (bins2[3:length(bins2)] - bins2[2:(length(bins2)-1)])/2, 1))
  )
  seu_diffexp@misc$bin_info1 <- bin_info1
  seu_diffexp@misc$bin_info2 <- bin_info2
  
  # add conditions
  conditions_oi <- c("D0", combination[[1]], combination[[2]], paste0(combination[[1]], "-", combination[[2]]))
  
  seu_diffexp@meta.data$condition <- "D0"
  cells_oi <- (seu_diffexp$readvector1 >= 5)
  seu_diffexp@meta.data[cells_oi, "condition"] <- conditions_oi[[2]]
  cells_oi <- (seu_diffexp$readvector2 >= 5)
  seu_diffexp@meta.data[cells_oi, "condition"] <- conditions_oi[[3]]
  cells_oi <- (seu_diffexp$readvector1 >= 5) & (seu_diffexp$readvector2 >= 5)
  seu_diffexp@meta.data[cells_oi, "condition"] <- conditions_oi[[4]]
  
  seu_diffexp@meta.data$condition <- factor(seu_diffexp@meta.data$condition, levels = conditions_oi)
  
  seu_diffexp
}

bc_tf_resetters <- list(
  c("Mycn", "Myog"),
  c("Mycn-1", "Myog"),
  c("Mycn-2", "Myog"),
  c("Mycn-4", "Myog"),
  
  c("Cebpa-1", "Myog"),
  c("Cebpa-4", "Myog"),
  
  c("Pparg-G6", "Runx2"),
  
  c("Cebpa-A6", "Pparg"),
  c("Cebpa", "Pparg"),
  
  c("Mycn-1", "Runx2"),
  c("Mycn-3", "Runx2"),
  c("Mycn-4", "Runx2")
)

for (bc_tf_resetter in bc_tf_resetters) {
  bc <- bc_tf_resetter[[1]]
  tf <- bc_tf_resetter[[2]]
  cells_oi <- rownames(seu@assays$bcreads)[apply(seu@assays$bcreads@counts, 2, which.max)] == bc
  print(sum(cells_oi))
  if (sum(cells_oi) > 0) {
    seu@assays$oereads@counts[tf, cells_oi] <- 0.
  }
  cells_oi <- rownames(seu@assays$bcumi)[apply(seu@assays$bcumi@counts, 2, which.max)] == bc
  print(sum(cells_oi))
  if (sum(cells_oi) > 0) {
    seu@assays$oeumi@counts[tf, cells_oi] <- 0.
  }
}

seu_diffexps <- map(combinations, prepare_data, seu = seu)
names(seu_diffexps) <- combinations %>% map_chr(paste0, collapse = "-")

write_rds(seu_diffexps, "output/12-combinations/seu_diffexps.rds")


### Create metadata
metadatas <- map2_dfr(
  names(seu_diffexps),
  seu_diffexps,
  function(combination_name, seu_diffexp) {seu_diffexp@meta.data %>%  as.data.frame() %>% mutate(combination = combination_name, condition = as.character(condition))
})
metadata <- metadatas %>% group_by(cell) %>% filter(row_number() == 1) %>% filter(condition != "D0")
metadata <- metadata %>% select(cell, condition, vector1, vector2, logvector1, logvector2, adiposcore1, combination, TF)

table(metadata$TF, metadata$condition)



write.table(metadata, file.path(output_folder, "metadata.tsv"), sep = "\t")

metadata2 <- metadatas %>% group_by(cell) %>% filter(row_number() == 1)
table(metadata2$TF, metadata2$condition)
