### Author: Wouter
### Date: 28.04.2023
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Calculate differential expression

library(Seurat)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(viridis)

setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-diffexp.R")

# scores_all <- read_rds(file.path("output/12-combinations/scores_all.rds"))

seu <- read_rds(file.path("output/12-combinations/seu.rds"))
seu <- read_rds(file.path("output/12-combinations/seu_cc.rds"))
data.annot <- seu@assays$RNA@meta.features

seu <- seu %>% RunPCA() %>% RunUMAP(dims = 1:20)

## Check out single combination
combination_ix <- 6# ??
# combination_ix <- 2 # ??
# combination_ix <- 7 # PPARG-MYCN
# combination_ix <- 1 # RUNX2-PPARG
# combination_ix <- 3 # CEBPA-MYCN
combination <- combinations[[combination_ix]]

diffexp <- read_rds(file.path("output/12-combinations/scores/discrete/", paste0(combination[[1]], "_", combination[[2]], ".rds")))
# diffexp <- scores_all[[combination_ix]]
print(combination)

TFsoi <- combination
seu_diffexp <- seu[, seu@meta.data %>% filter(TF %in% c(TFsoi[[1]], TFsoi[[2]], "D0", paste0(TFsoi[[1]], "-", TFsoi[[2]]), paste0(TFsoi[[2]], "-", TFsoi[[1]]))) %>% pull(cell)]

seu_diffexp$vector1 <- seu@assays[["oe"]]@counts[TFsoi[[1]],]
seu_diffexp$vector2 <- seu@assays[["oe"]]@counts[TFsoi[[2]],]

seu_diffexp$logvector1 <- log1p(seu_diffexp$vector1)
seu_diffexp$logvector2 <- log1p(seu_diffexp$vector2)

adipo.markers <-  c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36",
                    "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb",
                    "Cidec")
adipo.markers <-  data.annot[data.annot$gene_short_name %in% adipo.markers ,"ens_id"]
DefaultAssay(seu_diffexp) <- "RNA"
seu_diffexp <- AddModuleScore(seu_diffexp, features = list(adipo = adipo.markers), name = "adiposcore")

genes_oi <- diffexp %>% filter(logFC_2 > 0.1) %>% filter(FDR_2 < 0.05) %>% filter((logFC_1 < 0.1) & (logFC_1 > -0.1)) %>% pull(gene)
seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = paste0(combination[[2]], "_unique"))
genes_oi <- diffexp %>% filter(logFC_1 > 0.1) %>% filter(FDR_1 < 0.05) %>% filter((logFC_2 < 0.1) & (logFC_2 > -0.1)) %>% pull(gene)
seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = paste0(combination[[1]], "_unique"))
genes_oi <- diffexp %>% filter(logFC_2 > 0.1) %>% filter(FDR_2 < 0.05) %>% pull(gene)
seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = combination[[1]])
genes_oi <- diffexp %>% filter(logFC_1 > 0.1) %>% filter(FDR_1 < 0.05) %>% pull(gene)
seu_diffexp <- AddModuleScore(seu_diffexp, features = list(a = genes_oi), name = combination[[2]])

symbol_oi <- 'Myog'
symbol_oi <- 'Fabp4'
# symbol_oi <- 'Bglap2'
symbol_oi <- 'Ccnd1'
symbol_oi <- 'Fabp4'
gene_oi <- rownames(data.annot %>% filter(gene_short_name == !!symbol_oi))
seu_diffexp@meta.data[[paste0(symbol_oi)]] <- seu_diffexp@assays$RNA@data[gene_oi, ]
seu_diffexp@meta.data$lognCount_RNA <- log(seu_diffexp$nCount_RNA)

col <- "n"
col <- symbol_oi
# col <- paste0(combination[[1]], "1")
col <- paste0(combination[[2]], "1")
col <- paste0(combination[[1]], "_unique1")
col <- paste0(combination[[2]], "_unique1")
# col <- "adiposcore1"
# col <- "lognCount_RNA"

plot <- ggplot(seu_diffexp@meta.data) + 
  geom_point(aes(logvector1, logvector2, color = !!ensym(col)), size =1) +
  scale_x_continuous(name = combination[[1]]) +
  scale_y_continuous(name = combination[[2]]) +
  scale_color_viridis(name = col, option = "magma", limits = c(0, NA), na.value = "black")+
  coord_fixed() + 
  theme_light() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot

bins1 <- c(0, seq(0.001, max(seu_diffexp@meta.data$logvector1)+0.1, length.out = 5))
bins2 <- c(0, seq(0.001, max(seu_diffexp@meta.data$logvector2)+0.1, length.out = 5))

seu_diffexp@meta.data$bin1 <- cut(seu_diffexp$logvector1, bins1, right = FALSE)
seu_diffexp@meta.data$bin2 <- cut(seu_diffexp$logvector2, bins2, right = FALSE)

summarized <- seu_diffexp@meta.data %>% group_by(bin1, bin2) %>% dplyr::summarize(score = mean(!!ensym(col)), n = n())
n_cells_cutoff <- 1
summarized <- summarized %>% mutate(score = ifelse(n < n_cells_cutoff, NaN, score))

theme_heatmap <- function(plot, combination){
  plot +
  scale_x_discrete(expand = c(0, 0), name = combination[[1]],) +
    scale_y_discrete(expand = c(0, 0), name = combination[[2]]) +
    coord_fixed() + 
    theme_light() +
    # theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

plot <- ggplot(summarized, aes(bin1, bin2, fill= score)) + 
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_viridis(na.value = "#33333333")
plot %>% theme_heatmap(combination)


seu_diffexp@meta.data$bin1 <- seu_diffexp$logvector1 > 0
seu_diffexp@meta.data$bin2 <- seu_diffexp$logvector2 > 0

seu_diffexp_g1 <- seu_diffexp[, seu_diffexp$Phase_corrected == "G1"]
seu_diffexp_cc <- seu_diffexp[, seu_diffexp$Phase_corrected != "G1"]

summarized_g1 <- seu_diffexp_g1@meta.data %>% group_by(bin1, bin2) %>% dplyr::summarize(score = mean(!!ensym(col)), n = n())
summarized_cc <- seu_diffexp_cc@meta.data %>% group_by(bin1, bin2) %>% dplyr::summarize(score = mean(!!ensym(col)), n = n())

shared_limits <- c(min(summarized_g1$score, summarized_cc$score), max(summarized_g1$score, summarized_cc$score))

plot_g1 <- ggplot(summarized_g1, aes(bin1, bin2, fill= score)) + 
  geom_tile() +
  scale_x_discrete(expand = c(0, 0), name = combination[[1]]) +
  scale_y_discrete(expand = c(0, 0), name = combination[[2]]) +
  scale_fill_viridis(na.value = "#33333333", limits = shared_limits, name = col) +
  coord_fixed() + 
  theme_light() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("G1")
plot_g1


plot_cc <- ggplot(summarized_cc, aes(bin1, bin2, fill= score)) + 
  geom_tile() +
  scale_x_discrete(expand = c(0, 0), name = combination[[1]]) +
  scale_y_discrete(expand = c(0, 0), name = combination[[2]]) +
  scale_fill_viridis(na.value = "#33333333", limits = shared_limits, name = col) +
  coord_fixed() + 
  theme_light() +
  # theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cycling")
plot_cc


patchwork::wrap_plots(
  plot_g1, plot_cc
)


##
diffexp_oi <- diffexp %>% filter(gene %in% c(gene_oi))

patchwork::wrap_plots(
  ggplot(diffexp) + geom_point(aes(logFC_1, logFC_12, color = significant_12 & significant_1))+ theme_light() + theme(legend.position = "none") + scale_x_continuous(combination[[1]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
  ggplot(diffexp) + geom_point(aes(logFC_2, logFC_12, color = significant_12 & significant_2)) + theme_light()+ theme(legend.position = "none") + scale_x_continuous(combination[[2]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
  ggplot(diffexp, aes(logFC_1, logFC_2)) +
    geom_point(aes(color = significant_1 & significant_2)) +
    geom_point(data = diffexp_oi, color = "green") +
    theme_light() +
    theme(legend.position = "none") +
    scale_x_continuous(combination[[1]]) +
    scale_y_continuous(combination[[2]])
)
