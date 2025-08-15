### Author: Wouter
### Date: 28.04.2023
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Calculate differential expression

library(Seurat)
library(tidyverse)
library(edgeR)

setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-diffexp.R")

seu <- read_rds(file.path("output/12-combinations/seu.rds"))

## Score all combination
scores_folder <- file.path("output/12-combinations/scores/discrete")
dir.create(scores_folder, showWarnings = FALSE, recursive = TRUE)
scores_all <- map(combinations, function(combination) {
  print(combination)
  scores <- get_diffexp(combination, seu)
  write_rds(scores, file.path(scores_folder, paste0(combination[[1]], "_", combination[[2]], ".rds")))
})

# write_rds(scores_all, file.path("output/12-combinations/scores/scores_all.rds"))


# scores_all_discrete <- map(combinations, function(combination) {
#   print(combination)
#   get_diffexp_discrete(combination, seu)
# })
# 
# write_rds(scores_all_discrete, file.path("output/12-combinations/scores_all_discrete.rds"))











scores_all <- list()
for (combination in combinations) {
  scores_all[[length(scores_all)+1]] <- read_rds(file.path(scores_folder, paste0(combination[[1]], "_", combination[[2]], ".rds")))
}


plot <- patchwork::wrap_plots(map2(scores_all, combinations, function(scores, combination) {
  patchwork::wrap_plots(
    ggplot(scores) + geom_point(aes(logFC_1, logFC_12, color = significant_12 & significant_1))+ theme_light() + theme(legend.position = "none") + scale_x_continuous(combination[[1]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
    ggplot(scores) + geom_point(aes(logFC_2, logFC_12, color = significant_12 & significant_2)) + theme_light()+ theme(legend.position = "none") + scale_x_continuous(combination[[2]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
    ggplot(scores) + geom_point(aes(logFC_1, logFC_2, color = significant_1 & significant_2)) + theme_light() + theme(legend.position = "none") + scale_x_continuous(combination[[1]]) + scale_y_continuous(combination[[2]])
  )
  
}), ncol = 1)

plot








scores <- read_rds(file.path(scores_folder, paste0(combination[[1]], "_", combination[[2]], ".rds")))

# 
# 
# 
# 
# 
# 
# 
# 
# 
# diffexp$symbol <- diffexp %>% select(gene) %>% left_join(data.annot) %>% pull(gene_short_name)
# diffexp_oi <- diffexp %>% filter(gene %in% c(gene_oi))
# 
# patchwork::wrap_plots(
#   ggplot(diffexp) + geom_point(aes(logFC_1, logFC_12, color = significant_12 & significant_1))+ theme_light() + theme(legend.position = "none") + scale_x_continuous(combination[[1]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
#   ggplot(diffexp) + geom_point(aes(logFC_2, logFC_12, color = significant_12 & significant_2)) + theme_light()+ theme(legend.position = "none") + scale_x_continuous(combination[[2]]) + scale_y_continuous(paste0(combination[[1]], "*", combination[[2]])),
#   ggplot(diffexp, aes(logFC_1, logFC_2)) + 
#     geom_point(aes(color = significant_1 & significant_2)) + 
#     geom_point(data = diffexp_oi, color = "green") + 
#     theme_light() + 
#     theme(legend.position = "none") + 
#     scale_x_continuous(combination[[1]]) + 
#     scale_y_continuous(combination[[2]])
# )
# 
# # look at how the diffexp of one is affected by to the other
# # genes_oi <- diffexp %>% filter(logFC_1 > 0.1) %>% filter(FDR_1 < 0.05) %>% pull(gene)
# genes_oi <- diffexp %>% filter(logFC_2 > 0.1) %>% filter(FDR_2 < 0.05) %>% pull(gene)
# seu_diffexp <- AddModuleScore(seu_diffexp, features = list(oi=genes_oi), name = "oi")
# 
# col = "oi1"
# 
# summarized <- seu_diffexp@meta.data %>% group_by(bin1, bin2) %>% dplyr::summarize(mean = mean(!!ensym(col)), n = n())
# stacked_mean <- reshape2::acast(summarized, list("bin1", "bin2"), value.var = "mean")
# stacked_mean[stacked_n < 5] <- NaN
# pheatmap::pheatmap(
#   stacked_mean,
#   cluster_rows = F,
#   cluster_cols=F,
#   labels_row = combination[[1]],
#   labels_col = combination[[2]],
#   breaks = seq(-max(abs(stacked_mean), na.rm = T), max(abs(stacked_mean), na.rm = T),length.out = 100)
# )
# 
