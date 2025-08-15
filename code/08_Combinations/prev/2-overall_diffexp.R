################################################################
#                                                              #
#                     DEGs - linear fit                        #
#                                                              #
################################################################

### Author: Wangjie
### Date: 12.12.2022
### Datasets: Single cell RNA-seq from TF-seq EXP5-6-7-8-9-10-11-12-13
### Goal: use linear model (edgeRQLF) to extract differential gene expression

library(tidyverse)
library(ggplot2)
library(Seurat)
library(viridis)
library(cowplot)
library(arrow)


code_folder <- ("~/NAS2/TF-seq/Wangjie/TF_resource_paper/code/12-combinations/")
output_folder <- ("~/NAS2/TF-seq/Wangjie/TF_resource_paper/output/12-combinations/")
source(file.path(code_folder, "functions-overall_diffexp.R"))

## DATA ANNOT
data.annot <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Metadata/GRCm38.96_Vector_data.annot.txt", sep = "\t")
rownames(data.annot) <- data.annot$ens_id

seu_all <- readRDS("~/NAS2/TF-seq/Wangjie/TF_resource_paper/output/3-integration/Seurat_integration_functional_TFcells.Rds")
seu_all$cell <- colnames(seu_all)

condition_info <- seu_all@meta.data %>% group_by(TF, batch) %>% tally() %>% filter(TF != "D0") %>% filter(TF != "D0_confluent")
condition_info$condition <- paste0(condition_info$batch, "-", condition_info$TF)

condition_oi <- list(batch = "exp12-13", TF = "Runx2")
condition_oi <- list(batch = "exp10", TF = "Runx2")

batchgroups <- tribble(
  ~batch, ~batchgroup,
  "exp5", "exp5-6",
  "exp6", "exp5-6",
  "exp7", "exp7-8",
  "exp8", "exp7-8",
  "exp9", "exp9",
  "exp10", "exp10-11",
  "exp11", "exp10-11",
  "exp12-13", "exp12-13",
)

filter_seu_all_condition <- function(seu_all, condition_oi) {
  cells_oi_TF <- (seu_all$batch == condition_oi$batch) & (seu_all$TF == condition_oi$TF) & (seu_all$Phase_corrected == "G1")
  batchgroup <- batchgroups %>% filter(batch == condition_oi$batch) %>% pull(batchgroup)
  batches <- batchgroups %>% filter(batchgroup == !!batchgroup) %>% pull(batch)
  cells_oi_D0 <- (seu_all$batch %in% batches) & (seu_all$TF == "D0") & (seu_all$Phase_corrected == "G1")
  
  cells_oi <- cells_oi_TF | cells_oi_D0
  
  if (sum(cells_oi) == 0) {stop("No cells found")}
  if ((sum(cells_oi_TF)) == 0) {stop("No cells found of TF")}
  print(sum(cells_oi))
  seu_diffexp <- seu_all[, cells_oi]
  
  seu_diffexp$TF <- factor(seu_diffexp$TF, levels = c("D0", condition_oi$TF))
  
  seu_diffexp
}

seu_diffexp <- filter_seu_all_condition(seu_all, condition_oi)

scores <- calculate_diffexp_discrete(seu_diffexp)

scores_all <- list()
for (condition in condition_info$condition) {
  if (condition %in% names(scores_all)){
    continue
  }
  print(condition)
  condition_oi <- condition_info %>% filter(condition == !!condition) %>% as.list()
  
  seu_diffexp <- filter_seu_all_condition(seu_all, condition_oi)
  
  scores <- calculate_diffexp_discrete(seu_diffexp)
  
  scores_all[[condition]] <- scores
}

write_rds(scores_all, file.path(output_folder, "scores_all.rds"))

###
scores_all <- read_rds(file.path(output_folder, "scores_all.rds"))

scores_joined <- map2_dfr(scores_all, names(scores_all), function(scores, condition) {
  scores %>% mutate(condition = condition)
})

arrow::write_parquet(scores_joined, file.path(output_folder, "scores.parquet"))
write_csv(condition_info, file.path(output_folder, "conditions_info.csv"))
data.annot %>% rename(gene = ens_id, symbol = gene_short_name) %>% write_csv(file.path(output_folder, "genes_info.csv"))

# condition_oi <- condition_info %>% filter(TF == 'Runx2') %>% last() %>% as.list()
