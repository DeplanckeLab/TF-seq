### Author: Wouter
### Datasets: Single cell RNA-seq from TF-seq EXP12-13
### Goal: Create separate seurat objects for each combination

library(Seurat)
library(tidyverse)
library(edgeR)

source("functions-combinations.R")
source("functions-diffexp.R")

output_folder <- file.path("output")

seu <- read_rds(file.path(output_folder, "seu.rds"))

diffexp_folder <- file.path(output_folder, "diffexp")
if (!dir.exists(diffexp_folder)) {
  dir.create(diffexp_folder, recursive = TRUE)
}

scores_all <- map(combinations, function(combination) {
  print(combination)
  scores <- get_diffexp(combination, seu)
  write_rds(scores, file.path(diffexp_folder, paste0(combination[[1]], "_", combination[[2]], ".rds")))
})