```R
suppressPackageStartupMessages({
  suppressMessages({
    suppressWarnings({
        library(dplyr, quietly = TRUE)
        library(Seurat, quietly = TRUE)
        library(readr, quietly = TRUE)
        library(ggplot2, quietly = TRUE)
        library("org.Mm.eg.db", quietly = TRUE)

        library(Seurat)
        library(ggplot2)
        library(cowplot)
        library(tidyverse)
        library(dplyr)
        # convert symbols to ensembl ids
        convert_to_ensembl <- function(symbols) {
            ensembl_ids <- mapIds(org.Mm.eg.db, keys = symbols, keytype = "SYMBOL", column = "ENSEMBL")
            return(ensembl_ids)
        }
        convert_to_symbol <- function(ensembl_ids) {
            symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL")
            return(symbols)
        }
        library(purrr)

        plots_folder <- path.expand("plots")
        if (!dir.exists(plots_folder)) {
            dir.create(plots_folder, recursive = TRUE)
        }
    })
  })
})
