# Power analysis of TF capacity

For not every TF do we detect the same amount of cells. Moreover, the dose levels that each TF reaches differ by several order of magnitudes. This can lead to differences in interpretation, particularly when classifying cells as being high- or low capacity. A TF may be classified as low capacity simply because it does not reach a high-enough dose. To investigate this, we perform a power analysis to determine the minimum dose required to detect a TF as high capacity. We use TFs with high dose levels as reference, and reduce their maximal dose levels to simulate how their capacity classification would change with lower doses.

We use this to ultimately calculate the sensitivity of classifying TFs correctly as low capacity. The cutoff if minimal dose was set at 3.8 to align with a power of approx. 80%.


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
```


```R
data_folder <- path.expand("../../data")
```


```R
seu <- read_rds(file.path(data_folder, "C3H10_10X_all_exps_D0regressed10pc_50pc_integrated_dosealigned.rds"))
```


```R
df <- read_rds(file.path(data_folder, "df_allG1Cells_PhaseCorrected_allTFs_D0regressed10pc_50pc_integrated.rds"))
