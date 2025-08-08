```R
suppressPackageStartupMessages({
  suppressMessages({
    suppressWarnings({
        library(dplyr)
        library(Seurat)
        library(readr)
        library(ggplot2)
        library("org.Mm.eg.db")
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
    })
  })
})
```


```R
batch_to_label <- c(exp05 = "batch1", exp06 = "batch2", exp07 = "batch3", exp08 = "batch4", exp09 = "batch5", exp10 = "batch6", exp11 = "batch7", `exp12-13`="batch8", exp14 = "batch9")
# scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/4-Seurat_objects/C3H10_10X_all_exps_D0regressed10pc_50pc_integrated_dosealigned.rds ./
# scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/4-Seurat_objects/C3H10_10X_all_exps_merged_genefiltered_50pc_integrated.rds ./
# scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/4-Seurat_objects/Functional_metadata_C3H10_10X_all_exps_merged_genefiltered_50pc_integrated_functional.rds ./Functional_metadata_C3H10_10X_all_exps_merged_genefiltered_50pc_integrated_functional.rds
# scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/4-Seurat_objects/Metadata_dose_aligned.rds ./Metadata_dose_aligned.rds
# scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/4-Seurat_objects/df_allG1Cells_PhaseCorrected_allTFs_D0regressed10pc_50pc_integrated.rds ./df_allG1Cells_PhaseCorrected_allTFs_D0regressed10pc_50pc_integrated.rds
```


```R
plots_folder <- path.expand("plots")
if(!dir.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
}
```


```R
data_folder <- path.expand("../../data")
```


```R
seu <- read_rds(file.path(data_folder, "C3H10_10X_all_exps_D0regressed10pc_50pc_integrated_dosealigned.rds"))
seu2 <- read_rds(file.path(data_folder, "C3H10_10X_all_exps_merged_genefiltered_50pc_integrated.rds"))
```


```R
metadata_dose_aligned <- read_rds(file.path(data_folder, "Metadata_dose_aligned.rds"))
functional_metadata <- read_rds(file.path(data_folder, "Functional_metadata_C3H10_10X_all_exps_merged_genefiltered_50pc_integrated_functional.rds"))
```


```R
seu2$Dose <- seu$Log_Vector_UMI
seu2$Dose_aligned <- metadata_dose_aligned
seu2$functional_cells <- functional_metadata
```


```R
setNames(seu2$functional_cells, colnames(seu2))[colnames(seu)] -> seu$functional_cells
```
