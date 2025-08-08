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


```R
# !scp updepla_tf_seq_rwd@fts.epfl.ch:sv_fts-updepla/tf_seq/Wangjie/6-tmp/TF_categories_on_potency_capacity_dosealigned.csv ./
# 
functionality <- read_csv(file.path(data_folder, "TF_categories_on_potency_capacity_dosealigned.csv"))
colnames(functionality)[1] <- "TF_row"

# strip whitespace
functionality$TF <- map_chr(functionality$TF, stringr::str_trim)
functionality$TF_row <- map_chr(functionality$TF_row, stringr::str_trim)
```

    [1m[22mNew names:
    [36m•[39m `` -> `...1`
    [1mRows: [22m[34m232[39m [1mColumns: [22m[34m11[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m (3): ...1, TF, category
    [32mdbl[39m (8): Asym, xmid, scal, Dose.max, Dose.min, startpoint, endpoint, predict...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.


## Scaling


```R
seu2$batch_overall <- unname(c(exp05 = "batch1", exp06 = "batch1", exp07 = "batch2", exp08 = "batch2", exp09 = "batch3", exp10 = "batch4", exp11 = "batch4", `exp12-13`="batch5", exp14 = "batch6")[seu2$batch])
```


```R
max_doses <- seu2@meta.data |> arrange(-Dose) |> group_by(batch_overall, TF) |> summarize(Dose = quantile(Dose, 0.9)) |> ungroup() |> dplyr::select(batch_overall, TF, Dose) |> filter(!(TF %in% c("D0", "D0_confluent", "Adipo_ref", "Irf3", "Bhlhe40", "Vdr")))

# Only take TFs with high enough absolute dose
max_doses_selected <- max_doses |> filter(Dose > 3.5)
```

    [1m[22m`summarise()` has grouped output by 'batch_overall'. You can override using the
    `.groups` argument.



```R
# Only take TFs in batch2 and batch 6
tfs_oi <- max_doses_selected |> dplyr::filter(batch_overall %in% c("batch2", "batch6")) |> group_by(TF, batch_overall) |> summarise(n = n()) |> group_by(TF) |> summarise(n = n()) |> arrange(desc(n)) |> filter(n >= 2)  |> pull(TF)
```

    [1m[22m`summarise()` has grouped output by 'TF'. You can override using the `.groups`
    argument.

