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


```R
# tfs_oi <- c("Pparg", "Mycn")
# tfs_oi <- c("Klf4")
```


```R
select_tf_dataset <- function(seu, tf, batches = NULL) {
    filtered_tf <- seu@meta.data |> filter((TF == !!tf) & (Phase == "G1"))
    if (is.null(batches)) {
        batches_oi <- unique(filtered_tf$batch_overall)
    } else {
        batches_oi <- batches
    }
    filtered_control <- seu@meta.data |> filter((TF %in% c("D0", "D0_confluent")) & (Phase == "G1") & (batch_overall %in% batches_oi))
    filtered_tf <- filtered_tf |> filter(batch_overall %in% batches_oi)
    filtered <- rbind(filtered_control, filtered_tf)
    seu <- seu |> subset(cells = rownames(filtered))
    return(seu)
}


select_tf_datasets <- function(seu, tf, batches = NULL) {
    filtered_tf <- seu@meta.data |> filter((TF == !!tf) & (Phase == "G1"))
    if (is.null(batches)) {
        batches_oi <- unique(filtered_tf$batch_overall)
    } else {
        batches_oi <- batches
    }
    filtered_control <- seu@meta.data |> filter((TF %in% c("D0", "D0_confluent")) & (Phase == "G1") & (batch_overall %in% batches_oi))
    filtered_tf <- filtered_tf |> filter(batch_overall %in% batches_oi)
    filtered <- rbind(filtered_control, filtered_tf)

    datasets <- list()
    for (batch in batches_oi) {
        filtered_oi <- filtered |> filter(batch_overall == !!batch)
        dataset <- seu |> subset(cells = rownames(filtered_oi))
        datasets[[batch]] <- dataset
    }
    return(datasets)
}

select_tf_datasets_batch <- function(seu, tf, batches = NULL) {
    filtered_tf <- seu@meta.data |> filter((TF == !!tf) & (Phase == "G1"))
    if (is.null(batches)) {
        batches_oi <- unique(filtered_tf$batch)
    } else {
        batches_oi <- batches
    }
    filtered_control <- seu@meta.data |> filter((TF %in% c("D0", "D0_confluent")) & (Phase == "G1") & (batch %in% batches_oi))
    filtered_tf <- filtered_tf |> filter(batch %in% batches_oi)
    filtered <- rbind(filtered_control, filtered_tf)

    datasets <- list()
    for (batch in batches_oi) {
        filtered_oi <- filtered |> filter(batch == !!batch)
        dataset <- seu |> subset(cells = rownames(filtered_oi))
        datasets[[batch]] <- dataset
    }
    return(datasets)
}
```

```R
mappings_tfs <- list()

for (tf in tfs_oi) {
    print(tf)

    # filter on batches with high dose
    batches <- max_doses_selected |> filter(TF == tf) |> pull(batch_overall)
    seu_tf <- select_tf_dataset(seu, tf, batches = batches)
    seus_tf <- select_tf_datasets(seu, tf, batches = batches)

    DefaultAssay(seu_tf) <- "corrected"

    # take the union of 1000 variable genes
    features <- c()
    for (batch in batches) {
        stds <- apply(seu_tf[,seu_tf$batch_overall == batch]@assays$corrected@data, 1, sd) 
        features <- c(features, rownames(seu_tf@assays$corrected@data)[order(-stds)][1:1000])
    }
    VariableFeatures(seu_tf) <- unique(features)
    # VariableFeatures(seu_tf) <- rownames(seu_tf)

    layer <- "corrected"
    spline_values_batches <- map(seus_tf, function(seu) {
        Y <- t(seu@assays[[layer]]@data[VariableFeatures(seu_tf), ])
        x <- seu$Dose

        # xmax <- max(x)
        xmax <- quantile(x[x>0], 0.9)

        knots <- c(seq(0, xmax, 0.5), xmax)
        
        x_desired <- c(seq(0., xmax,by = 0.025), xmax)
        # x_desired <- c(seq(0., xmax,length.out=100), xmax)
        
        spline_values <- apply(Y, 2, function(y) {
            if (sum(y > 0) < 3) {
                return(rep(0, length(x_desired)))
            }
            spline <- smooth.spline(x, y, all.knots = knots, tol = 0.5)
            predict(spline, x_desired)$y
        })

        return(list(y = spline_values, x = x_desired))
    })

    # combinations = data.frame(t(combn(names(spline_values_batches), 2)))
    combinations = expand.grid(names(spline_values_batches), names(spline_values_batches))
    colnames(combinations) <- c("batch1", "batch2")
    combinations$oi <- map2_lgl(combinations$batch1, combinations$batch2, function(batch1, batch2) {
        # check whether batch1 is earlier than batch 2 in the levels
        return (which(names(spline_values_batches) == batch1) > which(names(spline_values_batches) == batch2))
    })
    combinations <- combinations |> dplyr::filter(oi)

    mapping <- function(scaling, x1, x2, y1, y2){
        x2_scaled <- x2 * scaling
        mapping <- map_int(x2_scaled, ~which.min(abs(x1 - .x)))
        dist1 <- sum(rowMeans((y1[mapping, ] - y2[, ])**2))

        x1_scaled <- x1 / scaling
        mapping <- map_int(x1_scaled, ~which.min(abs(x2 - .x)))
        dist2 <- sum(rowMeans((y2[mapping, ] - y1[, ])**2))

        return(dist1 + dist2)
    }

    combinations_oi <- combinations |> dplyr::filter(oi)

    mappings <- map2_dfr(combinations_oi$batch1, combinations_oi$batch2, function(batch1, batch2) {
        y1 <- spline_values_batches[[batch1]]$y
        y2 <- spline_values_batches[[batch2]]$y

        # center at Dose 0
        y1 <- sweep(y1, 2, y1[1, ], "-")
        y2 <- sweep(y2, 2, y2[1, ], "-")

        x1 <- spline_values_batches[[batch1]]$x
        x2 <- spline_values_batches[[batch2]]$x

        scalers <- 2**(seq(-2., 2., length.out = 100))
        distances <- purrr::map_dbl(scalers, mapping, x1 = x1, x2 = x2, y1 = y1, y2 = y2)
        plotdata <- tibble(scaling = scalers, distance = distances, batch1 = batch1, batch2 = batch2)
    })

    mappings_tfs[[tf]] <- mappings
}
```
