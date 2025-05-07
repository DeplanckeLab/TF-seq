# title: "Dose alignment across batches"
# author: "Wouter Saelens / Wangjie Liu"
# date: "2024-09-30"


setwd("./")

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
library(purrr)

seu <- read_rds("results/C3H10_10X_all_exps_D0regressed_integrated.rds")


select_tf_dataset <- function(seu, tf, batches = NULL) {
    filtered_tf <- seu@meta.data |> filter((TF == !!tf) & (Phase_corrected == "G1"))
    if (is.null(batches)) {
        batches_oi <- unique(filtered_tf$batch_overall)
    } else {
        batches_oi <- batches
    }
    filtered_control <- seu@meta.data |> filter((TF %in% c("D0", "D0_confluent")) & (Phase_corrected == "G1") & (batch_overall %in% batches_oi))
    filtered_tf <- filtered_tf |> filter(batch_overall %in% batches_oi)
    filtered <- rbind(filtered_control, filtered_tf)
    seu <- seu |> subset(cells = rownames(filtered))
    return(seu)
}


select_tf_datasets <- function(seu, tf, batches = NULL) {
    filtered_tf <- seu@meta.data |> filter((TF == !!tf) & (Phase_corrected == "G1"))
    if (is.null(batches)) {
        batches_oi <- unique(filtered_tf$batch_overall)
    } else {
        batches_oi <- batches
    }
    filtered_control <- seu@meta.data |> filter((TF %in% c("D0", "D0_confluent")) & (Phase_corrected == "G1") & (batch_overall %in% batches_oi))
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

# select TFs to use for alignment

max_doses <- seu@meta.data |> arrange(-Dose) |> group_by(batch_overall, TF) |> summarize(Dose = quantile(Dose, 0.9)) |> ungroup() |> dplyr::select(batch_overall, TF, Dose) |> filter(!(TF %in% c("D0", "D0_confluent", "Adipo_ref", "Irf3", "Bhlhe40", "Vdr")))

# Only take TFs with high enough absolute dose
max_doses_selected <- max_doses |> filter(Dose > 3.5)

# Only take TFs with at least 2 batches
tfs_oi <- max_doses_selected |> group_by(TF, batch_overall) |> summarise(n = n()) |> group_by(TF) |> summarise(n = n()) |> arrange(desc(n)) |> filter(n >= 2)  |> pull(TF)

# perform the actual alignment per TF
mappings_tfs <- list()

for (tf in tfs_oi) {
    print(tf)

    # filter on batches with high dose
    batches <- max_doses_selected |> filter(TF == tf) |> pull(batch_overall)
    seu_tf <- select_tf_dataset(seu, tf, batches = batches)
    seus_tf <- select_tf_datasets(seu, tf, batches = batches)

    DefaultAssay(seu_tf) <- "corrected"
    stds <- apply(seu_tf@assays$corrected@data, 1, sd)
    VariableFeatures(seu_tf) <- rownames(seu_tf@assays$corrected@data)[order(-stds)][1:2000]

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

        scalers <- 2**(seq(-1., 1., length.out = 100))
        distances <- purrr::map_dbl(scalers, mapping, x1 = x1, x2 = x2, y1 = y1, y2 = y2)
        plotdata <- tibble(scaling = scalers, distance = distances, batch1 = batch1, batch2 = batch2)
    })

    mappings_tfs[[tf]] <- mappings
}

mappings <- map2_dfr(mappings_tfs, names(mappings_tfs), function(df, name) {mutate(df, tf = name)}) |> filter(distance < 20)

mappings$pair <- paste(mappings$tf, mappings$batch1, mappings$batch2)

mappings_oi <- mappings
# mappings_oi <- mappings |> dplyr::filter(batch1 == "batch4" & batch2 == "batch3")
# mappings_oi <- mappings |> dplyr::filter(batch1 == "batch2" & batch2 == "batch1")
# mappings_oi <- mappings |> dplyr::filter(batch1 == "batch3" & batch2 == "batch1")
# mappings_oi <- mappings |> dplyr::filter(batch1 == "batch6" & batch2 == "batch3")
optimal_mappings <- mappings_oi |> arrange(distance) |> group_by(tf, batch1, batch2) |> dplyr::slice(1)

optimal_mapping <- optimal_mappings |> group_by(batch1, batch2) |> summarise(scaling = mean(scaling), mean_distance = mean(distance), n = n()) |> filter(n > 4) |> ungroup()
optimal_mapping <- bind_rows(optimal_mapping, optimal_mapping |> mutate(scaling = 1/scaling) |> dplyr::select(batch1 = batch2, batch2 = batch1, scaling, mean_distance))
optimal_mapping <- left_join(expand.grid(batch1 = unique(seu$batch_overall), batch2 = unique(seu$batch_overall)), optimal_mapping, by = c("batch1", "batch2")) |> arrange(batch1, batch2) |> mutate(scaling = ifelse(is.na(scaling), 1, scaling)) |> mutate(mean_distance = ifelse(is.na(mean_distance), 0, mean_distance))

optimal_mapping_matrix_raw <- reshape2::dcast(optimal_mapping, batch2 ~ batch1, value.var = "scaling")
optimal_mapping_matrix <- as.matrix(optimal_mapping_matrix_raw)[, 2:ncol(optimal_mapping_matrix_raw)]
rownames(optimal_mapping_matrix) <- optimal_mapping_matrix_raw$batch2
optimal_mapping_matrix <- optimal_mapping_matrix[rownames(optimal_mapping_matrix), ]
optimal_mapping_matrix <- matrix(as.numeric(optimal_mapping_matrix), nrow = nrow(optimal_mapping_matrix), ncol = ncol(optimal_mapping_matrix), dimnames = dimnames(optimal_mapping_matrix))
optimal_mapping_matrix[is.na(optimal_mapping_matrix)] <- 1.
optimal_mapping_matrix

# plot mapping for a specific pair
mappings_oi <- mappings |> dplyr::filter(batch1 == "batch6" & batch2 == "batch2")
options(repr.plot.width=10, repr.plot.height=8)
ggplot(mappings_oi) + 
    geom_line(aes(x = scaling, y = distance, color = pair)) + 
    scale_x_log10() + theme_minimal()

# correct the doses
seu$Dose_aligned <- seu$Dose
tfs <- unique(seu$TF)
tfs <- tfs[tfs != "D0" & tfs != "D0_confluent" & tfs != "Adipo_ref"]

plots <- list()

for (tf in tfs) {
    metadata_oi <- seu@meta.data |> filter(TF == tf)

    # determine correction
    batch_overalls <- unique(metadata_oi$batch_overall)
    if (length(batch_overalls) < 2) {
        next
    }

    batch_choice <- metadata_oi$batch_overall[which.max(metadata_oi$Dose)]

    corrections <- optimal_mapping_matrix[batch_choice, batch_overalls]

    print(tf)
    print(corrections)

    # divided by the correction factor
    metadata_oi$Dose_aligned <- metadata_oi$Dose / corrections[match(metadata_oi$batch_overall, batch_overalls)]

    plot <- patchwork::wrap_plots(
        ggplot(metadata_oi) + geom_boxplot(aes(x = batch_overall, y = Dose, fill = batch_overall)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
        ggplot(metadata_oi) + geom_boxplot(aes(x = batch_overall, y = Dose_aligned, fill = batch_overall)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )

    seu$Dose_aligned[rownames(metadata_oi)] <- metadata_oi$Dose_aligned

    plots[[tf]] <- plot
}


saveRDS(seu, file = "results/C3H10_10X_all_exps_D0regressed_integrated_dosealigned.rds")



## -------------------------- plot alignment for individual TFs, e.g., Runx2
plots[["Runx2"]]

seu2 <- select_tf_dataset(seu, "Runx2") |> FindVariableFeatures(nfeatures = 2000, selection.method = "mean.var.plot") |> ScaleData() |> RunPCA() |> RunUMAP(dims = 1:20)

plotdata <- seu2@meta.data
plotdata$expression <- Seurat::FetchData(seu2, vars = convert_to_ensembl(c("Runx2")))[[1]]
options(repr.plot.width=12, repr.plot.height=5)
patchwork::wrap_plots(
    ggplot(plotdata, aes(x = Dose, y = expression, color = batch_overall)) + geom_point() + geom_smooth() + theme_minimal(),
    ggplot(plotdata, aes(x = Dose_aligned, y = expression, color = batch_overall)) + geom_point() + geom_smooth() + theme_minimal()
)



