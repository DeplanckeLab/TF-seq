
combinations <- list(
  c("Pparg", "Runx2"), #1
  c("Cebpa", "Pparg"), #2
  c("Cebpa", "Mycn"), #3
  c("Cebpa", "Myog"), #4
  c("Mycn", "Myog"), #5
  c("Mycn", "Runx2"), #6
  c("Mycn", "Pparg") #7
)


# First, ensure you have the necessary libraries installed and loaded
# install.packages(c("Seurat", "edgeR"))
library(Seurat)
library(edgeR)

#' Perform Mutual Exclusion Differential Expression using edgeR-QLF
#'
#' This function performs two separate differential expression tests based on two metadata vectors.
#' 1. DE between cells with 'low' vs 'high' values in vector1, excluding cells with 'high' values in vector2.
#' 2. DE between cells with 'low' vs 'high' values in vector2, excluding cells with 'high' values in vector1.
#'
#' It uses the edgeR QL-F test, which is recommended for scRNA-seq data.
#'
#' @param seurat_object A Seurat object containing the count data and metadata.
#' @param assay The name of the assay to use for DE (e.g., "RNA", "oeumi"). Default is "RNA".
#' @param vector1_col The name of the metadata column for the first vector (e.g., "vector1").
#' @param vector2_col The name of the metadata column for the second vector (e.g., "vector2").
#' @param threshold The numeric threshold to define 'low' vs 'high' groups. Default is 5.
#'
#' @return A data frame containing merged differential expression results.
#'         Columns are suffixed with _1 or _2 to denote the comparison.
#'         Includes logFC, PValue, and FDR for each comparison.
#'
#' @examples
#' \dontrun{
#' # Load example data (pbmc_small is small and convenient)
#' library(SeuratData)
#' InstallData("pbmc3k")
#' data("pbmc3k")
#' pbmc3k_seu <- pbmc3k
#'
#' # Create the metadata vectors based on two markers, e.g., MS4A1 (B-cells) and CD3D (T-cells)
#' # Note: We use the raw counts from the "RNA" assay for this example
#' pbmc3k_seu$vector1 <- GetAssayData(pbmc3k_seu, assay = "RNA", slot = "counts")["MS4A1",]
#' pbmc3k_seu$vector2 <- GetAssayData(pbmc3k_seu, assay = "RNA", slot = "counts")["CD3D",]
#'
#' # Run the function with a threshold of 0 (present vs absent)
#' de_results <- performMutualExclusionDE(seurat_object = pbmc3k_seu,
#'                                        assay = "RNA",
#'                                        vector1_col = "vector1",
#'                                        vector2_col = "vector2",
#'                                        threshold = 0)
#'
#' # View the top results
#' head(de_results)
#'
#' # Find genes that are highly upregulated in MS4A1+ cells but not in CD3D+ cells
#' library(dplyr)
#' de_results %>%
#'   filter(FDR_1 < 0.05, logFC_1 > 1) %>%
#'   arrange(FDR_1) %>%
#'   head()
#' }
performMutualExclusionDE <- function(seurat_object,
                                              assay = "RNA",
                                              vector1_col = "vector1",
                                              vector2_col = "vector2",
                                              threshold = 5,
                                              filter_genes = TRUE) { # New parameter

  run_edger_qlf <- function(counts, group, filter) {
    if (length(levels(group)) < 2) {
      warning("Only one group found. Skipping DE analysis.")
      return(NULL)
    }
    
    dge <- DGEList(counts = counts, group = group)
    design <- model.matrix(~ group)

    # --- MODIFIED PART ---
    if (filter) {
      message("Filtering low-expression genes (recommended).")
      keep <- filterByExpr(dge, design, min.count = 1, min.prop=0.05, min.total.count = 10)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
    } else {
      message("Skipping gene filtering step (not recommended).")
    }
    # --- END MODIFIED PART ---

    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmQLFit(dge, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = 2)
    results <- topTags(qlf, n = Inf)$table
    return(results)
  }

  # ... (The rest of the main function logic is identical to before)
  # It will call the modified run_edger_qlf helper function
  
  message("Fetching metadata...")
  meta <- seurat_object@meta.data

  if (!all(c(vector1_col, vector2_col) %in% colnames(meta))) {
    stop("One or both specified vector columns are not in the Seurat object's metadata.")
  }

  # Comparison 1
  message(paste0("\n--- Starting Comparison 1: DE based on '", vector1_col, "' ---"))
  cells_for_comp1 <- rownames(meta[meta[[vector2_col]] <= threshold, ])
  seu_subset1 <- subset(seurat_object, cells = cells_for_comp1)
  group1 <- factor(ifelse(seu_subset1@meta.data[[vector1_col]] < threshold, "low", "high"), levels = c("low", "high"))
  counts1 <- GetAssayData(seu_subset1, assay = assay, slot = "counts")
  results1 <- run_edger_qlf(counts = counts1, group = group1, filter = filter_genes)
  if (!is.null(results1)) { colnames(results1) <- paste0(colnames(results1), "_1") }

  # Comparison 2
  message(paste0("\n--- Starting Comparison 2: DE based on '", vector2_col, "' ---"))
  cells_for_comp2 <- rownames(meta[meta[[vector1_col]] <= threshold, ])
  seu_subset2 <- subset(seurat_object, cells = cells_for_comp2)
  group2 <- factor(ifelse(seu_subset2@meta.data[[vector2_col]] < threshold, "low", "high"), levels = c("low", "high"))
  counts2 <- GetAssayData(seu_subset2, assay = assay, slot = "counts")
  results2 <- run_edger_qlf(counts = counts2, group = group2, filter = filter_genes)
  if (!is.null(results2)) { colnames(results2) <- paste0(colnames(results2), "_2") }

  # Merge Results
  message("\n--- Merging results ---")
  if (is.null(results1) && is.null(results2)) { return(data.frame()) }
  if (is.null(results1)) { results2$gene <- rownames(results2); return(results2) }
  if (is.null(results2)) { results1$gene <- rownames(results1); return(results1) }
  merged_df <- merge(results1, results2, by = "row.names", all = TRUE)
  colnames(merged_df)[1] <- "gene"
  
  message("Done!")
  return(merged_df)
}