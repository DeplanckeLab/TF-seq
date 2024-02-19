# scTF-seq
Material and source code for the scTF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.

## 1. Preprocessing of scTF-seq libraries
scTF-seq generates a library containing both the single-cell RNA-seq data and the TF barcodes. While it is possible to use this library only, we recommend making a second library, enriching the TF barcodes, to better assign cells to their corresponding TF barcodes (see the **Methods** section in the manuscript).

Therefore, a typical run output two libraries:
- A standard 10x library, to be analyzed with cellranger
- An enriched library, containing cell barcodes and TF barcodes, that is used for assigning TF to cells.

### 1.1. Running CellRanger on the 10x library
Here we follow standard protocols. We only add a trimming step before, using cutadapt, to get rid of TSO adapters and polyA fragments. The complete code is thus as follows:

```bash
cutadapt -m 20 -G ^AAGCAGTGGTATCAACGCAGAGTACATGGG -B "A{150}" -o tfseq_trimmed_L001_R1_001.fastq -p tfseq_trimmed_L001_R2_001.fastq tfseq_L001_R1_001.fastqz $tfseq_L001_R2_001.fastqz

gzip tfseq_trimmed_L001_R1_001.fastq
gzip tfseq_trimmed_L001_R2_001.fastq

cellranger count --id tfseq_trimmed --fastqs=./ --sample=tfseq_trimmed --transcriptome=${10x_genome} --nosecondary
```
**Note:** In the manuscript, we used the GRCm38 genome assembly (mm10) from Ensembl (release [96](https://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz)) created using `cellranger mkref`. Specifically, we created an assembly by merging the GRCm38.96 mouse genome with the [vector sequence](https://github.com/DeplanckeLab/TF-seq/blob/main/vector_sequence/pSIN-TRE-TFs-3-HA-puroR_BC_final.fa), which was then used to quantify the amount of vector integrated into the cell. This "Vector" abundance was then termed as "Dose" in the manuscript, as a proxy to TF overexpression abundance.

**Note 2:** Experiments 12 and 13 were actually two sequencing runs of the same library. So, after generating the cellranger outputs independently for exp12 and exp13, we merged them into a unique library using `cellranger aggr`, that we named **exp12-13**

### 1.2. Counting TF barcodes in the enriched library <sub>(see full code here: [[sh](code/1.2_Counting_TF_barcodes.sh)])</sub>
For this step, we implemented a Java tool called [TF-seq Tools](https://github.com/DeplanckeLab/TFseqTools/). Please check the dedicated GitHub page for more details. In short, we first align the R2 enriched fastq file on the vector, and then we use TF-seq Tools to count the TF-barcodes, and assign their corresponding cell barcodes.

```bash
# Aligning to the vector
STAR --runMode alignReads --genomeDir ${vector_genome} --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn TFEnrich_R2.fastq.gz
mv Aligned.out.bam TFEnrich_R2.bam
# Count the TF-barcodes
java -jar /software/TFseqTools-1.1.jar Counter --r1 TFEnrich_R1.fastq.gz --r2 TFEnrich_R2.bam --tf C3H10_10X_Metadata.txt -p BU --UMI 12 --BC 16
```
**Note:** We provide the vector genome (.fasta file) [here](vector_sequence/pSIN-TRE-TFs-3-HA-puroR_BC_final.fa). The C3H10_10X_Metadata_expXX.txt files are accessible [here](metadata/).

**Note 2:** Barcode length was kept as 16 for all experiments. However, UMI length was 10 for exp05 and exp06, and 12 for all other experiments.

### 1.3. Robustly assigning each cell to a TF <sub>(see full code here: [[Rmd](code/1.3_Assigning_TF.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.3_Assigning_TF.html)])</sub>
At this point, we get 
- a scRNA-seq count matrix from cellranger (**1.1**) with automated filtered cells based on cellranger cutoffs for detecting empty cells
- a TF-Cell mapping matrix (read count matrix or UMI count matrix) from TF-seq Tools (**1.2**) which detected TF barcodes and their corresponding cell barcodes. Although the read count matrix correlated well (Pearson's r = 0.9) with the UMI count matrix, the latter yielded ~40% less cells associated with TF-IDs. So we used the read count matrix not the UMI count matrix for **1.3** to keep more cells at this step. The UMI count matrix was used for calculating TF dose as described below (**1.5**)

So now we take the overlapping cell barcodes between these two matrices, in order to build the final scRNA-seq object (both containing the scRNA-seq counts from 10x, and the TF assignment as a metadata).
But we quickly realized that we needed an extra filtering step on this final TF-Cell mapping matrix, because even though most of cells were assigned a single TF barcode (as expected), some of them were aggregating reads from multiple TF barcodes (which should not happen). So we decided to clearly filter these latter out.

#### 1.3.1 Experiments 5-11
For this, we implemented a KneePoint algorithm, to only keep the cells where the main TF was clearly majoritarily present in the cell.

```R
# Cellranger matrix (output of 1.1)
h5f = rhdf5::H5Fopen(name = "filtered_feature_bc_matrix.h5", flags = "H5F_ACC_RDONLY")
data.barcodes_cellranger <- as.data.frame(as.matrix(h5f$matrix))$V1$barcodes
data.barcodes_cellranger <- gsub(pattern="-1", replacement = "", x = data.barcodes_cellranger)
# TF-cell mapping matrix (output of 1.2)
data.tf_read <- data.table::fread("C3H10_10X_enriched.read_matrix.txt", data.table = F)

# Filtering TF-cell mapping matrix, keeping only overlapping barcodes with 10x cellranger matrix
common.barcodes <- intersect(colnames(data.tf_read)[3:ncol(data.tf_read)], data.barcodes_cellranger)
data.tf_read <- data.tf_read[,common.barcodes] # Filtered Matrix
# Calculate rate of main TF
nCounts <- colSums(data.tf_read)
nMax <- apply(data.tf_read, 2, max) 
TFs.MaxRate <- nMax / nCounts
TFsBC.MaxRate.cutoff <- TFs.MaxRate[nCounts > 5] # Select cells with at least 5 reads
# Knee point cutoff detection
TFsBC.MaxRate.cutoff_order <- TFsBC.MaxRate.cutoff[order(TFsBC.MaxRate.cutoff, decreasing = T)] # order from 1 to 0 
cutoff <- SamSPECTRAL::kneepointDetection(TFsBC.MaxRate.cutoff_order)

# Cell barcodes passing filtering
cell_barcodes_filtered <- names(TFsBC.MaxRate.cutoff)[TFsBC.MaxRate.cutoff > TFsBC.MaxRate.cutoff_order[cutoff$MinIndex]]
```
After filtering out these cells with potentially ambiguous TF attribution, we robustly assigned each TF to their cells by selecting the most abundant TF.

#### 1.3.2 Experiments 12-13
Experiments 12-13 were a bit special since some cells could contain a combination of two TFs. Therefore the previously defined kneepoint algorithm cannot work.
In this case, we implemented another algorithm, in two parts, for identifying singlets, doublets and cell combinations. 

First, we sum the two top TFs for each cell barcode, to create a "pseudo TF", where we apply the kneepoint algorithm. Any non-selected TF is then filtered out. Next, we check which of the remaining TFs could be singlets, by focusing solely on the most abundant TF, using the same kneepoint algorithm. Finally, once singlets and "potential combination" are isolated, we check which of the "potential combination" were used in the experiment (see [here](metadata/C3H10_10X_Metadata_exp12-13_combination.txt)), and we filtered out all the "Doublets" that were not intended combinations.

```R
## Cellranger matrix (output of 1.1)
data.cell_ranger <- Seurat::Read10X_h5(filename = "filtered_feature_bc_matrix.h5", use.names = F)
data.barcodes <- unlist(data.cell_ranger@Dimnames[2])

# TF-cell mapping matrix (output of 1.2)
data.tf_read <- data.table::fread("C3H10_10X_enriched.read_matrix.txt", data.table = F)
rownames(data.tf_read) <- data.tf_read$TFName
data.barcodes_tf <- colnames(data.tf_read)[3:ncol(data.tf_read)]

# Overlap
common.barcodes <- intersect(data.barcodes_tf, data.barcodes)
data.tf_read <- data.tf_read[,common.barcodes] # Filtered TF Matrix

# Get combinations
meta_combinations <- data.table::fread("C3H10_10X_Metadata_exp12-13_combination.txt", data.table = F)
rownames(meta_combinations) <- paste0(meta_combinations$TF1,"-",meta_combinations$TF2)

# Precompute total TF read counts per cell barcode
nCounts <- colSums(data.tf_read)

# Calculate rate of main and second TF using Rfast package
all.top2 <- Rfast::colnth(as.matrix(data.tf_read), rep(2,ncol(data.tf_read)), descending = T, num.of.nths = 2)
first.indexes <- Rfast::colnth(as.matrix(data.tf_read), index.return = T, rep(1,ncol(data.tf_read)), descending = T, num.of.nths = 1) # Can't do both?
second.indexes <- Rfast::colnth(as.matrix(data.tf_read), index.return = T, rep(2,ncol(data.tf_read)), descending = T, num.of.nths = 1) # Can't do both?

# Prepare output dataframe with rate of main and second TF
data.tf_read.top2 <- data.frame(matrix(nrow = ncol(data.tf_read), ncol = 6), row.names = colnames(data.tf_read))
colnames(data.tf_read.top2) <- c("TFmax", "TFmax.rate", "nMax", "TF2nd", "TF2nd.rate", "n2nd")

# Fill it
data.tf_read.top2$nMax <- all.top2[1,]
data.tf_read.top2$TFmax.rate <- data.tf_read.top2$nMax / nCounts
data.tf_read.top2$n2nd <- all.top2[2,]
data.tf_read.top2$TF2nd.rate <- data.tf_read.top2$n2nd / nCounts
data.tf_read.top2$TFmax <- rownames(data.tf_read)[first.indexes]
data.tf_read.top2$TF2nd <- rownames(data.tf_read)[second.indexes]

# Generate pseudo-TF as sum of top 2 TFs
data.tf_read.top2$nPseudoTF <- data.tf_read.top2$nMax + data.tf_read.top2$n2nd
data.tf_read.top2$PseudoTF.rate <- (data.tf_read.top2$nMax + data.tf_read.top2$n2nd) / nCounts

# Remove top barcodes if n = 0
data.tf_read.top2$TFmax[data.tf_read.top2$nMax == 0] <- NA
data.tf_read.top2$TF2nd[data.tf_read.top2$n2nd == 0] <- NA

# Which pair of TF could be a combination, based on metadata
data.tf_read.top2$potential.combination <- paste0(data.tf_read.top2$TFmax,"-",data.tf_read.top2$TF2nd) %in% rownames(meta_combinations)

# Filter cells with too low number of counts
data.tf_read.top2 <- data.tf_read.top2[names(nCounts[nCounts > 5]),,drop=F]

# Find the cutoff of the "Pseudo TF" thanks to kneepoint detection 
PseudoTF.rate_order <- data.tf_read.top2[order(data.tf_read.top2$PseudoTF.rate, decreasing = T), "PseudoTF.rate", drop=F]
cutoff <- SamSPECTRAL::kneepointDetection(PseudoTF.rate_order$PseudoTF.rate)

# Filter cells based on cutoff kneepoint pseudo TF
data.tf_read.top2 <- subset(data.tf_read.top2, PseudoTF.rate >= PseudoTF.rate_order[cutoff$MinIndex, "PseudoTF.rate"])

# Now, find the cutoff of the "Main TF", thanks to kneepoint detection to identify singlets versus combinations (real and doublets)
TFmax.rate_order <- data.tf_read.top2[order(data.tf_read.top2$TFmax.rate, decreasing = T), "TFmax.rate", drop=F]
cutoff <- SamSPECTRAL::kneepointDetection(TFmax.rate_order$TFmax.rate)

# Select cells based on cutoff kneepoint TF max rate of selected cells 
data.tf_read.top2$singles <- data.tf_read.top2$TFmax.rate >= TFmax.rate_order[cutoff$MinIndex,"TFmax.rate"]

# Remove Doublets (keep only combination and single TFs)
data.tf_read.top2 <- subset(data.tf_read.top2, potential.combination == TRUE | singles == TRUE)

# Assign TFs and combinations
data.tf_read.top2$TF2nd[data.tf_read.top2$singles] <- NA
data.tf_read.top2$TF <- meta_combinations[paste0(data.tf_read.top2$TFmax,"-",data.tf_read.top2$TF2nd), "combination"]
data.tf_read.top2$TF[is.na(data.tf_read.top2$TF)] <- data.tf_read.top2$TFmax[is.na(data.tf_read.top2$TF)]
```

### 1.4. Filtering outlier cells <sub>(see full code here: [[Rmd](code/1.4_Filtering_outlier_cells.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.4_Filtering_outlier_cells.html)])</sub>
At step 1.3 we create a Seurat object containing the raw data counts and a cell metadata containing their assigned TFs. At this stage, following the previous pipeline, all cells have an assigned TF.<br/>
Now, following the standard Seurat pipeline, we aimed at removing outlier cells given the following criteria:
- Remove outlier cells using the `isOutlier` function of the `scater` package for library depth (nCounts) and number of detected genes (nFeature)
- Remove cells with too much mitochondrial RNA or ribosomal RNA
- Remove cells with not enough protein-coding RNA

Here is the function that is performing this filtering:

```R
filtering_outlierCells <- function(seurat_path, libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75){
data.seurat <- readRDS(seurat_path) # Read Seurat object generated at step 1.3

### ---Features and Library size
# Looking for outlier cells in the nCount_RNA distribution (returns TRUE/FALSE array)
libsize.drop <- scater::isOutlier(data.seurat$nCount_RNA, nmads=libsize_nmads, type="lower", log=TRUE) # nCount_RNA / colSums(data.seurat)
# Looking for outlier cells in the nFeature_RNA distribution (returns TRUE/FALSE array)
features.drop <- scater::isOutlier(data.seurat$nFeature_RNA, nmads=features_nmads, type="lower", log=TRUE) # nFeature_RNA / as.vector(colSums(data.seurat > 0))

### --- Gene annotation
data.gene_annot <- fread("metadata/GRCm38_gene_annot.tsv", data.table = F)
  
### --- Mitochondrial
mito.genes <- subset(data.gene_annot, is_mitochondrial)$ensembl_id
mito.genes <- mito.genes[mito.genes %in% rownames(data.seurat)]
data.seurat$percent.mito <- data.seurat[mito.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
### --- Ribosomal
ribo.genes <- subset(data.gene_annot, is_ribosomal)$ensembl_id
ribo.genes <- ribo.genes[ribo.genes %in% rownames(data.seurat)]
data.seurat$percent.rRNA <- data.seurat[ribo.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
### --- Protein Coding
protCod.genes <- subset(data.gene_annot, biotype == "protein_coding")$ensembl_id
protCod.genes <- protCod.genes[protCod.genes %in% rownames(data.seurat)]
data.seurat$percent.ProtCod <- data.seurat[protCod.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100

### Running the filtering
data.seurat <- data.seurat[,features.drop < 1 & libsize.drop < 1 & data.seurat$percent.mito < max_pc_mito & data.seurat$percent.rRNA < max_pc_rRNA & data.seurat$percent.ProtCod > min_pc_protCod]
saveRDS(data.seurat, file = "...") # Saving the filtered Seurat object
}
```

We did not use the same filtering thresholds for all experiments (exp05-12), because of different qualities and biological contexts across experiments. Here is a summary of our filterings:
```R
filtering_outlierCells("exp05", libsize_nmads = 4, features_nmads = 4, max_pc_mito = 10, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp06", libsize_nmads = 4, features_nmads = 4, max_pc_mito = 10, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp07", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 25, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp08", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp09", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 35,  min_pc_protCod = 75)
filtering_outlierCells("exp10", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp11", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 15, max_pc_rRNA = 40,  min_pc_protCod = 75)
filtering_outlierCells("exp12-13", libsize_nmads = 6, features_nmads = 6, max_pc_mito = 30, max_pc_rRNA = 60,  min_pc_protCod = 75)
```

### 1.5. Final dataset

#### 1.5.1 Cell cycle scoring and Dose calculation <sub>(see full code here: [[Rmd](code/1.5.1_Cell_cycle_and_vector.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.5.1_Cell_cycle_and_vector.html)])</sub>

For each of the previously created Seurat objects at **step 1.4**, we now calculate the **Cell Cycle Phase** using the `CellCycleScoring` function of Seurat, and the *mus musculus* Cell Cycle genes downloaded from [https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv](metadata/Mus_musculus.csv).
We also calculate the **Dose** using the TF UMI matrices from the enriched libraries generated at **step 1.2**.

**Of note:** We often call the Dose 'Vector', since we use the Vector-mapped read abundance as a proxy for the Dose.

```R
# This is run for 'exp_seurat' & 'exp_tf_cell' in [exp05, exp06, exp07, exp08, exp09, exp10, exp11, exp12-13], with the corresponding Seurat object ('exp_seurat'), and TF UMI matrix from enriched library ('exp_tf_cell')
cell_cycle <- function(exp_seurat, exp_tf_cell){
  # Read Seurat object of previous step (1.4)
  data.seurat <- readRDS(exp_seurat)
  
  # Putting Vector as metadata
  data.seurat$Vector_10X <- data.seurat@assays$RNA@counts["Vector",]
  data.seurat <- data.seurat[rownames(data.seurat) != "Vector",]
  
  # TF matrix
  data.tf_read <- data.table::fread(exp_tf_cell, data.table = F)
  rownames(data.tf_read) <- data.tf_read$TFName

  # Combinations in exp12-13
  meta_combinations <- data.table::fread("metadata/C3H10_10X_Metadata_exp12-13_combination.txt", data.table = F)
  meta_combinations <- meta_combinations[!duplicated(meta_combinations$combination),]
  rownames(meta_combinations) <- meta_combinations$combination

  # Cell cycle genes (downloaded from https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv)
  cellcycle_genes <- data.table::fread("metadata/Mus_musculus.csv", sep=",", data.table = F, stringsAsFactors = F)

  # Creating the Dose metadata (Vector_UMI) by extracting the UMI values from the TF barcode matrix
  data.seurat$Vector_UMI <- NA
  for(i in 1:ncol(data.tf_read)){
    data.seurat$Vector_UMI[i] <- data.tf_read[data.seurat$TF[i],i]
    if(is.na(data.seurat$Vector_UMI[i])){
      # Combination
      tf1 <- meta_combinations[data.seurat$TF[i],"TF1"]
      tf2 <- meta_combinations[data.seurat$TF[i],"TF2"]
      data.seurat$Vector_UMI[i] <- data.tf_read[tf1,i] + data.tf_read[tf2,i]
    }
  }
  data.seurat$Log_Vector_UMI <- log(1 + data.seurat$Vector_UMI)
  
  # Normalization (required for CellCycleScoring)
  data.seurat <- Seurat::NormalizeData(data.seurat, verbose = F)
  
  # Calculate Cell Cycle score
  data.seurat <- Seurat::CellCycleScoring(data.seurat, s.features = subset(cellcycle_genes, phase == "S")$geneID, g2m.features = subset(cellcycle_genes, phase == "G2/M")$geneID)
  
  # Calculate Cell Cycle score (corrected)
  # Change threshold for Phase, Seurat assigns Phase if the max score is bigger than 0 which seems a too low threshold
  data.seurat_tmp$Phase_corrected <- ifelse(pmax(data.seurat_tmp$G2M.Score, data.seurat_tmp$S.Score) > 0.1, as.character(data.seurat_tmp$Phase), "G1")
  
  # Saving Seurat object
  saveRDS(data.seurat, file = ""...")
}
```

#### 1.5.2 D0 assignment, Adipo_ref & Myo_ref assignment, and TF renaming <sub>(see full code here: [[Rmd](code/1.5.2_D0_Ref_TF_renaming.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.5.2_D0_Ref_TF_renaming.html)])</sub>

For the final released dataset, we kept only TFs with more than 8 cells, and we assigned the internal "mCherry-BCXX" TF barcodes to their corresponding lineage and timepoints (see the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file). Then we kept only the D0 cells (non-differentiated MSCs) as control cells and reference adipocyte (Adipo_ref) and myocyte (Myo_ref) cells.<br/>

#### 1.5.3 Integration <sub>(see full code here: [[Rmd](code/1.5.3_Seurat_Integration.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.5.3_Seurat_Integration.html)])</sub>

Finally, we integrate all Seurat objects into the final atlas by using the integration tool of the Seurat package.

```R
# Here we assume 'all_exps' contains all Seurat objects for exp in [exp05, exp06, exp07, exp08, exp09, exp10, exp11, exp12-13]

# First, Normalize and HVG on each Seurat object
data.seurat_for_mnn <- list()
for(exp in all_exps){
  data.seurat <- all_exps[[exp]]
  data.seurat$batch <- exp
  data.seurat <- NormalizeData(data.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  data.seurat_for_mnn[[exp]] <- FindVariableFeatures(data.seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
}

# Then running the Seurat integration
int_features <- SelectIntegrationFeatures(object.list = data.seurat_for_mnn, verbose = F)
int_anchors <- FindIntegrationAnchors(object.list = data.seurat_for_mnn, anchor.features = int_features, verbose = F)
data.seurat_integrated <- IntegrateData(anchorset = int_anchors, verbose = F)

# Run the default Seurat pipeline on the integrated object
data.seurat_integrated <- ScaleData(data.seurat_integrated, assay = "integrated", verbose = F)
data.seurat_integrated <- RunPCA(data.seurat_integrated, assay = "integrated", npcs = 200, verbose = F)
data.seurat_integrated <- RunTSNE(data.seurat_integrated, reduction = "pca", dims = 1:200)
data.seurat_integrated <- RunUMAP(data.seurat_integrated, reduction = "pca", dims = 1:200)
data.seurat_integrated <- SetIdent(data.seurat_integrated, value = "TF")

# Save Integrated Seurat object
saveRDS(data.seurat_integrated, "...")
```
### 1.6. Calculating functional cells <sub>(see full code here: [[Rmd](code/1.6_Functional_cells.Rmd)][[html](https://htmlpreview.github.io/?https://github.com/DeplanckeLab/TF-seq/blob/main/code/1.6_Functional_cells.html)])</sub>

We calculated **functional cells**, i.e. cells that are transcriptomically different from D0 control cells, by calculating the Euclidean distance of each TF cell to the D0 control cells, in the PCA space.

**Of note:** To avoid any bias, we've run the script below for each 1) TF, for each 2) batch, for each 3) cell cycle phase, i.e. 5880 combinations.

```R
# As mentioned before, we first filter the integrated atlas created in 1.5.3 based on one of the 5880 possible combination of 1) TF, 2) batch and 3) cell cycle phase, and create the **'seurat_object'**
find_functional_cells <- function(seurat.object, nPCs = 10, threshold = 0.8, threshold.non_con = 0.8, threshold.con = 0.8){
  # Normalizing and calculating PCA on the subset object
  seurat.object <- NormalizeData(seurat.object, verbose = F) 
  seurat.object <- suppressWarnings(FindVariableFeatures(seurat.object, verbose = F, selection.method = "vst"))
  seurat.object <- ScaleData(seurat.object, verbose = F)
  n_pcs <- min(nPCs + 1, ncol(seurat.object))-1
  seurat.object <- RunPCA(seurat.object, npcs = n_pcs, features = VariableFeatures(seurat.object), verbose = F)
  pca <- seurat.object@reductions$pca@cell.embeddings
  
  # Computing the 3 centroids for all D0 in PCA space
  mean_D0s <- colMeans(pca[seurat.object$TF %in% c("D0", "D0_confluent"),])
  mean_D0s_non_confluent <- colMeans(pca[seurat.object$TF == "D0",])
  mean_D0s_confluent <- colMeans(pca[seurat.object$TF == "D0_confluent",])

  # Distance of all cells to the 3 centroids
  Dist_to_meanD0s <- apply(pca, 1, function(x){
      pos.vector <- rbind(x, mean_D0s)
      length.vectors <- dist(pos.vector, method = "euclidean", diag = F, upper = F, p = 2)[1] # p = 2, euclidean distance
      return(length.vectors)
  })
  
  Dist_to_meanD0s_non_confluent <- apply(pca, 1, function(x){
      pos.vector_non_confluent <- rbind(x, mean_D0s_non_confluent)
      length.vectors_non_confluent <- dist(pos.vector_non_confluent, method = "euclidean", diag = F, upper = F, p = 2)[1] # p = 2, euclidean distance
      return(length.vectors_non_confluent)
  })
  
  Dist_to_meanD0s_confluent <- apply(pca, 1, function(x){
      pos.vector_confluent <- rbind(x, mean_D0s_confluent)
      length.vectors_confluent <- dist(pos.vector_confluent, method = "euclidean", diag = F, upper = F, p = 2)[1] # p = 2, euclidean distance
      return(length.vectors_confluent)
  })
  
  # Compute threshold as quantile of D0 cells
  distance_thresholded_D0 <- quantile(Dist_to_meanD0s[seurat.object$TF %in% c("D0", "D0_confluent")], threshold)
  distance_thresholded_non_confluent <- quantile(Dist_to_meanD0s_non_confluent[seurat.object$TF %in% c("D0")], threshold.non_con)
  distance_thresholded_confluent <- quantile(Dist_to_meanD0s_confluent[seurat.object$TF %in% c("D0_confluent")], threshold.con)
  
  # Extract FUNCTIONAL cells, i.e. cells that pass the threshold
  cells.oi <- names(Dist_to_meanD0s)[Dist_to_meanD0s > distance_thresholded_D0]
  cells.oi_non_confluent <- names(Dist_to_meanD0s_non_confluent)[Dist_to_meanD0s_non_confluent > distance_thresholded_non_confluent]
  cells.oi_confluent <- names(Dist_to_meanD0s_confluent)[Dist_to_meanD0s_confluent > distance_thresholded_confluent]
  cells.funct <- intersect(intersect(cells.oi, cells.oi_non_confluent), cells.oi_confluent)
    
  # FUNCTIONAL cells are the non-D0 cells
  cells.funct <- names(seurat.object$TF[cells.funct] )[!seurat.object$TF[cells.funct] %in% c("D0", "D0_confluent")]
  return(cells.funct)
  }
```

## 2. Manuscript Figures

All manuscript figures are reproducible through scripts deposited [here](figures/)
