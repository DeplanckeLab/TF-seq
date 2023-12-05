# TF-seq
Material and source code for the TF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.

## 1. Preprocessing of TF-seq libraries
TF-seq generates a library containing both the single-cell RNA-seq data and the TF barcodes. While it is possible to use this library only, we recommend making a second library, enriching the TF barcodes, to better assign cells to their corresponding TF barcodes (see the **Methods** section in the manuscript).

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

### 1.2. Counting TF barcodes in the enriched library
For this step, we implemented a Java tool called [TF-seq Tools](https://github.com/DeplanckeLab/TFseqTools/). Please check the dedicated GitHub page for more details. In short, we first align the R2 enriched fastq file on the vector, and then we use TF-seq Tools to count the TF-barcodes, and assign their corresponding cell barcodes.

```bash
# Aligning to the vector
STAR --runMode alignReads --genomeDir ${vector_genome} --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn TFEnrich_R2.fastq.gz
mv Aligned.out.bam TFEnrich_R2.bam
# Count the TF-barcodes
java -jar TFCounter-0.1.jar Counter -r1 TFEnrich_R1.fastq.gz -r2 TFEnrich_R2.bam -tf C3H10_10X_Metadata.txt -p BU -UMI 12 -BC 16
```
**Note:** We provide the vector genome (.fasta file) [here](vector_sequence/pSIN-TRE-TFs-3-HA-puroR_BC_final.fa). The C3H10_10X_Metadata.txt file is accessible [here](metadata/C3H10_10X_Metadata.txt). Note that for the manuscript, we did not use the whole C3H10_10X_Metadata.txt matrix, but subsetted it for each 10x library (exp5 to exp13).

### 1.3. Robustly assigning each cell to a TF
At this point, we get 
- a scRNA-seq count matrix from cellranger (**1.1**) with automated filtered cells based on cellranger cutoffs for detecting empty cells
- a TF-Cell mapping matrix from TF-seq Tools (**1.2**) which detected TF barcodes and their corresponding cell barcodes (we used the read count matrix, not the UMI count matrix).

So now we take the overlapping cell barcodes between these two matrices, in order to build the final scRNA-seq object (both containing the scRNA-seq counts from 10x, and the TF assignment as a metadata).
But we quickly realized that we needed an extra filtering step on this final TF-Cell mapping matrix, because even though most of cells were assigned a single TF barcode (as expected), some of them were aggregating reads from multiple TF barcodes (which should not happen). So we decided to clearly filter these latter out.
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

### 1.4. Filtering outlier cells

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
  
  ### --- Mitochondrial
  mito.genes <- read.table("~/SVRAW1/prainer/Files/Mouse/data.annot/Mus_musculus.GRCm38.96_mito.annot.txt")
  mito.genes <- mito.genes$ens_id[mito.genes$ens_id %in% rownames(data.seurat)]
  # Calculating the ratio of mitochondrial reads for each cell
  data.seurat$percent.mito <- data.seurat[mito.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
  ### --- Ribosomal
  ribo.genes <- read.table("~/SVRAW1/prainer/Files/Mouse/data.annot/Mus_musculus.GRCm38.91_rRNA.annot.txt")
  ribo.genes <- ribo.genes$ens_id[ribo.genes$ens_id %in% rownames(data.seurat)]
  # Calculating the ratio of ribosomal reads for each cell
  data.seurat$percent.rRNA <- data.seurat[ribo.genes, ]$nCount_RNA/data.seurat$nCount_RNA*100
  
  ### --- Protein Coding
  protCod.genes <- read.table("~/SVRAW1/prainer/TF_scRNAseq_04.2019/Metadata/GRCm38.96_Vector_data.annot.txt", sep = "\t")
  protCod.genes <- subset(protCod.genes, biotype == "protein_coding")
  protCod.genes <- protCod.genes$ens_id[protCod.genes$ens_id %in% rownames(data.seurat)]
  # Calculating the ratio of protein coding reads for each cell
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
```

### 1.5. Calculating functional cells

TODO

### 1.6. Final dataset

For the final released dataset, we kept only TFs with more than 8 cells, and we assigned the internal "mCherry-BCXX" TF barcodes to D0 cells (non differentiated MSCs) and MatureAdipo cells (based on adiposcore).

## 2. Manuscript Figures
### Figure 1 A
TODO
