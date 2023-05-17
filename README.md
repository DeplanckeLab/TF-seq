# TF-seq
Material and source code for TF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010)

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
**Note:** In the manuscript, we used the GRCm38 genome assembly (mm10) from Ensembl (release 96) created using `cellranger mkref`. But any 10x genome would work here.

### 1.2. Assign each cell to a TF
For this step, we implemented a Java tool called [TF-seq Tools](https://github.com/DeplanckeLab/TFseqTools/). Please check the dedicated GitHub page for more details. In short, we first align the R2 enriched fastq file on the vector, and then we use TF-seq Tools to count the TF-barcodes, and assign their corresponding cell barcodes.

```bash
# Aligning to the vector
STAR --runMode alignReads --genomeDir ${vector_genome} --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn TFEnrich_R2.fastq.gz
mv Aligned.out.bam TFEnrich_R2.bam
# Count the TF-barcodes
java -jar TFCounter-0.1.jar Counter -r1 TFEnrich_R1.fastq.gz -r2 TFEnrich_R2.bam -tf C3H10_10X_Metadata.txt -p BU -UMI 12 -BC 16
```
**Note:** We provide the vector genome (.fasta file) [here](vector_genome_assembly/pSIN-TRE-TFs-3-HA-puroR_BC_final.fa). The C3H10_10X_Metadata.txt file is accessible [here](metadata/C3H10_10X_Metadata.txt).

### 1.3. Filtering the TF matrix
**TODO:** This is the code of Pernille?

## 2. Manuscript Figures
### Figure 1 A
bla
