# TF-seq
Material and source code for TF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-XXXX](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXX)

## 1. Preprocessing
TF-seq generates a library containing both the single-cell RNA-seq data and the TF barcodes. While it is possible to run it only once, we recommend making a second library, enriching the TF barcodes, to better assign cells to their corresponding TF barcodes.

Therefore, you would end up with two libraries:
- A standard 10x library, to be analyzed with cellranger
- An enriched library, containing cell barcodes and TF barcodes, that you will use for assigning TF to cells.

### 1.1. Running CellRanger on the 10x library
Here we follow standard protocols. We only add a trimming step before, using cutadapt, to get rid of TSO adapters and polyA fragments. The complete code is thus as follows:

```bash
cutadapt -m 20 -G ^AAGCAGTGGTATCAACGCAGAGTACATGGG -B "A{150}" -o tfseq_trimmed_L001_R1_001.fastq -p tfseq_trimmed_L001_R2_001.fastq tfseq_L001_R1_001.fastqz $tfseq_L001_R2_001.fastqz

gzip tfseq_trimmed_L001_R1_001.fastq
gzip tfseq_trimmed_L001_R2_001.fastq

cellranger count --id tfseq_trimmed --fastqs=./ --sample=tfseq_trimmed --transcriptome=${whatever_10x_genome_plus_vector} --nosecondary
```

### 1.2. Assign each cell to a TF
For this we implemented a Java tool, called [TF-seq Tools](https://github.com/DeplanckeLab/TFseqTools/). Please check the dedicated GitHub page for this step.
In short, we align the R2 enriched fastq file on the vector, and then we use TF-seq Tools to count the TF-barcodes, and their mapping to the cell barcode.

```bash
# Aligning to the vector
STAR --runMode alignReads --genomeDir ${vector_genome} --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --readFilesIn TFEnrich_R2.fastq.gz
mv Aligned.out.bam TFEnrich_R2.bam
# Count the TF-barcodes
java -jar TFCounter-0.1.jar Counter -r1 TFEnrich_R1.fastq.gz -r2 TFEnrich_R2.bam -tf TF_barcodes.txt -p BU -UMI 12 -BC 16
```

Note: We provide the vector (.fasta file) here. The TF_barcodes.txt file is accessible from ArrayExpress.

### 1.3.
