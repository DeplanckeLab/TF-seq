# TF-seq
Material and source code for TF-seq manuscript

## 1. Preprocessing
TF-seq generates a library containing both the single-cell RNA-seq data and the TF barcodes. While it is possible to run it only once, we recommend making a second library, enriching the TF barcodes, to better assign cells to their corresponding TF barcodes.

Therefore, you would end up with two libraries:
- A standard 10x library, to be analyzed with cellranger
- An enriched library, containing cell barcodes and TF barcodes, that you will use for assigning TF to cells.

### 1.1. Running CellRanger on the 10x library
Here we follow standard protocols. We only add a trimming step before, using cutadapt, to get rid of some adapters and polyA fragments. The complete code is thus as follows:

```bash
cutadapt -m 20 -G ^AAGCAGTGGTATCAACGCAGAGTACATGGG -B "A{150}" -o tfseq_trimmed_L001_R1_001.fastq -p tfseq_trimmed_L001_R2_001.fastq tfseq_L001_R1_001.fastqz $tfseq_L001_R2_001.fastqz

gzip tfseq_trimmed_L001_R1_001.fastq
gzip tfseq_trimmed_L001_R2_001.fastq

cellranger count --id tfseq_trimmed --fastqs=./ --sample=tfseq_trimmed --transcriptome=${whatever_10x_genome_plus_vector} --nosecondary
```

### 1.2. Assign each cell to a TF
For this we implemented a Java tool, called [TF-seq Tools](https://github.com/DeplanckeLab/TFseqTools/). Please check the dedicated GitHub page for this step.

### 1.1. Download software
TF-seq preprocessing tool is provided as a [single executable jar file](../master/releases/TFseqTools-1.0.jar?raw=true).
The .jar file contains all required materials and can be run on any terminal.
