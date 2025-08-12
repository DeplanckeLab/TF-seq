# scTF-seq
Material and source code for the scTF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.


## 1. Each figure panel can be reproduced using the following scripts/notebooks

| Panel        | Script                                      |
| -----------: |  ------------------------------------------ |
| **Figure 1**                  |
| **Figure 3**                  |
| **3a**    | ?  |
| **3b**    | [ZZ_Endogenous/2-create_figures.ipynb](code/ZZ_Endogenous/2-create_figures.ipynb)  |
| **3c**    | [ZZ_Endogenous/2-create_figures.ipynb](code/ZZ_Endogenous/2-create_figures.ipynb)  |
| **3d**    | [ZZ_Endogenous/2-create_figures.ipynb](code/ZZ_Endogenous/2-create_figures.ipynb)  |
| **3e**    | [ZZ_Endogenous/2-create_figures.ipynb](code/ZZ_Endogenous/2-create_figures.ipynb)  |
| **3f**    | ?  |
| **3g**    | ?  |
| **Figure 4**                  |
| **Figure 7**                  |
| **7a**    | [ZZ_combinations/0-prepare_main_seurat.ipynb](code/ZZ_combinations/0-prepare_main_seurat.ipynb)  |
| **7b**    | [ZZ_combinations/4-state_bias.ipynb](code/ZZ_combinations/4-state_bias.ipynb)                  |
| **7c**    | [ZZ_combinations/4-state_bias.ipynb](code/ZZ_combinations/4-state_bias.ipynb)                  |
| **7d**    | [ZZ_combinations/3-pairwise_dosage](code/ZZ_combinations/3-pairwise_dosage.ipynb)                  |
| **7e**    | [ZZ_combinations/3-pairwise_dosage](code/ZZ_combinations/3-pairwise_dosage.ipynb)                  |
| **7e**    | [ZZ_combinations/5-unique_state_enrichment.ipynb](code/ZZ_combinations/5-unique_state_enrichment.ipynb)                  |
| **Extended Data Figure Figure 7**                  |
| **ed9a**    | [ZZ_combinations/5-unique_state_enrichment.ipynb](code/ZZ_combinations/5-unique_state_enrichment.ipynb)                  |
| **ed9b**    | [ZZ_combinations/3-pairwise_dosage.ipynb](code/ZZ_combinations/3-pairwise_dosage.ipynb)                  |
| **ed9e**    | ?                  |
| **ed9f**    | [ZZ_combinations/5-unique_state_enrichment.ipynb](code/ZZ_combinations/5-unique_state_enrichment.ipynb)                  |
