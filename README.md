# scTF-seq
Material and source code for the scTF-seq manuscript. 

📰Manuscript: TBD

📜Preprint: [10.1101/2024.01.30.577921v1](www.biorxiv.org/content/10.1101/2024.01.30.577921v1) 

📅Data: [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010)

📊Data explorer: [ASAP96](https://asap.epfl.ch/projects/ASAP96)


## Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.

**⏭ Preprocessing of raw data to a usuable Seurat object can be found in [code/01_Preprocessing_and_integration](code/01_Preprocessing_and_integration/).**


## Notebooks/scripts that reproduce each figure panel

<table border="1" cellspacing="100" cellpadding="0" width="100%">
  <thead>
    <tr>
      <th>Panel</th>
      <th>Script</th>
    </tr>
  </thead>
  <tbody>
    <tr><th colspan="2" border="1">Figure 1</th></tr>
    <tr><td>1<strong>e-f</strong></td><td><a href="code/01_Preprocessing_and_integration">01_Preprocessing_and_integration</a></td></tr>
    <tr><th colspan="2" border="1">Figure 2</th></tr>
    <tr><td>1<strong>a-d</strong></td><td><a href="code/02_Clustering_DE_and_GSEA">02_Clustering_DE_and_GSEA</a></td></tr>
    <tr><th colspan="2" border="1">Figure 3</th></tr>
    <tr><td>3<strong>a,f-g</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity">04_Dose_sensitivity_and_reprogramming_capacity</a></td></tr>
    <tr><td>3<strong>b-e</strong></td><td><a href="code/07_Endogenous_dose/2-create_figures.ipynb">07_Endogenous_dose/2-create_figures.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Figure 4</th></tr>
    <tr><td>4</td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/1.3_protein_features_enrichment.R">04_Dose_sensitivity_and_reprogramming_capacity/1.3_protein_features_enrichment.R</a></td></tr>
    <tr><th colspan="2" border="1">Figure 5</th></tr>
    <tr><td>5<strong>a-m</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity">04_Dose_sensitivity_and_reprogramming_capacity</a></td></tr>
    <tr><th colspan="2" border="1">Figure 6</th></tr>
    <tr><td>6<strong>a-e</strong></td><td><a href="code/09_Cell_cycle">09_Cell_cycle</a></td></tr>
    <tr><th colspan="2" border="1">Figure 7</th></tr>
    <tr><td>7<strong>a</strong></td><td><a href="code/10_Combinations/0-prepare_main_seurat.ipynb">10_Combinations/0-prepare_main_seurat.ipynb</a></td></tr>
    <tr><td>7<strong>b-c</strong></td><td><a href="code/10_Combinations/4-state_bias.ipynb">10_Combinations/4-state_bias.ipynb</a></td></tr>
    <tr><td>7<strong>d</strong></td><td><a href="code/10_Combinations/3-pairwise_dosage.ipynb">10_Combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td>7<strong>e</strong></td><td><a href="code/10_Combinations/5-unique_state_enrichment.ipynb">10_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 1</th></tr>
    <tr><td>ED1<strong>a-b</strong></td><td><a href="code/01_Preprocessing_and_integration/">01_Preprocessing_and_integration/</a></td></tr>
    <tr><td>ED1<strong>c</strong></td><td><a href="code/08_RNAscope_quantification/1.1_plot_qupath_measurements_mCherry.R">08_RNAscope_quantification/1.1_plot_qupath_measurements_mCherry.R</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 2</th></tr>
    <tr><td>ED2<strong>a-d</strong></td><td><a href="code/05_Other_overexpression_studies/dose_comparison.py">05_Other_overexpression_studies/dose_comparison.py</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 3</th></tr>
    <tr><td>ED3<strong>a-d</strong></td><td><a href="code/02_Clustering_DE_and_GSEA/">02_Clustering_DE_and_GSEA</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 4</th></tr>
    <tr><td>ED4<strong>a-b,d</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/1.3_protein_features_enrichment.R">04_Dose_sensitivity_and_reprogramming_capacity/1.3_protein_features_enrichment.R</a></td></tr>
    <tr><td>ED4<strong>c</strong></td><td><a href="code/07_Endogenous_dose/2-create_figures.ipynb">07_Endogenous_dose/2-create_figures.ipynb</a></td></tr>
    <tr><td>ED4<strong>e</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb">04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb</a></td></tr>
    <tr><td>ED4<strong>f</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb">04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb</a></td></tr>
    <tr><td>ED4<strong>g</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb">04_Dose_sensitivity_and_reprogramming_capacity/Power/power.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 6</th></tr>
    <tr><td>ED6<strong>a-f</strong></td><td><a href="code/02_Clustering_DE_and_GSEA">02_Clustering_DE_and_GSEA</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 7</th></tr>
    <tr><td>ED7<strong>a-l</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity">04_Dose_sensitivity_and_reprogramming_capacity</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 8</th></tr>
    <tr><td>ED7<strong>a-e</strong></td><td><a href="code/09_Cell_cycle">09_Cell_cycle</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 9</th></tr>
    <tr><td>ED9<strong>a</strong></td><td><a href="code/10_Combinations/5-unique_state_enrichment.ipynb">10_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><td>ED9<strong>b</strong></td><td><a href="code/10_Combinations/3-pairwise_dosage.ipynb">10_Combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td>ED9<strong>f</strong></td><td><a href="code/10_Combinations/5-unique_state_enrichment.ipynb">10_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Supplementary Figure 3</th></tr>
    <tr><td>S3<strong>a-g</strong></td><td><a href="code/03_Reproducibility/reproducibility.ipynb">03_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>h-k</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/Scaling/scaling.ipynb">04_Dose_sensitivity_and_reprogramming_capacity/Scaling/scaling.ipynb</a></td></tr>
    <tr><td>S3<strong>i-j</strong></td><td><a href="code/04_Dose_sensitivity_and_reprogramming_capacity/Scaling/scaling.ipynb">04_Dose_sensitivity_and_reprogramming_capacity/Scaling/scaling.ipynb</a></td></tr>
  </tbody>
</table>
