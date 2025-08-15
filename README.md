# scTF-seq
Material and source code for the scTF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.

**⏭ Preprocessing of raw data to a usuable Seurat object can be found in [code/01_Preprocessing_and_integration](code/01_Preprocessing_and_integration/).**


## 1. Notebooks/scripts that reproduce each figure panel

<table border="1" cellspacing="100" cellpadding="0" width="100%">
  <thead>
    <tr>
      <th>Panel</th>
      <th>Script</th>
    </tr>
  </thead>
  <tbody>
    <tr><th colspan="2" border="1">Figure 1</th></tr>
    <tr><th colspan="2" border="1">Figure 2</th></tr>
    <tr><th colspan="2" border="1">Figure 3</th></tr>
    <tr><td>3<strong>a</strong></td><td>?</td></tr>
    <tr><td>3<strong>b</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>3<strong>c</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>3<strong>d</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>3<strong>e</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>3<strong>e</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>3<strong>f</strong></td><td>?</td></tr>
    <tr><td>3<strong>g</strong></td><td>?</td></tr>
    <tr><th colspan="2" border="1">Figure 5</th></tr>
    <tr><th colspan="2" border="1">Figure 4</th></tr>
    <tr><th colspan="2" border="1">Figure 6</th></tr>
    <tr><th colspan="2" border="1">Figure 7</th></tr>
    <tr><td>7<strong>a</strong></td><td><a href="code/07_Combinations/0-prepare_main_seurat.ipynb">07_Combinations/0-prepare_main_seurat.ipynb</a></td></tr>
    <tr><td>7<strong>b</strong></td><td><a href="code/07_Combinations/4-state_bias.ipynb">07_Combinations/4-state_bias.ipynb</a></td></tr>
    <tr><td>7<strong>c</strong></td><td><a href="code/07_Combinations/4-state_bias.ipynb">07_Combinations/4-state_bias.ipynb</a></td></tr>
    <tr><td>7<strong>d</strong></td><td><a href="code/07_Combinations/3-pairwise_dosage.ipynb">07_Combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td>7<strong>e</strong></td><td><a href="code/07_Combinations/3-pairwise_dosage.ipynb">07_Combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td>7<strong>e</strong></td><td><a href="code/07_Combinations/5-unique_state_enrichment.ipynb">07_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 2</th></tr>
    <tr><td>ED2<strong>a</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td>ED2<strong>b</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td>ED2<strong>c</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td>ED2<strong>d</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td>ED2<strong>e</strong></td><td>?</td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 4</th></tr>
    <tr><td>ED4<strong>a</strong></td><td>?</td></tr>
    <tr><td>ED4<strong>b</strong></td><td>?</td></tr>
    <tr><td>ED4<strong>c</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td>ED4<strong>d</strong></td><td>?</td></tr>
    <tr><td>ED4<strong>e</strong></td><td><a href="code/ZZ_Power/power.ipynb">ZZ_Power/power.ipynb</a></td></tr>
    <tr><td>ED4<strong>f</strong></td><td><a href="code/ZZ_Power/power.ipynb">ZZ_Power/power.ipynb</a></td></tr>
    <tr><td>ED4<strong>g</strong></td><td><a href="code/ZZ_Power/power.ipynb">ZZ_Power/power.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 9</th></tr>
    <tr><td>ED9<strong>a</strong></td><td><a href="code/07_Combinations/5-unique_state_enrichment.ipynb">07_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><td>ED9<strong>b</strong></td><td><a href="code/07_Combinations/3-pairwise_dosage.ipynb">07_Combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td>ED9<strong>e</strong></td><td>?</td></tr>
    <tr><td>ED9<strong>f</strong></td><td><a href="code/07_Combinations/5-unique_state_enrichment.ipynb">07_Combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Supplementary Figure 3</th></tr>
    <tr><td>S3<strong>a</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>b</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>c</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>d</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>e</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>f</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>g</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td>S3<strong>h</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td>S3<strong>i</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td>S3<strong>j</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td>S3<strong>k</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
  </tbody>
</table>
