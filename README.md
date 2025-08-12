# scTF-seq
Material and source code for the scTF-seq manuscript

## 0. Raw sequencing data
All .fastq files were deposited on ArrayExpress under the accession number [E-MTAB-13010](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13010). It comprises nine 10x libraries (exp5-13), each of them containing different overexpressed TFs. Information about which TF is associated with which library can be found in the [metadata files](metadata/).<br/>
Of note, the "mCherry-BCXX" barcodes were used internally to design the cell lineage studies (adipocyte, myocytes, and chondrocytes differentiation). This information can be found in the [C3H10_10X_Metadata.xlsx](metadata/C3H10_10X_Metadata.xlsx) file.


## 1. Each figure panel can be reproduced using the following scripts/notebooks

<table border="1" cellspacing="0" cellpadding="5">
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
    <tr><td><strong>3a</strong></td><td>?</td></tr>
    <tr><td><strong>3b</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>3c</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>3d</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>3e</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>3e</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>3f</strong></td><td>?</td></tr>
    <tr><td><strong>3g</strong></td><td>?</td></tr>
    <tr><th colspan="2" border="1">Figure 5</th></tr>
    <tr><th colspan="2" border="1">Figure 4</th></tr>
    <tr><th colspan="2" border="1">Figure 6</th></tr>
    <tr><th colspan="2" border="1">Figure 7</th></tr>
    <tr><td><strong>7a</strong></td><td><a href="code/ZZ_combinations/0-prepare_main_seurat.ipynb">ZZ_combinations/0-prepare_main_seurat.ipynb</a></td></tr>
    <tr><td><strong>7b</strong></td><td><a href="code/ZZ_combinations/4-state_bias.ipynb">ZZ_combinations/4-state_bias.ipynb</a></td></tr>
    <tr><td><strong>7c</strong></td><td><a href="code/ZZ_combinations/4-state_bias.ipynb">ZZ_combinations/4-state_bias.ipynb</a></td></tr>
    <tr><td><strong>7d</strong></td><td><a href="code/ZZ_combinations/3-pairwise_dosage.ipynb">ZZ_combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td><strong>7e</strong></td><td><a href="code/ZZ_combinations/3-pairwise_dosage.ipynb">ZZ_combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td><strong>7e</strong></td><td><a href="code/ZZ_combinations/5-unique_state_enrichment.ipynb">ZZ_combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 2</th></tr>
    <tr><td><strong>ED2a</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td><strong>ED2b</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td><strong>ED2c</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td><strong>ED2d</strong></td><td><a href="code/ZZ_Dose_comparison/dose_comparison.py">ZZ_Dose_comparison/dose_comparison.py</a></td></tr>
    <tr><td><strong>ED2e</strong></td><td>?</td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 4</th></tr>
    <tr><td><strong>ED4a</strong></td><td>?</td></tr>
    <tr><td><strong>ED4b</strong></td><td>?</td></tr>
    <tr><td><strong>ED4c</strong></td><td><a href="code/ZZ_Endogenous/2-create_figures.ipynb">ZZ_Endogenous/2-create_figures.ipynb</a></td></tr>
    <tr><td><strong>ED4d</strong></td><td>?</td></tr>
    <tr><td><strong>ED4e</strong></td><td><a href="code/ZZ_power/power.ipynb">ZZ_power/power.ipynb</a></td></tr>
    <tr><td><strong>ED4f</strong></td><td><a href="code/ZZ_power/power.ipynb">ZZ_power/power.ipynb</a></td></tr>
    <tr><td><strong>ED4g</strong></td><td><a href="code/ZZ_power/power.ipynb">ZZ_power/power.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Extended Data Figure 9</th></tr>
    <tr><td><strong>ED9a</strong></td><td><a href="code/ZZ_combinations/5-unique_state_enrichment.ipynb">ZZ_combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><td><strong>ED9b</strong></td><td><a href="code/ZZ_combinations/3-pairwise_dosage.ipynb">ZZ_combinations/3-pairwise_dosage.ipynb</a></td></tr>
    <tr><td><strong>ED9e</strong></td><td>?</td></tr>
    <tr><td><strong>ED9f</strong></td><td><a href="code/ZZ_combinations/5-unique_state_enrichment.ipynb">ZZ_combinations/5-unique_state_enrichment.ipynb</a></td></tr>
    <tr><th colspan="2" border="1">Supplementary Figure 3</th></tr>
    <tr><td><strong>S3a</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3b</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3c</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3d</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3e</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3f</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3g</strong></td><td><a href="code/ZZ_Reproducibility/reproducibility.ipynb">ZZ_Reproducibility/reproducibility.ipynb</a></td></tr>
    <tr><td><strong>S3h</strong></td><td><a href="code/ZZ_Reproducibility/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td><strong>S3i</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td><strong>S3j</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
    <tr><td><strong>S3k</strong></td><td><a href="code/ZZ_Scaling/scaling.ipynb">ZZ_Scaling/scaling.ipynb</a></td></tr>
  </tbody>
</table>