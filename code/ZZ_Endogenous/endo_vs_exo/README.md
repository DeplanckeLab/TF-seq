# Compare endogenous and exogenous overexpression

It is often asserted that overexpression leads to very high "non-physiological" levels of gene expression.
In this notebook we wanted to test this by comparing levels of gene expression in public single-cell datasets
with the levels of overexpression in the scTF-seq dataset. Of particular interest is the level of exogenous expression at which you
observe a significant change in gene expression.
Note the obvious limitations of this approach.
1. We may underestimate physiological levels if:
- Some states of interest may not be present in the public atlas
- Clusters are defined too broad in the public atlas
2. We may underestimate exogenous levels if:
- The overexpression leads to mRNAs that are more easily translated


```python
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import tqdm.auto as tqdm

import polyptich as pp

pp.setup_ipython()

import scanpy as sc

import pickle
```


```python
data_folder = pp.paths.get_data()
```


```python
obs_alltfs = (
    pd.read_csv(data_folder / "obs.csv").rename(columns={"Unnamed: 0": "cell"}).set_index("cell")
)
mtx = scipy.io.mmread("matrix.mtx").T.tocsr()
obs_alltfs["ix"] = range(obs_alltfs.shape[0])

var_alltfs = pd.read_csv(data_folder / "var.csv", index_col=0)
```


```python
obs_alltfs["activity_to_D0"] = (
    pd.read_csv(
        data_folder / "df_allG1Cells_PhaseCorrected_allTFs_D0regressed10pc_50pc_integrated.csv"
    )
    .set_index("cell")["Activity_TF_to_D0"]
    .reindex(obs_alltfs.index)
)
```

Define a set of housekeeping genes to compare against.
We explicitely chose a set of genes from different essential functions such as cytoskeleton, ribosomal proteins, glycolysis, ubiquitinilation, respiration, and translation.


```python
# housekeeping was defined in the 1-download_census script
housekeeping = pickle.load(open("housekeeping.pkl", "rb"))
tfs = obs_alltfs["TF"].unique().tolist()
```


```python
from plot_endogenous import extract_tf_dataset
```


```python
tfcapacity = pd.read_csv(
    data_folder / "TF_categories_on_potency_capacity_dosealigned.csv", index_col=0
)
tfcapacity.index.name = "TF"
```

# Load data


```python
datasets = pickle.load(open("datasets.pkl", "rb"))

adata = pickle.load(open("adata2.pkl", "rb"))
adata.var.index = adata.var.feature_name

adata.obs["dataset_ix"] = pd.Series(datasets.index, index = datasets["dataset_id"])[adata.obs["dataset_id"]].values

adata.obs["n_counts"] = adata.X.sum(axis=1).A1
adata.obs["norm"] = 1 / (adata[:, housekeeping].X.sum(1).A1 + 1)
adata = adata[adata.obs["norm"] < 0.5]

adata.layers["norm"] = adata.X.copy()
adata.layers["norm"] = (
    adata.layers["norm"].multiply(adata.obs["norm"].values[:, None]).tocsr()
)
```

    /tmp/ipykernel_3886466/1375686858.py:12: ImplicitModificationWarning: Setting element `.layers['norm']` of view, initializing view as actual.
      adata.layers["norm"] = adata.X.copy()


We extract the datasets for Myo_ref and Adipo_ref, to see how enodogenous TF expression changes during induction of myogenesis and adipogenesis in our C3H10T1/2 cells.


```python
adatamyo = extract_tf_dataset(mtx, obs_alltfs, var_alltfs, "Myo_ref", housekeeping=housekeeping)
adatamyo = adatamyo[adatamyo.obs["Phase_corrected"] == "G1"]
adatamyo = adatamyo[adatamyo.obs["TF"] == "Myo_ref"]
adatamyo = adatamyo[adatamyo.obs.index[:200]]

adataadipo = extract_tf_dataset(mtx, obs_alltfs, var_alltfs, "Adipo_ref", housekeeping=housekeeping)
adataadipo = adataadipo[adataadipo.obs["Phase_corrected"] == "G1"]
adataadipo = adataadipo[adataadipo.obs["TF"] == "Adipo_ref"]
adataadipo = adataadipo[adataadipo.obs.index[:200]]
```

## Individual


```python
gene = "Hnf4g"
# gene = "Zeb1"
# gene = "Hoxa9"
# gene = "Myod1"
# gene = "Myog"
# gene = "Pparg"
# gene = "Cebpa"
gene = "Runx2"
tf = gene
adatatf = extract_tf_dataset(mtx, obs_alltfs, var_alltfs, tf, housekeeping=housekeeping)
adatatf = adatatf[adatatf.obs["Phase_corrected"] == "G1"]
adatatf_norm = adatatf.copy()

adatatf_oi = adatatf[adatatf.obs["TF"] == gene]
```


```python
adatatf.obs["target"] = adatatf.obs["activity_to_D0"]
```

    /tmp/ipykernel_3886466/1284872391.py:1: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



```python
from plot_endogenous import plot as plot_endogenous
```


```python
tf = "Pparg"
fig = plot_endogenous(
    mtx, obs_alltfs, var_alltfs, tf, tf, housekeeping, adata, datasets=datasets, adatamyo = adatamyo, adataadipo = adataadipo, n_top = 20, width = 1.5
)
fig.display()
```

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]


    R[write to console]: Loading required package: mgcv
    


    R[write to console]: Loading required package: nlme
    


    R[write to console]: This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.
    



    
![png](output_17_4.png)
    



```python
import pathlib
import tqdm

plot_folder = pathlib.Path("./plots_endogenous")
import shutil

# shutil.rmtree(plot_folder, ignore_errors=True)
plot_folder.mkdir(exist_ok=True)

tfs_oi = [tf for tf in tfs if tf in adata.var.index]
tfs_oi = ["Cebpa", "Pparg", "Runx2", "Myog", "Nfkb1", "Meis2", "Fosl2"]

for gene in tqdm.tqdm(tfs_oi):
    fig = plot_endogenous(
        mtx, obs_alltfs, var_alltfs, gene, gene, housekeeping, adata, datasets=datasets, adatamyo = adatamyo, adataadipo = adataadipo
    )
    fig.savefig(plot_folder / f"{gene}.pdf", display = False, transparent=True)
    fig.savefig(plot_folder / f"{gene}.png", display = True, transparent=True)
```

      0%|                                                                                           | 0/7 [00:00<?, ?it/s]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_2.png)
    


     14%|███████████▊                                                                       | 1/7 [00:01<00:08,  1.44s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_5.png)
    


     29%|███████████████████████▋                                                           | 2/7 [00:02<00:07,  1.48s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_8.png)
    


     43%|███████████████████████████████████▌                                               | 3/7 [00:04<00:06,  1.51s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_11.png)
    


     57%|███████████████████████████████████████████████▍                                   | 4/7 [00:05<00:04,  1.40s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_14.png)
    


     71%|███████████████████████████████████████████████████████████▎                       | 5/7 [00:06<00:02,  1.35s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_17.png)
    


     86%|███████████████████████████████████████████████████████████████████████▏           | 6/7 [00:08<00:01,  1.40s/it]

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_18_20.png)
    


    100%|███████████████████████████████████████████████████████████████████████████████████| 7/7 [00:09<00:00,  1.34s/it]

    100%|███████████████████████████████████████████████████████████████████████████████████| 7/7 [00:09<00:00,  1.39s/it]

    



```python
citations = {}
fig = plot_endogenous(
        mtx, obs_alltfs, var_alltfs, "Cdx2", "Cdx2", housekeeping, adata, datasets=datasets, adatamyo = adatamyo, adataadipo = adataadipo, width = 1.5, n_top = 10, add_citations = False, citations = citations
    )
fig.display()
```

    /srv/data/wouters/projects/tfseq/code/ZZ_Endogenous/endo_vs_exo/plot_endogenous.py:173: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
      adatatf.obs["target"] = adatatf.obs["activity_to_D0"]



    
![png](output_19_1.png)
    



```python
datasets_oi = datasets.loc[citations.keys()]
datasets_oi["footnote_symbol"] = [citations[x] for x in datasets_oi.index]
plotdata = datasets_oi.groupby("collection_doi_label").agg(
    {"footnote_symbol": tuple}
)
plotdata = plotdata.sort_values("footnote_symbol", ascending=False)
plotdata["ix"] = range(plotdata.shape[0])

fig, ax = plt.subplots(figsize = (2, len(plotdata) * 0.15))
for label, dataset in plotdata.iterrows():
    label = label.replace("The Tabula Sapiens Consortium* et al.", "Tabula Sapiens")
    ax.text(
        0,
        dataset["ix"],
        f"$^{{{','.join(dataset['footnote_symbol'])}}}$ {label}",
        ha="left",
        va="center",
        fontsize=7,
    )
ax.set_ylim(-0.5, len(plotdata) - 0.5)
ax.axis("off")
fig.savefig(plot_folder / f"small_footnotes.pdf", transparent=True)
```


    
![png](output_20_0.png)
    


## Global


```python
from plot_endogenous import smooth_spline_fit_se, abbreviate_sentence

tfs_oi = ["Runx2"]
tfs_oi = tfcapacity.query("category != 'low-capacity'").index
tfs_oi = [x for x in tfcapacity.index if x not in ["Pparg1", "Spert"]]
tfs_oi = [x for x in tfs_oi if x in adata.var.index]

tfscores_exogenous = []
for gene in tqdm.tqdm(tfs_oi):
    adatatf = extract_tf_dataset(mtx, obs_alltfs, var_alltfs, gene, housekeeping)
    adatatf = adatatf[adatatf.obs["Phase_corrected"] == "G1"].copy()
    adatatf.obs["target"] = adatatf.obs["activity_to_D0"]

    plotdata = pd.DataFrame(
        {
            "expression": np.log(adatatf.obs["Vector_10X_norm"] + 1),
            "dataset_cell_type": "scTF-seq",
            "target": adatatf.obs["target"],
        }
    )
    plotdata_0 = plotdata.loc[plotdata["expression"] == 0]
    plotdata_non0 = plotdata.loc[plotdata["expression"] != 0]

    extent = plotdata_non0["expression"].max()
    extent_left = -extent * 0.02
    plotdata_smooth = pd.DataFrame(
        {
            "expression": np.linspace(0, extent, 100),
        }
    )
    plotdata_smooth["target"], plotdata_smooth["target_se"] = smooth_spline_fit_se(
        plotdata["expression"], plotdata["target"], plotdata_smooth["expression"]
    ).T

    cutoff = 0.23
    if plotdata_smooth["target"].max() < cutoff:
        target = plotdata_smooth["expression"].loc[
                plotdata_smooth["target"].argmax()
            ]
    else:
        target = plotdata_smooth["expression"][
            plotdata_smooth["target"] > cutoff
        ].iloc[0]


    tfscores_exogenous.append(
        {
            "tf": gene,
            "target_03": target,
            "max_dose": plotdata_smooth["expression"].max(),
            "never_reached": plotdata_smooth["target"].max() < cutoff,
        }
    )

tfscores_exogenous = pd.DataFrame(tfscores_exogenous)
```

      0%|                                                                                         | 0/227 [00:00<?, ?it/s]

      0%|▎                                                                                | 1/227 [00:00<00:36,  6.27it/s]

      1%|█                                                                                | 3/227 [00:00<00:21, 10.38it/s]

      2%|█▊                                                                               | 5/227 [00:00<00:18, 11.89it/s]

      3%|██▍                                                                              | 7/227 [00:00<00:17, 12.50it/s]

      4%|███▏                                                                             | 9/227 [00:00<00:26,  8.30it/s]

      5%|███▉                                                                            | 11/227 [00:01<00:22,  9.47it/s]

      6%|████▌                                                                           | 13/227 [00:01<00:20, 10.61it/s]

      7%|█████▎                                                                          | 15/227 [00:01<00:18, 11.40it/s]

      7%|█████▉                                                                          | 17/227 [00:01<00:17, 12.15it/s]

      8%|██████▋                                                                         | 19/227 [00:01<00:17, 11.94it/s]

      9%|███████▍                                                                        | 21/227 [00:01<00:17, 11.46it/s]

     10%|████████                                                                        | 23/227 [00:02<00:16, 12.07it/s]

     11%|████████▊                                                                       | 25/227 [00:02<00:17, 11.51it/s]

     12%|█████████▌                                                                      | 27/227 [00:02<00:18, 10.92it/s]

     13%|██████████▏                                                                     | 29/227 [00:02<00:16, 11.65it/s]

     14%|██████████▉                                                                     | 31/227 [00:02<00:15, 12.29it/s]

     15%|███████████▋                                                                    | 33/227 [00:02<00:15, 12.80it/s]

     15%|████████████▎                                                                   | 35/227 [00:03<00:14, 13.17it/s]

     16%|█████████████                                                                   | 37/227 [00:03<00:14, 13.44it/s]

     17%|█████████████▋                                                                  | 39/227 [00:03<00:13, 13.61it/s]

     18%|██████████████▍                                                                 | 41/227 [00:03<00:13, 13.79it/s]

     19%|███████████████▏                                                                | 43/227 [00:03<00:13, 13.79it/s]

     20%|███████████████▊                                                                | 45/227 [00:03<00:13, 13.86it/s]

     21%|████████████████▌                                                               | 47/227 [00:03<00:12, 14.06it/s]

     22%|█████████████████▎                                                              | 49/227 [00:04<00:12, 14.14it/s]

     22%|█████████████████▉                                                              | 51/227 [00:04<00:13, 13.50it/s]

     23%|██████████████████▋                                                             | 53/227 [00:04<00:13, 12.71it/s]

     24%|███████████████████▍                                                            | 55/227 [00:04<00:14, 11.75it/s]

     25%|████████████████████                                                            | 57/227 [00:04<00:16, 10.15it/s]

     26%|████████████████████▊                                                           | 59/227 [00:05<00:16, 10.04it/s]

     27%|█████████████████████▍                                                          | 61/227 [00:05<00:16, 10.34it/s]

     28%|██████████████████████▏                                                         | 63/227 [00:05<00:15, 10.55it/s]

     29%|██████████████████████▉                                                         | 65/227 [00:05<00:15, 10.29it/s]

     30%|███████████████████████▌                                                        | 67/227 [00:05<00:15, 10.29it/s]

     30%|████████████████████████▎                                                       | 69/227 [00:05<00:15, 10.53it/s]

     31%|█████████████████████████                                                       | 71/227 [00:06<00:14, 10.60it/s]

     32%|█████████████████████████▋                                                      | 73/227 [00:06<00:14, 10.75it/s]

     33%|██████████████████████████▍                                                     | 75/227 [00:06<00:16,  9.10it/s]

     34%|███████████████████████████▏                                                    | 77/227 [00:06<00:15,  9.56it/s]

     35%|███████████████████████████▊                                                    | 79/227 [00:07<00:14,  9.97it/s]

     36%|████████████████████████████▌                                                   | 81/227 [00:07<00:14, 10.33it/s]

     37%|█████████████████████████████▎                                                  | 83/227 [00:07<00:13, 10.34it/s]

     37%|█████████████████████████████▉                                                  | 85/227 [00:07<00:13, 10.54it/s]

     38%|██████████████████████████████▋                                                 | 87/227 [00:07<00:14,  9.37it/s]

     39%|███████████████████████████████                                                 | 88/227 [00:07<00:14,  9.39it/s]

     40%|███████████████████████████████▋                                                | 90/227 [00:08<00:14,  9.54it/s]

     41%|████████████████████████████████▍                                               | 92/227 [00:08<00:13,  9.88it/s]

     41%|█████████████████████████████████▏                                              | 94/227 [00:08<00:13, 10.19it/s]

     42%|█████████████████████████████████▊                                              | 96/227 [00:08<00:12, 10.52it/s]

     43%|██████████████████████████████████▌                                             | 98/227 [00:08<00:12, 10.59it/s]

     44%|██████████████████████████████████▊                                            | 100/227 [00:09<00:11, 10.75it/s]

     45%|███████████████████████████████████▍                                           | 102/227 [00:09<00:11, 10.84it/s]

     46%|████████████████████████████████████▏                                          | 104/227 [00:09<00:11, 10.85it/s]

     47%|████████████████████████████████████▉                                          | 106/227 [00:09<00:11, 10.85it/s]

     48%|█████████████████████████████████████▌                                         | 108/227 [00:09<00:10, 10.94it/s]

     48%|██████████████████████████████████████▎                                        | 110/227 [00:09<00:10, 10.92it/s]

     49%|██████████████████████████████████████▉                                        | 112/227 [00:10<00:11, 10.13it/s]

     50%|███████████████████████████████████████▋                                       | 114/227 [00:10<00:10, 10.30it/s]

     51%|████████████████████████████████████████▎                                      | 116/227 [00:10<00:10, 10.14it/s]

     52%|█████████████████████████████████████████                                      | 118/227 [00:10<00:11,  9.35it/s]

     52%|█████████████████████████████████████████▍                                     | 119/227 [00:11<00:16,  6.55it/s]

     53%|█████████████████████████████████████████▊                                     | 120/227 [00:11<00:15,  7.01it/s]

     54%|██████████████████████████████████████████▍                                    | 122/227 [00:11<00:12,  8.08it/s]

     55%|███████████████████████████████████████████▏                                   | 124/227 [00:11<00:11,  9.02it/s]

     56%|███████████████████████████████████████████▊                                   | 126/227 [00:11<00:11,  9.18it/s]

     56%|████████████████████████████████████████████▌                                  | 128/227 [00:12<00:10,  9.65it/s]

     57%|█████████████████████████████████████████████▏                                 | 130/227 [00:12<00:12,  7.98it/s]

     58%|█████████████████████████████████████████████▉                                 | 132/227 [00:12<00:10,  8.82it/s]

     59%|██████████████████████████████████████████████▋                                | 134/227 [00:12<00:09,  9.51it/s]

     60%|███████████████████████████████████████████████▎                               | 136/227 [00:12<00:09, 10.03it/s]

     61%|████████████████████████████████████████████████                               | 138/227 [00:13<00:08, 10.46it/s]

     62%|████████████████████████████████████████████████▋                              | 140/227 [00:13<00:08, 10.27it/s]

     63%|█████████████████████████████████████████████████▍                             | 142/227 [00:13<00:07, 10.70it/s]

     63%|██████████████████████████████████████████████████                             | 144/227 [00:13<00:08,  9.73it/s]

     64%|██████████████████████████████████████████████████▊                            | 146/227 [00:13<00:08, 10.11it/s]

     65%|███████████████████████████████████████████████████▌                           | 148/227 [00:14<00:07, 10.51it/s]

     66%|████████████████████████████████████████████████████▏                          | 150/227 [00:14<00:10,  7.28it/s]

     67%|████████████████████████████████████████████████████▉                          | 152/227 [00:14<00:09,  8.21it/s]

     68%|█████████████████████████████████████████████████████▌                         | 154/227 [00:14<00:08,  8.96it/s]

     69%|██████████████████████████████████████████████████████▎                        | 156/227 [00:15<00:07,  9.52it/s]

     70%|██████████████████████████████████████████████████████▉                        | 158/227 [00:15<00:06, 10.05it/s]

     70%|███████████████████████████████████████████████████████▋                       | 160/227 [00:15<00:06, 10.45it/s]

     71%|████████████████████████████████████████████████████████▍                      | 162/227 [00:15<00:07,  8.25it/s]

     72%|█████████████████████████████████████████████████████████                      | 164/227 [00:15<00:07,  8.96it/s]

     73%|█████████████████████████████████████████████████████████▊                     | 166/227 [00:16<00:06,  9.60it/s]

     74%|██████████████████████████████████████████████████████████▍                    | 168/227 [00:16<00:05,  9.97it/s]

     75%|███████████████████████████████████████████████████████████▏                   | 170/227 [00:16<00:06,  9.22it/s]

     75%|███████████████████████████████████████████████████████████▌                   | 171/227 [00:16<00:06,  9.25it/s]

     76%|████████████████████████████████████████████████████████████▏                  | 173/227 [00:16<00:05,  9.95it/s]

     77%|████████████████████████████████████████████████████████████▉                  | 175/227 [00:17<00:05,  9.74it/s]

     78%|█████████████████████████████████████████████████████████████▌                 | 177/227 [00:17<00:04, 10.24it/s]

     79%|██████████████████████████████████████████████████████████████▎                | 179/227 [00:17<00:04, 10.45it/s]

     80%|██████████████████████████████████████████████████████████████▉                | 181/227 [00:17<00:04, 10.80it/s]

     81%|███████████████████████████████████████████████████████████████▋               | 183/227 [00:17<00:04, 10.98it/s]

     81%|████████████████████████████████████████████████████████████████▍              | 185/227 [00:17<00:03, 11.10it/s]

     82%|█████████████████████████████████████████████████████████████████              | 187/227 [00:18<00:04,  8.61it/s]

     83%|█████████████████████████████████████████████████████████████████▍             | 188/227 [00:18<00:04,  8.80it/s]

     84%|██████████████████████████████████████████████████████████████████             | 190/227 [00:18<00:03,  9.43it/s]

     85%|██████████████████████████████████████████████████████████████████▊            | 192/227 [00:18<00:03,  9.86it/s]

     85%|███████████████████████████████████████████████████████████████████▌           | 194/227 [00:18<00:03,  9.95it/s]

     86%|████████████████████████████████████████████████████████████████████▏          | 196/227 [00:19<00:03, 10.26it/s]

     87%|████████████████████████████████████████████████████████████████████▉          | 198/227 [00:19<00:02, 10.41it/s]

     88%|█████████████████████████████████████████████████████████████████████▌         | 200/227 [00:19<00:02, 10.61it/s]

     89%|██████████████████████████████████████████████████████████████████████▎        | 202/227 [00:19<00:02, 10.66it/s]

     90%|██████████████████████████████████████████████████████████████████████▉        | 204/227 [00:19<00:02, 10.53it/s]

     91%|███████████████████████████████████████████████████████████████████████▋       | 206/227 [00:20<00:02,  7.66it/s]

     92%|████████████████████████████████████████████████████████████████████████▍      | 208/227 [00:20<00:02,  8.51it/s]

     93%|█████████████████████████████████████████████████████████████████████████      | 210/227 [00:20<00:01,  9.27it/s]

     93%|█████████████████████████████████████████████████████████████████████████▊     | 212/227 [00:20<00:01,  9.79it/s]

     94%|██████████████████████████████████████████████████████████████████████████▍    | 214/227 [00:20<00:01, 10.33it/s]

     95%|███████████████████████████████████████████████████████████████████████████▏   | 216/227 [00:21<00:01, 10.73it/s]

     96%|███████████████████████████████████████████████████████████████████████████▊   | 218/227 [00:21<00:00, 10.89it/s]

     97%|████████████████████████████████████████████████████████████████████████████▌  | 220/227 [00:21<00:00, 11.04it/s]

     98%|█████████████████████████████████████████████████████████████████████████████▎ | 222/227 [00:21<00:00, 11.15it/s]

     99%|█████████████████████████████████████████████████████████████████████████████▉ | 224/227 [00:21<00:00, 11.18it/s]

    100%|██████████████████████████████████████████████████████████████████████████████▋| 226/227 [00:22<00:00, 11.03it/s]

    100%|███████████████████████████████████████████████████████████████████████████████| 227/227 [00:22<00:00, 10.25it/s]

    



```python
tfscores_endogenous = []

for gene in tqdm.tqdm(tfs_oi):
    plotdata = pd.DataFrame(
        {
            "expression": np.log1p(sc.get.obs_df(adata, gene, layer="norm").values[:, 0]),
            "dataset_ix": adata.obs["dataset_ix"],
            "cell_type": adata.obs["cell_type"],
        }
    )
    plotdata["dataset_cell_type"] = (
        plotdata["cell_type"].astype(str).apply(abbreviate_sentence, max_length=50)
        + " "
        + plotdata["dataset_ix"].astype(str)
    )

    dataset_cell_type_oi = (
        plotdata.groupby("dataset_cell_type")["expression"]
        .mean()
        .sort_values(ascending=False)
        .index[0]
    )
    plotdata_oi = plotdata.query("dataset_cell_type == @dataset_cell_type_oi")

    tfscores_endogenous.append(
        {
            "q90": plotdata_oi["expression"].quantile(0.90),
            "q95": plotdata_oi["expression"].quantile(0.94),
            "q99": plotdata_oi["expression"].quantile(0.99),
            "q10": plotdata_oi["expression"].quantile(0.10),
            "q05": plotdata_oi["expression"].quantile(0.05),
            "mean": plotdata_oi["expression"].mean(),
            "tf": gene,
            "med": plotdata_oi["expression"].median(),
            "cell_type": ' '.join(dataset_cell_type_oi.split(" ")[:-1]),
        }
    )
tfscores_endogenous = pd.DataFrame(tfscores_endogenous)
```

      0%|                                                                                         | 0/227 [00:00<?, ?it/s]

      0%|▎                                                                                | 1/227 [00:00<01:04,  3.51it/s]

      1%|▋                                                                                | 2/227 [00:00<01:04,  3.48it/s]

      1%|█                                                                                | 3/227 [00:00<01:02,  3.58it/s]

      2%|█▍                                                                               | 4/227 [00:01<01:01,  3.61it/s]

      2%|█▊                                                                               | 5/227 [00:01<01:01,  3.64it/s]

      3%|██▏                                                                              | 6/227 [00:01<01:00,  3.63it/s]

      3%|██▍                                                                              | 7/227 [00:01<01:00,  3.66it/s]

      4%|██▊                                                                              | 8/227 [00:02<00:59,  3.65it/s]

      4%|███▏                                                                             | 9/227 [00:02<00:59,  3.65it/s]

      4%|███▌                                                                            | 10/227 [00:02<00:59,  3.66it/s]

      5%|███▉                                                                            | 11/227 [00:03<00:58,  3.67it/s]

      5%|████▏                                                                           | 12/227 [00:03<01:00,  3.56it/s]

      6%|████▌                                                                           | 13/227 [00:03<00:59,  3.59it/s]

      6%|████▉                                                                           | 14/227 [00:03<00:58,  3.62it/s]

      7%|█████▎                                                                          | 15/227 [00:04<00:58,  3.64it/s]

      7%|█████▋                                                                          | 16/227 [00:04<00:57,  3.64it/s]

      7%|█████▉                                                                          | 17/227 [00:04<00:57,  3.63it/s]

      8%|██████▎                                                                         | 18/227 [00:04<00:57,  3.65it/s]

      8%|██████▋                                                                         | 19/227 [00:05<00:56,  3.66it/s]

      9%|███████                                                                         | 20/227 [00:05<00:57,  3.59it/s]

      9%|███████▍                                                                        | 21/227 [00:05<00:57,  3.58it/s]

     10%|███████▊                                                                        | 22/227 [00:06<00:56,  3.61it/s]

     10%|████████                                                                        | 23/227 [00:06<00:56,  3.63it/s]

     11%|████████▍                                                                       | 24/227 [00:06<00:55,  3.64it/s]

     11%|████████▊                                                                       | 25/227 [00:06<00:55,  3.65it/s]

     11%|█████████▏                                                                      | 26/227 [00:07<00:54,  3.66it/s]

     12%|█████████▌                                                                      | 27/227 [00:07<00:54,  3.66it/s]

     12%|█████████▊                                                                      | 28/227 [00:07<00:54,  3.63it/s]

     13%|██████████▏                                                                     | 29/227 [00:07<00:54,  3.63it/s]

     13%|██████████▌                                                                     | 30/227 [00:08<00:54,  3.64it/s]

     14%|██████████▉                                                                     | 31/227 [00:08<00:53,  3.65it/s]

     14%|███████████▎                                                                    | 32/227 [00:08<00:53,  3.66it/s]

     15%|███████████▋                                                                    | 33/227 [00:09<00:53,  3.65it/s]

     15%|███████████▉                                                                    | 34/227 [00:09<00:52,  3.66it/s]

     15%|████████████▎                                                                   | 35/227 [00:09<00:52,  3.67it/s]

     16%|████████████▋                                                                   | 36/227 [00:09<00:52,  3.66it/s]

     16%|█████████████                                                                   | 37/227 [00:10<00:52,  3.63it/s]

     17%|█████████████▍                                                                  | 38/227 [00:10<00:51,  3.64it/s]

     17%|█████████████▋                                                                  | 39/227 [00:10<00:51,  3.65it/s]

     18%|██████████████                                                                  | 40/227 [00:11<00:50,  3.67it/s]

     18%|██████████████▍                                                                 | 41/227 [00:11<00:51,  3.64it/s]

     19%|██████████████▊                                                                 | 42/227 [00:11<00:50,  3.65it/s]

     19%|███████████████▏                                                                | 43/227 [00:11<00:50,  3.65it/s]

     19%|███████████████▌                                                                | 44/227 [00:12<00:50,  3.65it/s]

     20%|███████████████▊                                                                | 45/227 [00:12<00:49,  3.66it/s]

     20%|████████████████▏                                                               | 46/227 [00:12<00:49,  3.66it/s]

     21%|████████████████▌                                                               | 47/227 [00:12<00:49,  3.65it/s]

     21%|████████████████▉                                                               | 48/227 [00:13<00:49,  3.65it/s]

     22%|█████████████████▎                                                              | 49/227 [00:13<00:48,  3.65it/s]

     22%|█████████████████▌                                                              | 50/227 [00:13<00:48,  3.65it/s]

     22%|█████████████████▉                                                              | 51/227 [00:14<00:48,  3.65it/s]

     23%|██████████████████▎                                                             | 52/227 [00:14<00:48,  3.64it/s]

     23%|██████████████████▋                                                             | 53/227 [00:14<00:47,  3.65it/s]

     24%|███████████████████                                                             | 54/227 [00:14<00:48,  3.59it/s]

     24%|███████████████████▍                                                            | 55/227 [00:15<00:47,  3.61it/s]

     25%|███████████████████▋                                                            | 56/227 [00:15<00:47,  3.63it/s]

     25%|████████████████████                                                            | 57/227 [00:15<00:46,  3.65it/s]

     26%|████████████████████▍                                                           | 58/227 [00:15<00:46,  3.66it/s]

     26%|████████████████████▊                                                           | 59/227 [00:16<00:45,  3.67it/s]

     26%|█████████████████████▏                                                          | 60/227 [00:16<00:45,  3.67it/s]

     27%|█████████████████████▍                                                          | 61/227 [00:16<00:45,  3.67it/s]

     27%|█████████████████████▊                                                          | 62/227 [00:17<00:45,  3.66it/s]

     28%|██████████████████████▏                                                         | 63/227 [00:17<00:44,  3.67it/s]

     28%|██████████████████████▌                                                         | 64/227 [00:17<00:44,  3.67it/s]

     29%|██████████████████████▉                                                         | 65/227 [00:17<00:44,  3.65it/s]

     29%|███████████████████████▎                                                        | 66/227 [00:18<00:44,  3.63it/s]

     30%|███████████████████████▌                                                        | 67/227 [00:18<00:43,  3.64it/s]

     30%|███████████████████████▉                                                        | 68/227 [00:18<00:43,  3.64it/s]

     30%|████████████████████████▎                                                       | 69/227 [00:18<00:43,  3.64it/s]

     31%|████████████████████████▋                                                       | 70/227 [00:19<00:43,  3.65it/s]

     31%|█████████████████████████                                                       | 71/227 [00:19<00:42,  3.63it/s]

     32%|█████████████████████████▎                                                      | 72/227 [00:19<00:42,  3.64it/s]

     32%|█████████████████████████▋                                                      | 73/227 [00:20<00:42,  3.65it/s]

     33%|██████████████████████████                                                      | 74/227 [00:20<00:41,  3.66it/s]

     33%|██████████████████████████▍                                                     | 75/227 [00:20<00:41,  3.67it/s]

     33%|██████████████████████████▊                                                     | 76/227 [00:20<00:41,  3.65it/s]

     34%|███████████████████████████▏                                                    | 77/227 [00:21<00:41,  3.66it/s]

     34%|███████████████████████████▍                                                    | 78/227 [00:21<00:40,  3.67it/s]

     35%|███████████████████████████▊                                                    | 79/227 [00:21<00:40,  3.65it/s]

     35%|████████████████████████████▏                                                   | 80/227 [00:21<00:40,  3.64it/s]

     36%|████████████████████████████▌                                                   | 81/227 [00:22<00:40,  3.65it/s]

     36%|████████████████████████████▉                                                   | 82/227 [00:22<00:39,  3.65it/s]

     37%|█████████████████████████████▎                                                  | 83/227 [00:22<00:39,  3.65it/s]

     37%|█████████████████████████████▌                                                  | 84/227 [00:23<00:39,  3.66it/s]

     37%|█████████████████████████████▉                                                  | 85/227 [00:23<00:38,  3.66it/s]

     38%|██████████████████████████████▎                                                 | 86/227 [00:23<00:38,  3.65it/s]

     38%|██████████████████████████████▋                                                 | 87/227 [00:23<00:38,  3.63it/s]

     39%|███████████████████████████████                                                 | 88/227 [00:24<00:38,  3.64it/s]

     39%|███████████████████████████████▎                                                | 89/227 [00:24<00:38,  3.62it/s]

     40%|███████████████████████████████▋                                                | 90/227 [00:24<00:37,  3.62it/s]

     40%|████████████████████████████████                                                | 91/227 [00:24<00:37,  3.63it/s]

     41%|████████████████████████████████▍                                               | 92/227 [00:25<00:37,  3.64it/s]

     41%|████████████████████████████████▊                                               | 93/227 [00:25<00:36,  3.64it/s]

     41%|█████████████████████████████████▏                                              | 94/227 [00:25<00:36,  3.64it/s]

     42%|█████████████████████████████████▍                                              | 95/227 [00:26<00:36,  3.64it/s]

     42%|█████████████████████████████████▊                                              | 96/227 [00:26<00:36,  3.64it/s]

     43%|██████████████████████████████████▏                                             | 97/227 [00:26<00:35,  3.63it/s]

     43%|██████████████████████████████████▌                                             | 98/227 [00:26<00:35,  3.64it/s]

     44%|██████████████████████████████████▉                                             | 99/227 [00:27<00:35,  3.65it/s]

     44%|██████████████████████████████████▊                                            | 100/227 [00:27<00:34,  3.65it/s]

     44%|███████████████████████████████████▏                                           | 101/227 [00:27<00:34,  3.65it/s]

     45%|███████████████████████████████████▍                                           | 102/227 [00:28<00:34,  3.65it/s]

     45%|███████████████████████████████████▊                                           | 103/227 [00:28<00:34,  3.63it/s]

     46%|████████████████████████████████████▏                                          | 104/227 [00:28<00:34,  3.62it/s]

     46%|████████████████████████████████████▌                                          | 105/227 [00:28<00:33,  3.60it/s]

     47%|████████████████████████████████████▉                                          | 106/227 [00:29<00:33,  3.59it/s]

     47%|█████████████████████████████████████▏                                         | 107/227 [00:29<00:33,  3.59it/s]

     48%|█████████████████████████████████████▌                                         | 108/227 [00:29<00:33,  3.61it/s]

     48%|█████████████████████████████████████▉                                         | 109/227 [00:29<00:32,  3.62it/s]

     48%|██████████████████████████████████████▎                                        | 110/227 [00:30<00:32,  3.63it/s]

     49%|██████████████████████████████████████▋                                        | 111/227 [00:30<00:31,  3.63it/s]

     49%|██████████████████████████████████████▉                                        | 112/227 [00:30<00:31,  3.61it/s]

     50%|███████████████████████████████████████▎                                       | 113/227 [00:31<00:31,  3.62it/s]

     50%|███████████████████████████████████████▋                                       | 114/227 [00:31<00:31,  3.63it/s]

     51%|████████████████████████████████████████                                       | 115/227 [00:31<00:30,  3.63it/s]

     51%|████████████████████████████████████████▎                                      | 116/227 [00:31<00:30,  3.64it/s]

     52%|████████████████████████████████████████▋                                      | 117/227 [00:32<00:30,  3.62it/s]

     52%|█████████████████████████████████████████                                      | 118/227 [00:32<00:30,  3.62it/s]

     52%|█████████████████████████████████████████▍                                     | 119/227 [00:32<00:29,  3.62it/s]

     53%|█████████████████████████████████████████▊                                     | 120/227 [00:32<00:29,  3.61it/s]

     53%|██████████████████████████████████████████                                     | 121/227 [00:33<00:29,  3.62it/s]

     54%|██████████████████████████████████████████▍                                    | 122/227 [00:33<00:29,  3.62it/s]

     54%|██████████████████████████████████████████▊                                    | 123/227 [00:33<00:28,  3.61it/s]

     55%|███████████████████████████████████████████▏                                   | 124/227 [00:34<00:28,  3.62it/s]

     55%|███████████████████████████████████████████▌                                   | 125/227 [00:34<00:28,  3.61it/s]

     56%|███████████████████████████████████████████▊                                   | 126/227 [00:34<00:27,  3.62it/s]

     56%|████████████████████████████████████████████▏                                  | 127/227 [00:34<00:27,  3.61it/s]

     56%|████████████████████████████████████████████▌                                  | 128/227 [00:35<00:27,  3.60it/s]

     57%|████████████████████████████████████████████▉                                  | 129/227 [00:35<00:27,  3.60it/s]

     57%|█████████████████████████████████████████████▏                                 | 130/227 [00:35<00:26,  3.61it/s]

     58%|█████████████████████████████████████████████▌                                 | 131/227 [00:36<00:26,  3.59it/s]

     58%|█████████████████████████████████████████████▉                                 | 132/227 [00:36<00:26,  3.60it/s]

     59%|██████████████████████████████████████████████▎                                | 133/227 [00:36<00:26,  3.60it/s]

     59%|██████████████████████████████████████████████▋                                | 134/227 [00:36<00:25,  3.61it/s]

     59%|██████████████████████████████████████████████▉                                | 135/227 [00:37<00:25,  3.61it/s]

     60%|███████████████████████████████████████████████▎                               | 136/227 [00:37<00:25,  3.61it/s]

     60%|███████████████████████████████████████████████▋                               | 137/227 [00:37<00:24,  3.60it/s]

     61%|████████████████████████████████████████████████                               | 138/227 [00:37<00:24,  3.61it/s]

     61%|████████████████████████████████████████████████▎                              | 139/227 [00:38<00:24,  3.60it/s]

     62%|████████████████████████████████████████████████▋                              | 140/227 [00:38<00:24,  3.61it/s]

     62%|█████████████████████████████████████████████████                              | 141/227 [00:38<00:23,  3.61it/s]

     63%|█████████████████████████████████████████████████▍                             | 142/227 [00:39<00:23,  3.60it/s]

     63%|█████████████████████████████████████████████████▊                             | 143/227 [00:39<00:23,  3.60it/s]

     63%|██████████████████████████████████████████████████                             | 144/227 [00:39<00:23,  3.60it/s]

     64%|██████████████████████████████████████████████████▍                            | 145/227 [00:39<00:22,  3.59it/s]

     64%|██████████████████████████████████████████████████▊                            | 146/227 [00:40<00:22,  3.60it/s]

     65%|███████████████████████████████████████████████████▏                           | 147/227 [00:40<00:22,  3.60it/s]

     65%|███████████████████████████████████████████████████▌                           | 148/227 [00:40<00:22,  3.58it/s]

     66%|███████████████████████████████████████████████████▊                           | 149/227 [00:41<00:21,  3.58it/s]

     66%|████████████████████████████████████████████████████▏                          | 150/227 [00:41<00:21,  3.59it/s]

     67%|████████████████████████████████████████████████████▌                          | 151/227 [00:41<00:21,  3.59it/s]

     67%|████████████████████████████████████████████████████▉                          | 152/227 [00:41<00:20,  3.59it/s]

     67%|█████████████████████████████████████████████████████▏                         | 153/227 [00:42<00:20,  3.59it/s]

     68%|█████████████████████████████████████████████████████▌                         | 154/227 [00:42<00:20,  3.60it/s]

     68%|█████████████████████████████████████████████████████▉                         | 155/227 [00:42<00:20,  3.60it/s]

     69%|██████████████████████████████████████████████████████▎                        | 156/227 [00:42<00:19,  3.59it/s]

     69%|██████████████████████████████████████████████████████▋                        | 157/227 [00:43<00:19,  3.61it/s]

     70%|██████████████████████████████████████████████████████▉                        | 158/227 [00:43<00:19,  3.62it/s]

     70%|███████████████████████████████████████████████████████▎                       | 159/227 [00:43<00:18,  3.62it/s]

     70%|███████████████████████████████████████████████████████▋                       | 160/227 [00:44<00:18,  3.61it/s]

     71%|████████████████████████████████████████████████████████                       | 161/227 [00:44<00:18,  3.60it/s]

     71%|████████████████████████████████████████████████████████▍                      | 162/227 [00:44<00:18,  3.61it/s]

     72%|████████████████████████████████████████████████████████▋                      | 163/227 [00:44<00:17,  3.62it/s]

     72%|█████████████████████████████████████████████████████████                      | 164/227 [00:45<00:17,  3.60it/s]

     73%|█████████████████████████████████████████████████████████▍                     | 165/227 [00:45<00:17,  3.60it/s]

     73%|█████████████████████████████████████████████████████████▊                     | 166/227 [00:45<00:16,  3.61it/s]

     74%|██████████████████████████████████████████████████████████                     | 167/227 [00:46<00:16,  3.59it/s]

     74%|██████████████████████████████████████████████████████████▍                    | 168/227 [00:46<00:16,  3.59it/s]

     74%|██████████████████████████████████████████████████████████▊                    | 169/227 [00:46<00:16,  3.59it/s]

     75%|███████████████████████████████████████████████████████████▏                   | 170/227 [00:46<00:15,  3.59it/s]

     75%|███████████████████████████████████████████████████████████▌                   | 171/227 [00:47<00:15,  3.57it/s]

     76%|███████████████████████████████████████████████████████████▊                   | 172/227 [00:47<00:15,  3.58it/s]

     76%|████████████████████████████████████████████████████████████▏                  | 173/227 [00:47<00:15,  3.59it/s]

     77%|████████████████████████████████████████████████████████████▌                  | 174/227 [00:47<00:14,  3.60it/s]

     77%|████████████████████████████████████████████████████████████▉                  | 175/227 [00:48<00:14,  3.59it/s]

     78%|█████████████████████████████████████████████████████████████▎                 | 176/227 [00:48<00:14,  3.60it/s]

     78%|█████████████████████████████████████████████████████████████▌                 | 177/227 [00:48<00:13,  3.61it/s]

     78%|█████████████████████████████████████████████████████████████▉                 | 178/227 [00:49<00:13,  3.61it/s]

     79%|██████████████████████████████████████████████████████████████▎                | 179/227 [00:49<00:13,  3.62it/s]

     79%|██████████████████████████████████████████████████████████████▋                | 180/227 [00:49<00:12,  3.62it/s]

     80%|██████████████████████████████████████████████████████████████▉                | 181/227 [00:49<00:12,  3.61it/s]

     80%|███████████████████████████████████████████████████████████████▎               | 182/227 [00:50<00:12,  3.60it/s]

     81%|███████████████████████████████████████████████████████████████▋               | 183/227 [00:50<00:12,  3.59it/s]

     81%|████████████████████████████████████████████████████████████████               | 184/227 [00:50<00:11,  3.59it/s]

     81%|████████████████████████████████████████████████████████████████▍              | 185/227 [00:51<00:11,  3.56it/s]

     82%|████████████████████████████████████████████████████████████████▋              | 186/227 [00:51<00:11,  3.56it/s]

     82%|█████████████████████████████████████████████████████████████████              | 187/227 [00:51<00:11,  3.56it/s]

     83%|█████████████████████████████████████████████████████████████████▍             | 188/227 [00:51<00:10,  3.58it/s]

     83%|█████████████████████████████████████████████████████████████████▊             | 189/227 [00:52<00:10,  3.58it/s]

     84%|██████████████████████████████████████████████████████████████████             | 190/227 [00:52<00:10,  3.59it/s]

     84%|██████████████████████████████████████████████████████████████████▍            | 191/227 [00:52<00:10,  3.60it/s]

     85%|██████████████████████████████████████████████████████████████████▊            | 192/227 [00:52<00:09,  3.60it/s]

     85%|███████████████████████████████████████████████████████████████████▏           | 193/227 [00:53<00:09,  3.60it/s]

     85%|███████████████████████████████████████████████████████████████████▌           | 194/227 [00:53<00:09,  3.60it/s]

     86%|███████████████████████████████████████████████████████████████████▊           | 195/227 [00:53<00:08,  3.60it/s]

     86%|████████████████████████████████████████████████████████████████████▏          | 196/227 [00:54<00:08,  3.59it/s]

     87%|████████████████████████████████████████████████████████████████████▌          | 197/227 [00:54<00:08,  3.60it/s]

     87%|████████████████████████████████████████████████████████████████████▉          | 198/227 [00:54<00:08,  3.61it/s]

     88%|█████████████████████████████████████████████████████████████████████▎         | 199/227 [00:54<00:07,  3.61it/s]

     88%|█████████████████████████████████████████████████████████████████████▌         | 200/227 [00:55<00:07,  3.61it/s]

     89%|█████████████████████████████████████████████████████████████████████▉         | 201/227 [00:55<00:07,  3.60it/s]

     89%|██████████████████████████████████████████████████████████████████████▎        | 202/227 [00:55<00:06,  3.59it/s]

     89%|██████████████████████████████████████████████████████████████████████▋        | 203/227 [00:56<00:06,  3.58it/s]

     90%|██████████████████████████████████████████████████████████████████████▉        | 204/227 [00:56<00:06,  3.57it/s]

     90%|███████████████████████████████████████████████████████████████████████▎       | 205/227 [00:56<00:06,  3.56it/s]

     91%|███████████████████████████████████████████████████████████████████████▋       | 206/227 [00:56<00:05,  3.57it/s]

     91%|████████████████████████████████████████████████████████████████████████       | 207/227 [00:57<00:05,  3.58it/s]

     92%|████████████████████████████████████████████████████████████████████████▍      | 208/227 [00:57<00:05,  3.56it/s]

     92%|████████████████████████████████████████████████████████████████████████▋      | 209/227 [00:57<00:05,  3.56it/s]

     93%|█████████████████████████████████████████████████████████████████████████      | 210/227 [00:58<00:04,  3.56it/s]

     93%|█████████████████████████████████████████████████████████████████████████▍     | 211/227 [00:58<00:04,  3.54it/s]

     93%|█████████████████████████████████████████████████████████████████████████▊     | 212/227 [00:58<00:04,  3.54it/s]

     94%|██████████████████████████████████████████████████████████████████████████▏    | 213/227 [00:58<00:03,  3.54it/s]

     94%|██████████████████████████████████████████████████████████████████████████▍    | 214/227 [00:59<00:03,  3.55it/s]

     95%|██████████████████████████████████████████████████████████████████████████▊    | 215/227 [00:59<00:03,  3.55it/s]

     95%|███████████████████████████████████████████████████████████████████████████▏   | 216/227 [00:59<00:03,  3.56it/s]

     96%|███████████████████████████████████████████████████████████████████████████▌   | 217/227 [00:59<00:02,  3.57it/s]

     96%|███████████████████████████████████████████████████████████████████████████▊   | 218/227 [01:00<00:02,  3.56it/s]

     96%|████████████████████████████████████████████████████████████████████████████▏  | 219/227 [01:00<00:02,  3.57it/s]

     97%|████████████████████████████████████████████████████████████████████████████▌  | 220/227 [01:00<00:01,  3.57it/s]

     97%|████████████████████████████████████████████████████████████████████████████▉  | 221/227 [01:01<00:01,  3.58it/s]

     98%|█████████████████████████████████████████████████████████████████████████████▎ | 222/227 [01:01<00:01,  3.57it/s]

     98%|█████████████████████████████████████████████████████████████████████████████▌ | 223/227 [01:01<00:01,  3.56it/s]

     99%|█████████████████████████████████████████████████████████████████████████████▉ | 224/227 [01:01<00:00,  3.56it/s]

     99%|██████████████████████████████████████████████████████████████████████████████▎| 225/227 [01:02<00:00,  3.53it/s]

    100%|██████████████████████████████████████████████████████████████████████████████▋| 226/227 [01:02<00:00,  3.55it/s]

    100%|███████████████████████████████████████████████████████████████████████████████| 227/227 [01:02<00:00,  3.56it/s]

    100%|███████████████████████████████████████████████████████████████████████████████| 227/227 [01:02<00:00,  3.61it/s]

    



```python
tfscores = tfscores_exogenous.merge(tfscores_endogenous, on="tf")
```


```python
(tfscores["target_03"] < tfscores["q95"]).mean()
```




    np.float64(0.2555066079295154)




```python
tfscores.query("target_03 < q95").sort_values("target_03")
# tfscores.query("target_03 > q95")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tf</th>
      <th>target_03</th>
      <th>max_dose</th>
      <th>never_reached</th>
      <th>q90</th>
      <th>q95</th>
      <th>q99</th>
      <th>q10</th>
      <th>q05</th>
      <th>mean</th>
      <th>med</th>
      <th>cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>155</th>
      <td>Yap1</td>
      <td>0.072797</td>
      <td>1.201152</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.700772</td>
      <td>0.828469</td>
      <td>0.012021</td>
      <td>0.000000</td>
      <td>0.350047</td>
      <td>0.318454</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>26</th>
      <td>Myod1</td>
      <td>0.080944</td>
      <td>1.335579</td>
      <td>False</td>
      <td>0.351850</td>
      <td>0.461144</td>
      <td>0.826536</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.091071</td>
      <td>0.000000</td>
      <td>skeletal muscle satellite</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Myog</td>
      <td>0.083126</td>
      <td>2.057374</td>
      <td>False</td>
      <td>0.293836</td>
      <td>0.366375</td>
      <td>0.471864</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.126268</td>
      <td>0.095310</td>
      <td>myofibroblast</td>
    </tr>
    <tr>
      <th>148</th>
      <td>Runx2</td>
      <td>0.089475</td>
      <td>2.214504</td>
      <td>False</td>
      <td>2.297560</td>
      <td>2.469549</td>
      <td>3.152662</td>
      <td>0.604156</td>
      <td>0.388217</td>
      <td>1.354819</td>
      <td>1.280934</td>
      <td>plasmacytoid dendritic</td>
    </tr>
    <tr>
      <th>158</th>
      <td>Pou3f2</td>
      <td>0.094864</td>
      <td>0.284591</td>
      <td>False</td>
      <td>0.266785</td>
      <td>0.323203</td>
      <td>0.452526</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.112960</td>
      <td>0.071038</td>
      <td>unknown</td>
    </tr>
    <tr>
      <th>179</th>
      <td>T</td>
      <td>0.097016</td>
      <td>0.960462</td>
      <td>False</td>
      <td>1.229361</td>
      <td>1.330191</td>
      <td>1.444994</td>
      <td>0.386110</td>
      <td>0.365783</td>
      <td>0.755934</td>
      <td>0.693147</td>
      <td>ciliated</td>
    </tr>
    <tr>
      <th>138</th>
      <td>Gata4</td>
      <td>0.101853</td>
      <td>0.373462</td>
      <td>False</td>
      <td>0.782097</td>
      <td>0.981183</td>
      <td>1.271637</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.387520</td>
      <td>0.336472</td>
      <td>pancreatic acinar</td>
    </tr>
    <tr>
      <th>122</th>
      <td>Tfeb</td>
      <td>0.106412</td>
      <td>1.053480</td>
      <td>False</td>
      <td>0.510826</td>
      <td>0.583754</td>
      <td>0.847298</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.204121</td>
      <td>0.182322</td>
      <td>kidney collecting duct intercalated</td>
    </tr>
    <tr>
      <th>169</th>
      <td>Klf4</td>
      <td>0.112454</td>
      <td>2.226584</td>
      <td>False</td>
      <td>0.256528</td>
      <td>0.275765</td>
      <td>0.422939</td>
      <td>0.021746</td>
      <td>0.013360</td>
      <td>0.131299</td>
      <td>0.121885</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>157</th>
      <td>Foxa1</td>
      <td>0.136052</td>
      <td>1.036092</td>
      <td>False</td>
      <td>0.287682</td>
      <td>0.326306</td>
      <td>0.455033</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.098589</td>
      <td>0.000000</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>221</th>
      <td>Krba1</td>
      <td>0.138016</td>
      <td>0.546544</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.287682</td>
      <td>0.549113</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.074853</td>
      <td>0.000000</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>171</th>
      <td>Ascl1</td>
      <td>0.145811</td>
      <td>1.031091</td>
      <td>False</td>
      <td>0.581701</td>
      <td>0.605166</td>
      <td>0.634497</td>
      <td>0.098541</td>
      <td>0.089642</td>
      <td>0.326500</td>
      <td>0.292447</td>
      <td>lung neuroendocrine</td>
    </tr>
    <tr>
      <th>119</th>
      <td>Zkscan1</td>
      <td>0.148801</td>
      <td>0.409203</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.315483</td>
      <td>0.455033</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.089902</td>
      <td>0.000000</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>121</th>
      <td>Tgif2</td>
      <td>0.154046</td>
      <td>0.311236</td>
      <td>True</td>
      <td>0.161050</td>
      <td>0.287682</td>
      <td>0.385457</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.039897</td>
      <td>0.000000</td>
      <td>primordial germ</td>
    </tr>
    <tr>
      <th>172</th>
      <td>Pdx1</td>
      <td>0.157709</td>
      <td>1.951646</td>
      <td>False</td>
      <td>0.126203</td>
      <td>0.221482</td>
      <td>0.479188</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.063994</td>
      <td>0.030481</td>
      <td>type B pancreatic</td>
    </tr>
    <tr>
      <th>178</th>
      <td>Tbx5</td>
      <td>0.163624</td>
      <td>1.079920</td>
      <td>False</td>
      <td>0.367725</td>
      <td>0.502546</td>
      <td>0.791305</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.119887</td>
      <td>0.000000</td>
      <td>cardiac muscle</td>
    </tr>
    <tr>
      <th>167</th>
      <td>Sox2</td>
      <td>0.165698</td>
      <td>0.964948</td>
      <td>False</td>
      <td>0.307908</td>
      <td>0.336472</td>
      <td>0.510826</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.120792</td>
      <td>0.095310</td>
      <td>astrocyte</td>
    </tr>
    <tr>
      <th>211</th>
      <td>Zbtb33</td>
      <td>0.177042</td>
      <td>1.095448</td>
      <td>True</td>
      <td>0.108518</td>
      <td>0.193752</td>
      <td>0.485516</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.031892</td>
      <td>0.000000</td>
      <td>mesothelial</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Pax9</td>
      <td>0.191140</td>
      <td>1.720261</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.948560</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.289472</td>
      <td>0.232153</td>
      <td>epithelial of parathyroid gland</td>
    </tr>
    <tr>
      <th>182</th>
      <td>Rogdi</td>
      <td>0.204452</td>
      <td>0.632523</td>
      <td>True</td>
      <td>0.669839</td>
      <td>0.781650</td>
      <td>1.074266</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.187035</td>
      <td>0.000000</td>
      <td>microglial</td>
    </tr>
    <tr>
      <th>110</th>
      <td>Grhl2</td>
      <td>0.216725</td>
      <td>1.950526</td>
      <td>False</td>
      <td>0.587787</td>
      <td>0.693147</td>
      <td>0.754807</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.245749</td>
      <td>0.223144</td>
      <td>kidney collecting duct intercalated</td>
    </tr>
    <tr>
      <th>87</th>
      <td>Hoxa9</td>
      <td>0.247128</td>
      <td>2.038804</td>
      <td>False</td>
      <td>0.458520</td>
      <td>0.510826</td>
      <td>0.651213</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.129027</td>
      <td>0.000000</td>
      <td>mesangial</td>
    </tr>
    <tr>
      <th>137</th>
      <td>Tead2</td>
      <td>0.253186</td>
      <td>1.089802</td>
      <td>True</td>
      <td>0.447333</td>
      <td>0.509307</td>
      <td>0.691324</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.179339</td>
      <td>0.174688</td>
      <td>iris pigment epithelial</td>
    </tr>
    <tr>
      <th>203</th>
      <td>Arnt2</td>
      <td>0.280690</td>
      <td>0.793952</td>
      <td>True</td>
      <td>0.747214</td>
      <td>0.805871</td>
      <td>1.156715</td>
      <td>0.064539</td>
      <td>0.000000</td>
      <td>0.373159</td>
      <td>0.321938</td>
      <td>inhibitory interneuron</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Hoxb6</td>
      <td>0.281746</td>
      <td>1.162203</td>
      <td>False</td>
      <td>0.606136</td>
      <td>0.693147</td>
      <td>0.836988</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.286272</td>
      <td>0.287682</td>
      <td>epithelial</td>
    </tr>
    <tr>
      <th>80</th>
      <td>Foxp2</td>
      <td>0.286871</td>
      <td>0.617395</td>
      <td>False</td>
      <td>2.332040</td>
      <td>2.531740</td>
      <td>3.008194</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.340013</td>
      <td>1.420405</td>
      <td>Purkinje</td>
    </tr>
    <tr>
      <th>135</th>
      <td>Hoxd10</td>
      <td>0.288097</td>
      <td>0.559247</td>
      <td>False</td>
      <td>0.294581</td>
      <td>0.405465</td>
      <td>0.521072</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.084409</td>
      <td>0.000000</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>214</th>
      <td>Dcun1d3</td>
      <td>0.288620</td>
      <td>0.952445</td>
      <td>True</td>
      <td>1.124428</td>
      <td>1.199922</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.504067</td>
      <td>0.470004</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Ets1</td>
      <td>0.297916</td>
      <td>1.282333</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.916291</td>
      <td>1.299283</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.361277</td>
      <td>0.287682</td>
      <td>endothelial of hepatic sinusoid</td>
    </tr>
    <tr>
      <th>164</th>
      <td>Etv6</td>
      <td>0.316011</td>
      <td>0.744883</td>
      <td>True</td>
      <td>1.568518</td>
      <td>1.690255</td>
      <td>1.934213</td>
      <td>0.463753</td>
      <td>0.304695</td>
      <td>1.021060</td>
      <td>0.990399</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>77</th>
      <td>Gata2</td>
      <td>0.324433</td>
      <td>1.036092</td>
      <td>False</td>
      <td>0.609121</td>
      <td>0.775114</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.317286</td>
      <td>0.287682</td>
      <td>mast</td>
    </tr>
    <tr>
      <th>177</th>
      <td>Lhx3</td>
      <td>0.344928</td>
      <td>0.853698</td>
      <td>False</td>
      <td>0.312077</td>
      <td>0.405465</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.130032</td>
      <td>0.075775</td>
      <td>cochlea auditory hair</td>
    </tr>
    <tr>
      <th>88</th>
      <td>Foxp1</td>
      <td>0.346856</td>
      <td>1.560854</td>
      <td>False</td>
      <td>2.327625</td>
      <td>2.411096</td>
      <td>2.755557</td>
      <td>0.845859</td>
      <td>0.421831</td>
      <td>1.540615</td>
      <td>1.531663</td>
      <td>medium spiny neuron</td>
    </tr>
    <tr>
      <th>149</th>
      <td>Pitx2</td>
      <td>0.353527</td>
      <td>1.129004</td>
      <td>False</td>
      <td>0.287682</td>
      <td>0.444379</td>
      <td>0.565369</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.106693</td>
      <td>0.000000</td>
      <td>pituitary gland</td>
    </tr>
    <tr>
      <th>128</th>
      <td>Fos</td>
      <td>0.359276</td>
      <td>2.540597</td>
      <td>False</td>
      <td>1.615536</td>
      <td>1.733731</td>
      <td>2.159150</td>
      <td>0.067067</td>
      <td>0.022001</td>
      <td>0.879291</td>
      <td>0.868055</td>
      <td>keratinocyte stem</td>
    </tr>
    <tr>
      <th>41</th>
      <td>Otx2</td>
      <td>0.381766</td>
      <td>1.453647</td>
      <td>False</td>
      <td>1.055107</td>
      <td>1.130220</td>
      <td>1.265305</td>
      <td>0.202733</td>
      <td>0.154151</td>
      <td>0.613967</td>
      <td>0.649642</td>
      <td>photoreceptor</td>
    </tr>
    <tr>
      <th>17</th>
      <td>Meox2</td>
      <td>0.388429</td>
      <td>0.801136</td>
      <td>False</td>
      <td>0.854641</td>
      <td>1.029619</td>
      <td>1.498894</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.314086</td>
      <td>0.168236</td>
      <td>vein endothelial</td>
    </tr>
    <tr>
      <th>145</th>
      <td>Onecut2</td>
      <td>0.390056</td>
      <td>1.485213</td>
      <td>False</td>
      <td>1.098612</td>
      <td>1.193846</td>
      <td>1.394299</td>
      <td>0.287682</td>
      <td>0.188647</td>
      <td>0.640803</td>
      <td>0.587787</td>
      <td>retina horizontal</td>
    </tr>
    <tr>
      <th>186</th>
      <td>Egr1</td>
      <td>0.429298</td>
      <td>2.656283</td>
      <td>False</td>
      <td>0.442853</td>
      <td>0.472550</td>
      <td>0.509672</td>
      <td>0.237715</td>
      <td>0.234438</td>
      <td>0.317721</td>
      <td>0.261314</td>
      <td>lung neuroendocrine</td>
    </tr>
    <tr>
      <th>202</th>
      <td>Arnt</td>
      <td>0.463743</td>
      <td>0.646627</td>
      <td>True</td>
      <td>0.772667</td>
      <td>0.847298</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.378255</td>
      <td>0.287682</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>104</th>
      <td>Klf12</td>
      <td>0.466078</td>
      <td>1.098612</td>
      <td>False</td>
      <td>1.349071</td>
      <td>1.451965</td>
      <td>1.861387</td>
      <td>0.223144</td>
      <td>0.170172</td>
      <td>0.779011</td>
      <td>0.742407</td>
      <td>inhibitory interneuron</td>
    </tr>
    <tr>
      <th>131</th>
      <td>Tbx3</td>
      <td>0.472053</td>
      <td>0.994325</td>
      <td>False</td>
      <td>0.591457</td>
      <td>0.693147</td>
      <td>0.998457</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.306270</td>
      <td>0.287682</td>
      <td>kidney collecting duct principal</td>
    </tr>
    <tr>
      <th>160</th>
      <td>Cebpa</td>
      <td>0.499775</td>
      <td>2.356083</td>
      <td>False</td>
      <td>0.447609</td>
      <td>0.510826</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.148580</td>
      <td>0.000000</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>46</th>
      <td>Hmg20a</td>
      <td>0.521596</td>
      <td>1.324052</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.847298</td>
      <td>1.092401</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.305195</td>
      <td>0.287682</td>
      <td>epithelial</td>
    </tr>
    <tr>
      <th>123</th>
      <td>Tbx15</td>
      <td>0.547133</td>
      <td>0.547133</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.846121</td>
      <td>1.430029</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.223232</td>
      <td>0.000000</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>33</th>
      <td>Gata3</td>
      <td>0.613205</td>
      <td>1.517683</td>
      <td>False</td>
      <td>1.271078</td>
      <td>1.511091</td>
      <td>1.795766</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.654521</td>
      <td>0.559616</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>142</th>
      <td>Zim1</td>
      <td>0.617815</td>
      <td>0.617815</td>
      <td>True</td>
      <td>0.591457</td>
      <td>0.693147</td>
      <td>0.975666</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.188605</td>
      <td>0.060625</td>
      <td>intrahepatic cholangiocyte</td>
    </tr>
    <tr>
      <th>73</th>
      <td>Etv1</td>
      <td>0.620985</td>
      <td>2.195625</td>
      <td>False</td>
      <td>1.081210</td>
      <td>1.098612</td>
      <td>1.299283</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.525133</td>
      <td>0.510826</td>
      <td>cerebellar granule</td>
    </tr>
    <tr>
      <th>30</th>
      <td>Sall3</td>
      <td>0.624030</td>
      <td>1.188057</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.713958</td>
      <td>1.109966</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.232762</td>
      <td>0.182322</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>140</th>
      <td>Trerf1</td>
      <td>0.643550</td>
      <td>0.643550</td>
      <td>True</td>
      <td>1.697032</td>
      <td>1.805817</td>
      <td>2.386003</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.789126</td>
      <td>0.693147</td>
      <td>monocyte</td>
    </tr>
    <tr>
      <th>139</th>
      <td>Tfap2b</td>
      <td>0.723630</td>
      <td>0.863126</td>
      <td>False</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>1.307984</td>
      <td>0.154151</td>
      <td>0.000000</td>
      <td>0.575194</td>
      <td>0.510826</td>
      <td>retina horizontal</td>
    </tr>
    <tr>
      <th>40</th>
      <td>Tead1</td>
      <td>0.761400</td>
      <td>1.177790</td>
      <td>False</td>
      <td>2.436508</td>
      <td>2.527240</td>
      <td>2.732840</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.224665</td>
      <td>1.299283</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>79</th>
      <td>Zfp800</td>
      <td>0.794269</td>
      <td>0.947381</td>
      <td>False</td>
      <td>0.801948</td>
      <td>0.862520</td>
      <td>1.005754</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.307149</td>
      <td>0.223144</td>
      <td>mast</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Nfia</td>
      <td>0.934744</td>
      <td>1.779609</td>
      <td>False</td>
      <td>2.761583</td>
      <td>2.879069</td>
      <td>3.056618</td>
      <td>0.947824</td>
      <td>0.638451</td>
      <td>1.886274</td>
      <td>1.945910</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>53</th>
      <td>Bach2</td>
      <td>1.029138</td>
      <td>1.756633</td>
      <td>True</td>
      <td>2.711835</td>
      <td>2.780516</td>
      <td>3.208268</td>
      <td>0.311239</td>
      <td>0.135593</td>
      <td>1.378633</td>
      <td>1.252763</td>
      <td>B</td>
    </tr>
    <tr>
      <th>185</th>
      <td>Meis2</td>
      <td>1.068696</td>
      <td>3.648308</td>
      <td>False</td>
      <td>1.898362</td>
      <td>2.030316</td>
      <td>2.414269</td>
      <td>0.890606</td>
      <td>0.772997</td>
      <td>1.362172</td>
      <td>1.336012</td>
      <td>medium spiny neuron</td>
    </tr>
    <tr>
      <th>47</th>
      <td>Zfp462</td>
      <td>1.420589</td>
      <td>1.435085</td>
      <td>False</td>
      <td>1.299283</td>
      <td>1.495538</td>
      <td>1.792492</td>
      <td>0.356675</td>
      <td>0.223144</td>
      <td>0.827878</td>
      <td>0.806638</td>
      <td>oligodendrocyte precursor</td>
    </tr>
    <tr>
      <th>52</th>
      <td>Rbpjl</td>
      <td>1.850600</td>
      <td>1.850600</td>
      <td>True</td>
      <td>2.696246</td>
      <td>2.991554</td>
      <td>4.729635</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.098516</td>
      <td>0.829417</td>
      <td>pancreatic acinar</td>
    </tr>
  </tbody>
</table>
</div>




```python
# plotdata = tfscores.loc[~plotdata["never_reached"]]
plotdata = tfscores.loc[tfscores["tf"].isin(tfcapacity.query("category != 'low-capacity'").index)]

plotdata = plotdata.sort_values(["q95"])
plotdata["ix"] = range(len(plotdata))
fig, ax = plt.subplots(figsize = (2.5, len(plotdata) * 0.15))
for _, row in plotdata.iterrows():
    rect = plt.Rectangle(
        (row["q05"], row["ix"]),
        row["q95"] - row["q05"],
        0.95,
        color="#555",
        lw = 0,
    )
    ax.add_patch(rect)
    rect2 = plt.Rectangle(
        (0, row["ix"]),
        row["q05"],
        0.95,
        color="#CCC",
        lw = 0,
    )
    ax.add_patch(rect2)

    # med
    # plt.plot([row["mean"], row["mean"]], [row["ix"], row["ix"] + 1], color="#555")

    if not row["never_reached"]:
        color = "#2ECC40" if row["target_03"] < row["q95"] else "#FF4136"
        plt.scatter(row["target_03"], row["ix"]+0.5, color=color, marker = "o", zorder = 10)
        # connect dot to q95 if outside of range
        if row["target_03"] > row["q95"]:
            plt.plot(
                [row["q95"], row["target_03"]],
                [row["ix"] + 0.5, row["ix"] + 0.5],
                color="#444",linestyle = "-",
                zorder = -1
            )

    color = "#0074D9" if row["max_dose"] > row["q95"] else "#FF851B"
    # marker = "." if not row["never_reached"] else "x"
    marker = "|"
    ax.scatter(
        row["max_dose"], row["ix"] + 0.5, color=color, marker=marker, zorder=10, s = 15
    )
    ax.plot(
        [row["max_dose"], row["q95"]],
        [row["ix"] + 0.5, row["ix"] + 0.5],
        color="#0074D9",
        linestyle="-",
        alpha = 0.2,
        zorder=-1,
    )

    label = row["cell_type"]
    if len(label) > 27:
        label = "..." + label[-25:]
    ax.text(1, row["ix"] + 0.5, label, va="center", ha="left", transform=mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData), fontsize = 6, color = "#888")
# for ix in plotdata["ix"].iloc[::2]:
#     plt.axhspan(ix, ix+1, color="#EEE", zorder=-10)

ax.set_yticks(plotdata["ix"] + 0.5)
ax.set_yticklabels(plotdata["tf"], fontstyle = "italic", fontsize = 9)
ax.tick_params(axis="y", which="major", pad=0.2, length=0)
ax.set_xlim(0, 3.8)
ax.set_ylim(0., len(plotdata))
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlabel("TF/housekeeping expression")

import pathlib
plot_folder = pathlib.Path("plots_endogenous")
fig.savefig(plot_folder / "endogenous_global.pdf", bbox_inches="tight")
fig.savefig(plot_folder / "endogenous_global.png", bbox_inches="tight")
```


    
![png](output_27_0.png)
    



```python
tfscores_high_capacity = tfscores.loc[tfscores["tf"].isin(tfcapacity.query("category != 'low-capacity'").index)]
```


```python
print(f"""
Among all TFs (after pre-filtering, so 270), {(tfscores["target_03"] < tfscores["q95"]).mean():.0%}. Among all high-capacity TFs (76), {(tfscores_high_capacity["target_03"] < tfscores_high_capacity["q95"]).mean():.0%}.
""")
```

    
    Among all TFs (after pre-filtering, so 270), 26%. Among all high-capacity TFs (76), 47%.
    



```python
tfscores.query("target_03 < q95").sort_values("target_03")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tf</th>
      <th>target_03</th>
      <th>max_dose</th>
      <th>never_reached</th>
      <th>q90</th>
      <th>q95</th>
      <th>q99</th>
      <th>q10</th>
      <th>q05</th>
      <th>mean</th>
      <th>med</th>
      <th>cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>155</th>
      <td>Yap1</td>
      <td>0.072797</td>
      <td>1.201152</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.700772</td>
      <td>0.828469</td>
      <td>0.012021</td>
      <td>0.000000</td>
      <td>0.350047</td>
      <td>0.318454</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>26</th>
      <td>Myod1</td>
      <td>0.080944</td>
      <td>1.335579</td>
      <td>False</td>
      <td>0.351850</td>
      <td>0.461144</td>
      <td>0.826536</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.091071</td>
      <td>0.000000</td>
      <td>skeletal muscle satellite</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Myog</td>
      <td>0.083126</td>
      <td>2.057374</td>
      <td>False</td>
      <td>0.293836</td>
      <td>0.366375</td>
      <td>0.471864</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.126268</td>
      <td>0.095310</td>
      <td>myofibroblast</td>
    </tr>
    <tr>
      <th>148</th>
      <td>Runx2</td>
      <td>0.089475</td>
      <td>2.214504</td>
      <td>False</td>
      <td>2.297560</td>
      <td>2.469549</td>
      <td>3.152662</td>
      <td>0.604156</td>
      <td>0.388217</td>
      <td>1.354819</td>
      <td>1.280934</td>
      <td>plasmacytoid dendritic</td>
    </tr>
    <tr>
      <th>158</th>
      <td>Pou3f2</td>
      <td>0.094864</td>
      <td>0.284591</td>
      <td>False</td>
      <td>0.266785</td>
      <td>0.323203</td>
      <td>0.452526</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.112960</td>
      <td>0.071038</td>
      <td>unknown</td>
    </tr>
    <tr>
      <th>179</th>
      <td>T</td>
      <td>0.097016</td>
      <td>0.960462</td>
      <td>False</td>
      <td>1.229361</td>
      <td>1.330191</td>
      <td>1.444994</td>
      <td>0.386110</td>
      <td>0.365783</td>
      <td>0.755934</td>
      <td>0.693147</td>
      <td>ciliated</td>
    </tr>
    <tr>
      <th>138</th>
      <td>Gata4</td>
      <td>0.101853</td>
      <td>0.373462</td>
      <td>False</td>
      <td>0.782097</td>
      <td>0.981183</td>
      <td>1.271637</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.387520</td>
      <td>0.336472</td>
      <td>pancreatic acinar</td>
    </tr>
    <tr>
      <th>122</th>
      <td>Tfeb</td>
      <td>0.106412</td>
      <td>1.053480</td>
      <td>False</td>
      <td>0.510826</td>
      <td>0.583754</td>
      <td>0.847298</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.204121</td>
      <td>0.182322</td>
      <td>kidney collecting duct intercalated</td>
    </tr>
    <tr>
      <th>169</th>
      <td>Klf4</td>
      <td>0.112454</td>
      <td>2.226584</td>
      <td>False</td>
      <td>0.256528</td>
      <td>0.275765</td>
      <td>0.422939</td>
      <td>0.021746</td>
      <td>0.013360</td>
      <td>0.131299</td>
      <td>0.121885</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>157</th>
      <td>Foxa1</td>
      <td>0.136052</td>
      <td>1.036092</td>
      <td>False</td>
      <td>0.287682</td>
      <td>0.326306</td>
      <td>0.455033</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.098589</td>
      <td>0.000000</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>221</th>
      <td>Krba1</td>
      <td>0.138016</td>
      <td>0.546544</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.287682</td>
      <td>0.549113</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.074853</td>
      <td>0.000000</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>171</th>
      <td>Ascl1</td>
      <td>0.145811</td>
      <td>1.031091</td>
      <td>False</td>
      <td>0.581701</td>
      <td>0.605166</td>
      <td>0.634497</td>
      <td>0.098541</td>
      <td>0.089642</td>
      <td>0.326500</td>
      <td>0.292447</td>
      <td>lung neuroendocrine</td>
    </tr>
    <tr>
      <th>119</th>
      <td>Zkscan1</td>
      <td>0.148801</td>
      <td>0.409203</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.315483</td>
      <td>0.455033</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.089902</td>
      <td>0.000000</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>121</th>
      <td>Tgif2</td>
      <td>0.154046</td>
      <td>0.311236</td>
      <td>True</td>
      <td>0.161050</td>
      <td>0.287682</td>
      <td>0.385457</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.039897</td>
      <td>0.000000</td>
      <td>primordial germ</td>
    </tr>
    <tr>
      <th>172</th>
      <td>Pdx1</td>
      <td>0.157709</td>
      <td>1.951646</td>
      <td>False</td>
      <td>0.126203</td>
      <td>0.221482</td>
      <td>0.479188</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.063994</td>
      <td>0.030481</td>
      <td>type B pancreatic</td>
    </tr>
    <tr>
      <th>178</th>
      <td>Tbx5</td>
      <td>0.163624</td>
      <td>1.079920</td>
      <td>False</td>
      <td>0.367725</td>
      <td>0.502546</td>
      <td>0.791305</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.119887</td>
      <td>0.000000</td>
      <td>cardiac muscle</td>
    </tr>
    <tr>
      <th>167</th>
      <td>Sox2</td>
      <td>0.165698</td>
      <td>0.964948</td>
      <td>False</td>
      <td>0.307908</td>
      <td>0.336472</td>
      <td>0.510826</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.120792</td>
      <td>0.095310</td>
      <td>astrocyte</td>
    </tr>
    <tr>
      <th>211</th>
      <td>Zbtb33</td>
      <td>0.177042</td>
      <td>1.095448</td>
      <td>True</td>
      <td>0.108518</td>
      <td>0.193752</td>
      <td>0.485516</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.031892</td>
      <td>0.000000</td>
      <td>mesothelial</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Pax9</td>
      <td>0.191140</td>
      <td>1.720261</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.948560</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.289472</td>
      <td>0.232153</td>
      <td>epithelial of parathyroid gland</td>
    </tr>
    <tr>
      <th>182</th>
      <td>Rogdi</td>
      <td>0.204452</td>
      <td>0.632523</td>
      <td>True</td>
      <td>0.669839</td>
      <td>0.781650</td>
      <td>1.074266</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.187035</td>
      <td>0.000000</td>
      <td>microglial</td>
    </tr>
    <tr>
      <th>110</th>
      <td>Grhl2</td>
      <td>0.216725</td>
      <td>1.950526</td>
      <td>False</td>
      <td>0.587787</td>
      <td>0.693147</td>
      <td>0.754807</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.245749</td>
      <td>0.223144</td>
      <td>kidney collecting duct intercalated</td>
    </tr>
    <tr>
      <th>87</th>
      <td>Hoxa9</td>
      <td>0.247128</td>
      <td>2.038804</td>
      <td>False</td>
      <td>0.458520</td>
      <td>0.510826</td>
      <td>0.651213</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.129027</td>
      <td>0.000000</td>
      <td>mesangial</td>
    </tr>
    <tr>
      <th>137</th>
      <td>Tead2</td>
      <td>0.253186</td>
      <td>1.089802</td>
      <td>True</td>
      <td>0.447333</td>
      <td>0.509307</td>
      <td>0.691324</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.179339</td>
      <td>0.174688</td>
      <td>iris pigment epithelial</td>
    </tr>
    <tr>
      <th>203</th>
      <td>Arnt2</td>
      <td>0.280690</td>
      <td>0.793952</td>
      <td>True</td>
      <td>0.747214</td>
      <td>0.805871</td>
      <td>1.156715</td>
      <td>0.064539</td>
      <td>0.000000</td>
      <td>0.373159</td>
      <td>0.321938</td>
      <td>inhibitory interneuron</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Hoxb6</td>
      <td>0.281746</td>
      <td>1.162203</td>
      <td>False</td>
      <td>0.606136</td>
      <td>0.693147</td>
      <td>0.836988</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.286272</td>
      <td>0.287682</td>
      <td>epithelial</td>
    </tr>
    <tr>
      <th>80</th>
      <td>Foxp2</td>
      <td>0.286871</td>
      <td>0.617395</td>
      <td>False</td>
      <td>2.332040</td>
      <td>2.531740</td>
      <td>3.008194</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.340013</td>
      <td>1.420405</td>
      <td>Purkinje</td>
    </tr>
    <tr>
      <th>135</th>
      <td>Hoxd10</td>
      <td>0.288097</td>
      <td>0.559247</td>
      <td>False</td>
      <td>0.294581</td>
      <td>0.405465</td>
      <td>0.521072</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.084409</td>
      <td>0.000000</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>214</th>
      <td>Dcun1d3</td>
      <td>0.288620</td>
      <td>0.952445</td>
      <td>True</td>
      <td>1.124428</td>
      <td>1.199922</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.504067</td>
      <td>0.470004</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Ets1</td>
      <td>0.297916</td>
      <td>1.282333</td>
      <td>False</td>
      <td>0.693147</td>
      <td>0.916291</td>
      <td>1.299283</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.361277</td>
      <td>0.287682</td>
      <td>endothelial of hepatic sinusoid</td>
    </tr>
    <tr>
      <th>164</th>
      <td>Etv6</td>
      <td>0.316011</td>
      <td>0.744883</td>
      <td>True</td>
      <td>1.568518</td>
      <td>1.690255</td>
      <td>1.934213</td>
      <td>0.463753</td>
      <td>0.304695</td>
      <td>1.021060</td>
      <td>0.990399</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>77</th>
      <td>Gata2</td>
      <td>0.324433</td>
      <td>1.036092</td>
      <td>False</td>
      <td>0.609121</td>
      <td>0.775114</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.317286</td>
      <td>0.287682</td>
      <td>mast</td>
    </tr>
    <tr>
      <th>177</th>
      <td>Lhx3</td>
      <td>0.344928</td>
      <td>0.853698</td>
      <td>False</td>
      <td>0.312077</td>
      <td>0.405465</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.130032</td>
      <td>0.075775</td>
      <td>cochlea auditory hair</td>
    </tr>
    <tr>
      <th>88</th>
      <td>Foxp1</td>
      <td>0.346856</td>
      <td>1.560854</td>
      <td>False</td>
      <td>2.327625</td>
      <td>2.411096</td>
      <td>2.755557</td>
      <td>0.845859</td>
      <td>0.421831</td>
      <td>1.540615</td>
      <td>1.531663</td>
      <td>medium spiny neuron</td>
    </tr>
    <tr>
      <th>149</th>
      <td>Pitx2</td>
      <td>0.353527</td>
      <td>1.129004</td>
      <td>False</td>
      <td>0.287682</td>
      <td>0.444379</td>
      <td>0.565369</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.106693</td>
      <td>0.000000</td>
      <td>pituitary gland</td>
    </tr>
    <tr>
      <th>128</th>
      <td>Fos</td>
      <td>0.359276</td>
      <td>2.540597</td>
      <td>False</td>
      <td>1.615536</td>
      <td>1.733731</td>
      <td>2.159150</td>
      <td>0.067067</td>
      <td>0.022001</td>
      <td>0.879291</td>
      <td>0.868055</td>
      <td>keratinocyte stem</td>
    </tr>
    <tr>
      <th>41</th>
      <td>Otx2</td>
      <td>0.381766</td>
      <td>1.453647</td>
      <td>False</td>
      <td>1.055107</td>
      <td>1.130220</td>
      <td>1.265305</td>
      <td>0.202733</td>
      <td>0.154151</td>
      <td>0.613967</td>
      <td>0.649642</td>
      <td>photoreceptor</td>
    </tr>
    <tr>
      <th>17</th>
      <td>Meox2</td>
      <td>0.388429</td>
      <td>0.801136</td>
      <td>False</td>
      <td>0.854641</td>
      <td>1.029619</td>
      <td>1.498894</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.314086</td>
      <td>0.168236</td>
      <td>vein endothelial</td>
    </tr>
    <tr>
      <th>145</th>
      <td>Onecut2</td>
      <td>0.390056</td>
      <td>1.485213</td>
      <td>False</td>
      <td>1.098612</td>
      <td>1.193846</td>
      <td>1.394299</td>
      <td>0.287682</td>
      <td>0.188647</td>
      <td>0.640803</td>
      <td>0.587787</td>
      <td>retina horizontal</td>
    </tr>
    <tr>
      <th>186</th>
      <td>Egr1</td>
      <td>0.429298</td>
      <td>2.656283</td>
      <td>False</td>
      <td>0.442853</td>
      <td>0.472550</td>
      <td>0.509672</td>
      <td>0.237715</td>
      <td>0.234438</td>
      <td>0.317721</td>
      <td>0.261314</td>
      <td>lung neuroendocrine</td>
    </tr>
    <tr>
      <th>202</th>
      <td>Arnt</td>
      <td>0.463743</td>
      <td>0.646627</td>
      <td>True</td>
      <td>0.772667</td>
      <td>0.847298</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.378255</td>
      <td>0.287682</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>104</th>
      <td>Klf12</td>
      <td>0.466078</td>
      <td>1.098612</td>
      <td>False</td>
      <td>1.349071</td>
      <td>1.451965</td>
      <td>1.861387</td>
      <td>0.223144</td>
      <td>0.170172</td>
      <td>0.779011</td>
      <td>0.742407</td>
      <td>inhibitory interneuron</td>
    </tr>
    <tr>
      <th>131</th>
      <td>Tbx3</td>
      <td>0.472053</td>
      <td>0.994325</td>
      <td>False</td>
      <td>0.591457</td>
      <td>0.693147</td>
      <td>0.998457</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.306270</td>
      <td>0.287682</td>
      <td>kidney collecting duct principal</td>
    </tr>
    <tr>
      <th>160</th>
      <td>Cebpa</td>
      <td>0.499775</td>
      <td>2.356083</td>
      <td>False</td>
      <td>0.447609</td>
      <td>0.510826</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.148580</td>
      <td>0.000000</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>46</th>
      <td>Hmg20a</td>
      <td>0.521596</td>
      <td>1.324052</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.847298</td>
      <td>1.092401</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.305195</td>
      <td>0.287682</td>
      <td>epithelial</td>
    </tr>
    <tr>
      <th>123</th>
      <td>Tbx15</td>
      <td>0.547133</td>
      <td>0.547133</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.846121</td>
      <td>1.430029</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.223232</td>
      <td>0.000000</td>
      <td>fibroblast</td>
    </tr>
    <tr>
      <th>33</th>
      <td>Gata3</td>
      <td>0.613205</td>
      <td>1.517683</td>
      <td>False</td>
      <td>1.271078</td>
      <td>1.511091</td>
      <td>1.795766</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.654521</td>
      <td>0.559616</td>
      <td>mammary gland epithelial</td>
    </tr>
    <tr>
      <th>142</th>
      <td>Zim1</td>
      <td>0.617815</td>
      <td>0.617815</td>
      <td>True</td>
      <td>0.591457</td>
      <td>0.693147</td>
      <td>0.975666</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.188605</td>
      <td>0.060625</td>
      <td>intrahepatic cholangiocyte</td>
    </tr>
    <tr>
      <th>73</th>
      <td>Etv1</td>
      <td>0.620985</td>
      <td>2.195625</td>
      <td>False</td>
      <td>1.081210</td>
      <td>1.098612</td>
      <td>1.299283</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.525133</td>
      <td>0.510826</td>
      <td>cerebellar granule</td>
    </tr>
    <tr>
      <th>30</th>
      <td>Sall3</td>
      <td>0.624030</td>
      <td>1.188057</td>
      <td>True</td>
      <td>0.693147</td>
      <td>0.713958</td>
      <td>1.109966</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.232762</td>
      <td>0.182322</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>140</th>
      <td>Trerf1</td>
      <td>0.643550</td>
      <td>0.643550</td>
      <td>True</td>
      <td>1.697032</td>
      <td>1.805817</td>
      <td>2.386003</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.789126</td>
      <td>0.693147</td>
      <td>monocyte</td>
    </tr>
    <tr>
      <th>139</th>
      <td>Tfap2b</td>
      <td>0.723630</td>
      <td>0.863126</td>
      <td>False</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>1.307984</td>
      <td>0.154151</td>
      <td>0.000000</td>
      <td>0.575194</td>
      <td>0.510826</td>
      <td>retina horizontal</td>
    </tr>
    <tr>
      <th>40</th>
      <td>Tead1</td>
      <td>0.761400</td>
      <td>1.177790</td>
      <td>False</td>
      <td>2.436508</td>
      <td>2.527240</td>
      <td>2.732840</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.224665</td>
      <td>1.299283</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>79</th>
      <td>Zfp800</td>
      <td>0.794269</td>
      <td>0.947381</td>
      <td>False</td>
      <td>0.801948</td>
      <td>0.862520</td>
      <td>1.005754</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.307149</td>
      <td>0.223144</td>
      <td>mast</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Nfia</td>
      <td>0.934744</td>
      <td>1.779609</td>
      <td>False</td>
      <td>2.761583</td>
      <td>2.879069</td>
      <td>3.056618</td>
      <td>0.947824</td>
      <td>0.638451</td>
      <td>1.886274</td>
      <td>1.945910</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>53</th>
      <td>Bach2</td>
      <td>1.029138</td>
      <td>1.756633</td>
      <td>True</td>
      <td>2.711835</td>
      <td>2.780516</td>
      <td>3.208268</td>
      <td>0.311239</td>
      <td>0.135593</td>
      <td>1.378633</td>
      <td>1.252763</td>
      <td>B</td>
    </tr>
    <tr>
      <th>185</th>
      <td>Meis2</td>
      <td>1.068696</td>
      <td>3.648308</td>
      <td>False</td>
      <td>1.898362</td>
      <td>2.030316</td>
      <td>2.414269</td>
      <td>0.890606</td>
      <td>0.772997</td>
      <td>1.362172</td>
      <td>1.336012</td>
      <td>medium spiny neuron</td>
    </tr>
    <tr>
      <th>47</th>
      <td>Zfp462</td>
      <td>1.420589</td>
      <td>1.435085</td>
      <td>False</td>
      <td>1.299283</td>
      <td>1.495538</td>
      <td>1.792492</td>
      <td>0.356675</td>
      <td>0.223144</td>
      <td>0.827878</td>
      <td>0.806638</td>
      <td>oligodendrocyte precursor</td>
    </tr>
    <tr>
      <th>52</th>
      <td>Rbpjl</td>
      <td>1.850600</td>
      <td>1.850600</td>
      <td>True</td>
      <td>2.696246</td>
      <td>2.991554</td>
      <td>4.729635</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.098516</td>
      <td>0.829417</td>
      <td>pancreatic acinar</td>
    </tr>
  </tbody>
</table>
</div>




```python
(tfscores["target_03"] < tfscores["q95"]).mean()
```




    np.float64(0.2555066079295154)




```python
tfscores
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>tf</th>
      <th>target_03</th>
      <th>max_dose</th>
      <th>never_reached</th>
      <th>q90</th>
      <th>q95</th>
      <th>q99</th>
      <th>q10</th>
      <th>q05</th>
      <th>mean</th>
      <th>med</th>
      <th>cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Esr2</td>
      <td>0.161522</td>
      <td>2.665112</td>
      <td>False</td>
      <td>0.068993</td>
      <td>0.084786</td>
      <td>0.156541</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.025791</td>
      <td>0.014528</td>
      <td>granulosa</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Nkx2-6</td>
      <td>0.849551</td>
      <td>2.273123</td>
      <td>False</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.001806</td>
      <td>0.000000</td>
      <td>enteroendocrine</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Zscan18</td>
      <td>0.765531</td>
      <td>1.353350</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.287682</td>
      <td>0.462032</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.064622</td>
      <td>0.000000</td>
      <td>pinealocyte</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Nfia</td>
      <td>0.934744</td>
      <td>1.779609</td>
      <td>False</td>
      <td>2.761583</td>
      <td>2.879069</td>
      <td>3.056618</td>
      <td>0.947824</td>
      <td>0.638451</td>
      <td>1.886274</td>
      <td>1.945910</td>
      <td>adipocyte</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Bcl6b</td>
      <td>1.013195</td>
      <td>1.013195</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.332311</td>
      <td>0.510826</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.073665</td>
      <td>0.000000</td>
      <td>endothelial</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>222</th>
      <td>Zfpl1</td>
      <td>1.298489</td>
      <td>1.298489</td>
      <td>True</td>
      <td>0.130369</td>
      <td>0.149125</td>
      <td>0.172570</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.049554</td>
      <td>0.010479</td>
      <td>lung neuroendocrine</td>
    </tr>
    <tr>
      <th>223</th>
      <td>Zscan29</td>
      <td>3.139915</td>
      <td>3.139915</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.287682</td>
      <td>0.427591</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.070275</td>
      <td>0.000000</td>
      <td>kidney d. conv. tubule epithelial</td>
    </tr>
    <tr>
      <th>224</th>
      <td>Mxd1</td>
      <td>2.585711</td>
      <td>2.585711</td>
      <td>True</td>
      <td>0.540836</td>
      <td>0.648874</td>
      <td>0.783921</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.340120</td>
      <td>0.346574</td>
      <td>neutrophil</td>
    </tr>
    <tr>
      <th>225</th>
      <td>Zfp811</td>
      <td>1.734601</td>
      <td>1.734601</td>
      <td>True</td>
      <td>0.470004</td>
      <td>0.540100</td>
      <td>0.860651</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.192195</td>
      <td>0.154151</td>
      <td>retina horizontal</td>
    </tr>
    <tr>
      <th>226</th>
      <td>A430033K04Rik</td>
      <td>2.341661</td>
      <td>2.341661</td>
      <td>True</td>
      <td>0.287682</td>
      <td>0.336472</td>
      <td>0.510826</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.062872</td>
      <td>0.000000</td>
      <td>epithelial</td>
    </tr>
  </tbody>
</table>
<p>227 rows × 12 columns</p>
</div>




```python

```
