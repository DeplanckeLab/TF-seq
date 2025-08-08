# %% [markdown]
# # Dose levels reached in scTF-seq compared to other overexpression datasets

# In scTF-seq, we specifically chose a construct that can reach a variety of dose levels. A key challenge in many overexpression or CRISPRa studies is that dose levels are not high enough to reach physiological levels. To investigate this, we compare the dose levels reached in scTF-seq to those reached in other datasets, such as the MORF dataset (Joung et al. 2023).

# Given that we are working in another context then Joung et al. (mESCs vs C3H10), this comparison is not straightforward. Moreover, the overlap between overexpressed genes at high-enough power levels between both datasets is sadly very small. In the end, only one TF led to a useful comparison: Fos, which has a strong response in both mESCs and C3H10 cells. We will focus on DE genes that are shared between both datasets, and compare their differential expression at various dose levels.


# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import tqdm.auto as tqdm

import jax
import latenta as la
import polyptich as pp

pp.setup_ipython()

import eyck
import seaborn as sns

# %%
import pathlib
plots_folder = pathlib.Path("plots")
if not plots_folder.exists():
    plots_folder.mkdir(parents=True)
data_folder = pp.paths.get_data()

# %% [markdown]
# ## Download and load MORF

# We use data from the MORF dataset (Joung et al. 2023) to compare exogenous TF expression.

# %%
import scanpy as sc

# %%
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216595/suppl/GSE216595%5F180124%5Fperturb%2Eh5ad%2Egz"
if not pathlib.Path("GSE216481/GSE216595_180124_perturb.h5ad").exists():
    pathlib.Path("GSE216481").mkdir(exist_ok=True, parents=True)
    !wget {url} -O GSE216481/GSE216595_180124_perturb.h5ad.gz
    !gunzip GSE216481/GSE216595_180124_perturb.h5ad.gz

# %%
adata_morf = sc.read_h5ad("GSE216481/GSE216595_180124_perturb.h5ad")

# %%
adata2_morf = adata_morf[
    adata_morf.obs["TF"].str.contains("FOS-")
    | adata_morf.obs["TF"].str.contains("STAT3")
]
adata2_morf.obs["TF"].value_counts()

# %%
eyck.m.t.plot_umap(adata2_morf, color=["batch", "TF", "MKI67", "EMP3"]).display()

# %%
diffexp_morf = eyck.m.t.diffexp.compare_two_groups(
    adata2_morf, adata2_morf.obs["TF"].str.contains("FOS").values
)
# diffexp.loc["CAV1"]
diffexp_morf.loc["S100A6"]

# %% [markdown]
# Because not all cells that have high Fos barcode expression are acatually perturbed, we will only use the cells that cluster in a separate Fos-enriched cluster.

# %%
sc.pp.pca(adata2_morf)
sc.pp.neighbors(adata2_morf)
sc.tl.umap(adata2_morf)
sc.tl.leiden(adata2_morf, resolution=0.5)
adata2_morf = adata2_morf[adata2_morf.obs["leiden"].isin(["0", "1"])]

# %%
sc.pp.pca(adata2_morf)
sc.pp.neighbors(adata2_morf)
sc.tl.umap(adata2_morf)

# %%
sc.pl.umap(adata2_morf, color=["batch", "TF", "leiden", "MKI67", "EMP3", "CDK1"])

# %%
adata2_morf.obs["TF2"] = ["FOS" if x == "0" else "mCherry" for x in adata2_morf.obs["leiden"]]

# %%
diffexp_morf2 = eyck.m.t.diffexp.compare_two_groups(
    adata2_morf,
    adata2_morf.obs["TF2"].isin(["FOS"]),
    adata2_morf.obs["TF2"].isin(["mCherry"]),
)
diffexp_morf2["lfc"] = (
    np.array(
        np.array(adata2_morf.raw.to_adata().X[adata2_morf.obs["leiden"].isin(["0"]).values].mean(0))[0]
    )
    - np.array(
        np.array(adata2_morf.raw.to_adata().X[adata2_morf.obs["leiden"].isin(["1"]).values].mean(0))[0]
    )
)
diffexp_morf2 = diffexp_morf2.groupby("gene").first()


# %% [markdown]
# ## Load scTF-seq data

# %%
obs = pd.read_csv(data_folder / "obs.csv").rename(columns={"Unnamed: 0": "cell"}).set_index("cell")
mtx = scipy.io.mmread(data_folder / "matrix.mtx").T.tocsr()
obs["ix"] = range(obs.shape[0])

var = pd.read_csv(data_folder / "var.csv", index_col=0)

# %%
def extract_tf_dataset(mtx, obs, tf, batches=None):
    if isinstance(tf, str):
        tf = [tf]
    obs_oi = obs[obs["TF"].isin(tf)]
    if batches is not None:
        obs_oi = obs_oi[obs_oi["batch_overall"].isin(batches)]
    batches = obs_oi["batch_overall"].unique()
    obs_ref = obs[
        (obs["batch_overall"].isin(batches)) & (obs["TF"].isin(["D0_confluent"]))
    ]

    obs = pd.concat([obs_ref, obs_oi])
    obs = obs.loc[obs["Phase_corrected"] == "G1"]
    counts = mtx.tocsr()[obs["ix"].values, :]

    return counts, obs

# %%
import scanpy as sc

counts_oi, obs_oi = extract_tf_dataset(mtx, obs, ["Fos"], batches=["batch6"])

genes_oi = np.array(counts_oi.sum(0))[0] > 5
counts_oi = counts_oi[:, genes_oi]
var_oi = var.iloc[genes_oi]
len(var_oi), len(obs_oi)

# %%
adata = sc.AnnData(X=counts_oi, obs=obs_oi, var=var_oi)

# %%
counts = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# %%
oi = adata.obs["TF"] == "Fos"
noi = ~oi
diffexp = eyck.m.t.diffexp.compare_two_groups(adata, oi)
diffexp["lfc"] = np.array(adata.X[oi.values].mean(0))[0] - np.array(adata.X[noi.values].mean(0))[0]
diffexp = diffexp.groupby("symbol").first()

# %% [markdown]
# ## Compare fold-changes

# %%
fig, ax = plt.subplots()
sns.histplot(diffexp.query("pvals_adj < 0.05")["logfoldchanges"])
sns.histplot(diffexp_morf.query("pvals_adj < 0.05")["logfoldchanges"])
ax.set_xlim(-5, 5)

# %%
fig, ax = plt.subplots()
sns.histplot(diffexp.query("pvals_adj < 0.05")["lfc"], bins=np.linspace(-5, 5, 100))
sns.histplot(
    diffexp_morf2.query("pvals_adj < 0.05")["lfc"], bins=np.linspace(-5, 5, 100)
)
ax.set_xlim(-5, 5)

# %% [markdown]
# ## Extract shared DE genes

# Given that we work in different contexts, we cannot simply compare one set of DE genes with another. Instead, we will focus on the intersection of DE genes.

# %%
gene_mapping = pd.DataFrame({"gene_mouse": diffexp.index})
gene_mapping["gene_human"] = gene_mapping["gene_mouse"].str.upper()
gene_mapping["shared"] = gene_mapping["gene_human"].isin(diffexp_morf2.index)
gene_mapping = gene_mapping.loc[gene_mapping["shared"]]

diffexp_shared = pd.concat(
    [
        diffexp
        .loc[gene_mapping["gene_mouse"]]
        .rename(
            columns={
                "logfoldchanges": "logfoldchanges_mouse",
                "pvals_adj": "pvals_adj_mouse",
                "scores": "scores_mouse",
            }
        )[["logfoldchanges_mouse", "pvals_adj_mouse", "scores_mouse"]]
        .reset_index(),
        diffexp_morf2.loc[gene_mapping["gene_human"]]
        .rename(
            columns={
                "logfoldchanges": "logfoldchanges_human",
                "pvals_adj": "pvals_adj_human",
                "scores": "scores_human",
            }
        )[["logfoldchanges_human", "pvals_adj_human", "scores_human"]]
        .reset_index(),
    ],
    axis=1,
    ignore_index=True,
)
diffexp_shared.columns = [
    "gene_mouse",
    "logfoldchanges_mouse",
    "pvals_adj_mouse",
    "scores_mouse",
    "gene_human",
    "logfoldchanges_human",
    "pvals_adj_human",
    "scores_human",
]
diffexp_shared_diff = diffexp_shared.query(
    "abs(scores_mouse) > 5 & abs(scores_human) > 5"
).dropna()
diffexp_shared_diff.shape

# %%
fig, ax = plt.subplots()
ax.scatter(
    diffexp_shared_diff["scores_mouse"],
    diffexp_shared_diff["scores_human"],
    s=1,
)

# %% [markdown]
# Shared DE gene fold-changes are moderately correlated.

# %%
np.corrcoef(
    diffexp_shared_diff["scores_mouse"],
    diffexp_shared_diff["scores_human"],
)

# %%
adata.obs["nCount_RNA"].mean(), adata2_morf.obs["n_counts"].mean()

# %%
diffexp_shared.query("abs(logfoldchanges_mouse) > log(2) & pvals_adj_mouse < 0.05").shape, diffexp_shared.query("abs(logfoldchanges_human) > log(2) & pvals_adj_human < 0.05").shape

# %%
diffexp_shared_diff.sort_values("pvals_adj_human", ascending = True).query("logfoldchanges_human > 1").query("logfoldchanges_mouse > 1")

# %%
# label = ["Cav1", "S100a6", "Rspo2"]
label = ["Emp3", "Cav1", "S100a6", "S100a4", "Anxa1", "Tgfb1"]

fig = pp.Figure()
ax = fig.main.add_right(pp.Panel((3, 2)))
plotdata = diffexp_shared.copy()
plotdata["x"] = plotdata["logfoldchanges_mouse"]
plotdata["y"] = -np.log(plotdata["pvals_adj_mouse"])
ax.scatter(
    plotdata["x"], plotdata["y"],
    s=1,
)
texts = []
for gene in label:
    ix = diffexp_shared["gene_mouse"] == gene
    texts.append(ax.text(
        plotdata["x"][ix],
        plotdata["y"][ix],
        gene,
        fontstyle = "italic",
        # bbox=dict(facecolor='white', pad=0.4),
    ))
    texts[-1].set_path_effects([mpl.patheffects.withStroke(linewidth=2, foreground='white')])
ax.set_xlim(-5, 5)
ax.set_xticks(np.log([1/16, 1/4, 1, 4, 16]))
ax.set_xticklabels(["1/16", "1/4", 1, 4, 16], fontsize = 10, rotation = 90)
ax.set_ylabel("-log(q-value)")
ax.set_ylim(0, 250)
# ax.axhline(-np.log(0.05), color='grey', linestyle='--')
ax.set_xlabel("fold-change")

fig.plot()
import adjustText
adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='red'), ax = ax, min_arrow_len =0)


ax = fig.main.add_right(pp.Panel((3, 2)))
plotdata = diffexp_shared.copy()
plotdata["x"] = plotdata["logfoldchanges_human"]
plotdata["y"] = -np.log(plotdata["pvals_adj_human"])
s = ax.scatter(
    plotdata["x"], plotdata["y"],
    s=1,
)
texts = []
for gene in label:
    ix = diffexp_shared["gene_human"] == gene.upper()
    texts.append(ax.text(
        plotdata["x"][ix],
        plotdata["y"][ix],
        gene.upper(),
        fontstyle = "italic",
        # bbox=dict(facecolor='white', pad=0.4),
    ))
    texts[-1].set_path_effects([mpl.patheffects.withStroke(linewidth=2, foreground='white')])
ax.set_xticks(np.log([1/16, 1/4, 1, 4, 16]))
ax.set_xticklabels(["1/16", "1/4", 1, 4, 16], fontsize = 10, rotation = 90)
ax.set_xlim(-5, 5)
ax.set_ylim(0, 250)
# ax.axhline(-np.log(0.05), color='grey', linestyle='--')

fig.plot()
import adjustText
adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='red'), ax = ax, min_arrow_len =0)

# ax = fig.main.add_right(pp.Panel((2, 2)))
# plotdata = diffexp_shared_diff.copy()
# plotdata["x"] = plotdata["scores_human"]
# plotdata["y"] = plotdata["scores_mouse"]
# s = ax.scatter(
#     plotdata["x"], plotdata["y"],
#     s=1,
# )
# texts = []
# for gene in label:
#     ix = diffexp_shared["gene_human"] == gene.upper()
#     texts.append(ax.text(
#         plotdata["x"][ix],
#         plotdata["y"][ix],
#         gene,
#         # bbox=dict(facecolor='white', pad=0.4),
#     ))
#     texts[-1].set_path_effects([mpl.patheffects.withStroke(linewidth=2, foreground='white')])
# ax.set_xticks(np.log([1/16, 1/4, 1, 4, 16]))
# ax.set_xticklabels(["1/16", "1/4", 1, 4, 16], fontsize = 10, rotation = 90)
# ax.axhline(-np.log(0.05), color='grey', linestyle='--')

# fig.plot()
# import adjustText
# adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='red'), ax = ax, min_arrow_len =0)
fig.display()

fig.savefig(plots_folder / "dose_diffexp_comparison.pdf")

# %%
fig = pp.Figure(pp.grid.Grid(padding_width = 0.1, padding_height = 0.4))

ax = fig.main.add_right(pp.Panel((1, 1)), row = 0)
adata.obs["TF2"] = ["FOS" if x == "Fos" else "mCherry" for x in adata.obs["TF"]]
eyck.m.t.plot_umap(adata, color="Dose", ax=ax, legend = False, cmap = "viridis", show_norm = False)
ax.set_title("scTF-seq")

ax = fig.main.add_right(pp.Panel((1, 1)), row = 1)
eyck.m.t.plot_umap(adata2_morf, color="TF2", ax=ax, legend = False)
ax.set_title("Joung et al. 2023")

for symbol in ["Emp3", "Cav1", "S100a6", "S100a4", "Anxa1", "Tgfb1"]:
# for symbol in ["Cav1", "Serpine1", "Marcks", "Kcnn4"]:
    ax = fig.main.add_right(pp.Panel((1, 1)), row = 0)
    eyck.m.t.plot_umap(adata, color=symbol, ax=ax, norms = ("q.05", "q.98"), show_norm = False)
    ax.set_title(f"{symbol}", fontstyle = "italic")
    ax = fig.main.add_right(pp.Panel((1, 1)), row = 1)
    eyck.m.t.plot_umap(adata2_morf.raw.to_adata(), color=symbol.upper(), ax=ax, norms = ("q.05", "q.98"), show_norm = False)
    ax.set_title(f"{symbol.upper()}", fontstyle = "italic")
fig.display()


fig.savefig(plots_folder / "dose_diffexp_examples.pdf")

# %%
(
    np.exp(diffexp_shared_diff["logfoldchanges_mouse"].abs().mean()),
    np.exp(diffexp_shared_diff["logfoldchanges_human"].abs().mean()),
)

# %%
diffexp_shared_diff.sort_values("scores_human")

# %% [markdown]
# ## Kinetic over various dose levels

# %%
distributions = []
ns = []
genes_oi = diffexp_shared_diff["gene_mouse"]
# genes_oi = diffexp.query("pvals_adj < 0.05 & abs(logfoldchanges) > log(2)").index
for dose_cutoff in tqdm.tqdm(np.arange(2.5, 8, 0.25)):
    oi = (adata.obs["TF"] == "Fos") & (np.log1p(adata.obs["Vector_UMI"]) < dose_cutoff)
    noi = adata.obs["TF"] != "Fos"
    diffexp_oi = eyck.m.t.diffexp.compare_two_groups(adata, oi)
    distributions.append(
        diffexp_oi.set_index("symbol").loc[genes_oi, "logfoldchanges"].values
    )
    ns.append(diffexp_oi.query("pvals_adj < 0.05 & abs(logfoldchanges) > log(2)").shape[0])

# %%
fig, ax = plt.subplots(figsize = (6, 2))

means = []
for x, dist in enumerate(distributions):
    dist_oi = np.abs(np.clip(dist, -3, 7))
    # vln = ax.violinplot(dist_oi, positions=[x], showextrema=False, showmedians = True)
    # for pc in vln['bodies']:
    #     pc.set_facecolor('#333')
    sns.boxplot(x=x, y=dist_oi, ax=ax, color = "tomato", showfliers = False, )
    sns.stripplot(x=x, y=dist_oi, ax=ax, s=1, color = "grey", zorder = 10, alpha = 0.5)
    means.append(np.median(dist_oi))
# ax.plot(means, color = "red")
ax.set_ylim(0, 2)

ax.set_yticks(np.log([1, 2, 4, 8]))
ax.set_yticklabels(["1", "2", "4", "8"], fontsize = 10)

x = x+1
dist = diffexp_shared_diff["logfoldchanges_human"].values
dist_oi = np.abs(np.clip(dist, -3, 7))
# ax.violinplot(dist_oi, positions=[x], showextrema=False, showmedians = True)
sns.boxplot(x=x+3, y=dist_oi, ax=ax, color = "cyan", showfliers = False)
sns.stripplot(x=x+3, y=dist_oi, ax=ax, s=1, zorder = 10, color = "grey", alpha = 0.5)
ax.axhline(np.median(dist_oi), zorder = 10, color = "cyan")
# ax.set_xlabel("scTF-seq", color = "tomato")

ticks = pd.DataFrame({"x": np.arange(0, len(means) + 1)})
ticks["label"] = (np.arange(2.5, 8, 0.25)).tolist() + ["Joung\n2023"]
ticks = ticks.iloc[::2]
ax.set_xticks(ticks["x"])
ax.set_xticklabels(ticks["label"], fontsize = 10)
ax.get_xticklabels()[-1].set_color("cyan")
ax.get_xticklabels()[-1].set_rotation(90)

fig.savefig(plots_folder / "dose_diffexp_kinetics.pdf", bbox_inches = "tight")

# %%
