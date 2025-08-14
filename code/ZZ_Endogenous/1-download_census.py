# %% [markdown]
# # Download census data

# %%
import polyptich as pp

import pickle

pp.setup_ipython()

import cellxgene_census

import pandas as pd

# %%
data_folder = pp.paths.get_data()

# %%
obs_alltfs = (
    pd.read_csv(data_folder / "obs.csv").rename(columns={"Unnamed: 0": "cell"}).set_index("cell")
)

# %% [markdown]
# Define a set of housekeeping genes to compare against.
# We explicitely chose a set of genes from different essential functions such as cytoskeleton, ribosomal proteins, glycolysis, ubiquitinilation, respiration, and translation.

# %%
housekeeping = [
    "Actb",
    "Gapdh",
    "Ppia",
    "Ppib",
    "Tubb5",
    "Rpl13",
    "Rpl13a",
    "Rpl19",
    "Ubc",
    "Gusb",
    "Ywhaz",
    "Eef1a1",
    "Pgk1",
]
tfs = obs_alltfs["TF"].unique().tolist()

pickle.dump(housekeeping, open("housekeeping.pkl", "wb"))
pickle.dump(tfs, open("tfs.pkl", "wb"))

# %%
census = cellxgene_census.open_soma()
mouse = census["census_data"]["mus_musculus"]
datasets = census["census_info"]["datasets"].read().concat().to_pandas()

# %%
pickle.dump(datasets, open("datasets.pkl", "wb"))

# %%
obs_df = (
    mouse.obs.read(
        column_names=["cell_type_ontology_term_id", "cell_type", "dataset_id", "is_primary_data"]
    )
    .concat()
    .to_pandas()
)
celltypes = obs_df.groupby(
    by=["cell_type_ontology_term_id", "cell_type", "dataset_id"],
    as_index=False,
    observed=True,
).size()

# %% [markdown]
# Select cell types of interest.
# By default, all!

# %%
celltypes_oi = celltypes

# %% [markdown]
# Select cells of interest. To speed up downloading, we take 200 random cells per dataset and cell type.

# %%
obs = cellxgene_census.get_obs(
    census,
    organism="mus_musculus",
    column_names=[
        "dataset_id",
        "cell_type_ontology_term_id",
        "cell_type",
        "is_primary_data",
        "soma_joinid",
    ],
    value_filter=f"(is_primary_data == True) & (cell_type_ontology_term_id in {celltypes_oi['cell_type_ontology_term_id'].tolist()})",
)
obs = (
    obs.groupby(["dataset_id", "cell_type"])
    .sample(200, replace=True, random_state=42)
    .drop_duplicates()
)

# %%
obs.groupby(["cell_type"], observed=True).size()

# %% [markdown]
# Download the data for only genes (=TFs) of interest.

# %%
genes = housekeeping + tfs
var = cellxgene_census.get_var(
    census,
    organism="mus_musculus",
    column_names=["feature_name", "feature_id", "soma_joinid"],
    value_filter=f"feature_name in {genes}",
)

# %%
adata = cellxgene_census.get_anndata(
    census=census,
    organism="mus_musculus",
    measurement_name="RNA",
    obs_coords=obs["soma_joinid"].values,
    obs_column_names=obs.columns,
    var_coords=var["soma_joinid"].values,
)

import pickle

pickle.dump(adata, open("adata2.pkl", "wb"))