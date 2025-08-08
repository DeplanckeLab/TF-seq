# Helper functions to plot endogenous vs exogenous data from scTF-seq
# Author: Wouter Saelens https://orcid.org/0000-0002-7114-6248

import scanpy as sc
import numpy as np
import pandas as pd
import eyck
import polyptich as pp


def smooth_spline_fit_se(x, y, x_smooth):
    import rpy2.robjects as ro

    ro.globalenv["x"] = ro.FloatVector(x)
    ro.globalenv["y"] = ro.FloatVector(y)
    ro.globalenv["x_pred"] = ro.FloatVector(x_smooth)
    script = """
    # Install and load the mgcv package if not yet done
    if (!require(mgcv)) {
    install.packages('mgcv')
    library(mgcv)
    }

    # Fit a GAM with a smoothing spline, just like smooth.spline
    gam_model <- gam(y ~ s(x, sp = 0.2), method = 'REML')

    # Make predictions
    y_pred <- predict(gam_model, newdata = data.frame(x = x_pred), type = 'response', se.fit = TRUE)

    # Extract predicted values and standard errors
    # fit <- y_pred
    # se <- y_pred
    fit <- y_pred$fit
    se <- y_pred$se.fit

    list(fit = fit, se = se)
    """
    out = ro.r(script)
    fit = np.array(out[0])
    se = np.array(out[1])
    return np.stack([fit, se]).T


def abbreviate_sentence(sentence: str, max_length: int) -> str:
    abbreviations = {
        # "extratelencephalic projecting": "ep.",
        "extratelencephalic ": "",
        "intratelencephalic ": "",
        "glutamatergic": "glutamat.",
        "convoluted": "conv.",
        "distal": "d.",
        "3-month-old stage": "3m",
        "20-month-old stage and over": ">20m",
        "18-month-old stage": "18m",
        "projecting ": "",
        " NucSeq": "",
        "meis2 expressing": "Meis2+",
        "pulmonary": "pulm.",
        " cell": "",
        "cortical": "cort.",
        "pvalb ": "",
        "GABAergic ": "",
        "metanephric ": "meta.",
        "lymphatic vessel ": "lymphatics",
    }

    for k, v in abbreviations.items():
        sentence = sentence.replace(k, v)
    shortened_sentence = sentence

    if len(shortened_sentence) <= max_length:
        return shortened_sentence

    shortened_words = sentence.split(" ")

    # If still too long, remove short words like "the", "is", etc.
    filtered_words = [
        word
        for word in shortened_words
        if len(word) > 3 or word in abbreviations.values()
    ]
    shortened_sentence = " ".join(filtered_words)

    if len(shortened_sentence) <= max_length:
        return shortened_sentence

    # If still too long, progressively remove words from the beginning
    while len(shortened_sentence) > max_length and filtered_words:
        filtered_words = filtered_words[1:]
        shortened_sentence = " ".join(filtered_words)

    return shortened_sentence


def extract_tf_dataset(mtx, obs, var, tf, housekeeping, batches=None):
    if isinstance(tf, str):
        tf = [tf]
    obs_oi = obs[obs["TF"].isin(tf)]
    if batches is not None:
        obs_oi = obs_oi[obs_oi["batch_overall"].isin(batches)]
    batches = obs_oi["batch_overall"].unique()
    obs_ref = obs[
        (obs["batch_overall"].isin(batches)) & (obs["TF"].isin(["D0", "D0_confluent"]))
    ]

    obs = pd.concat([obs_ref, obs_oi])
    obs = obs.loc[obs["Phase_corrected"] == "G1"]
    counts = mtx.tocsr()[obs["ix"].values, :]

    adata = sc.AnnData(X=counts, obs=obs, var=var)

    adata.layers["norm"] = adata.X.copy()
    adata.obs["norm"] = 1e1 / (
        adata[:, eyck.m.t.gene_id(adata.var, housekeeping)].X.sum(1).A1 + 1
    )
    adata.layers["norm"] = (
        adata.layers["norm"].multiply(adata.obs["norm"].values[:, None]).tocsr()
    )
    adata.obs["Vector_10X"] = (
        adata.obs["Vector_10X"] * (adata.obs["TF"].isin(tf)).values
    )
    adata.obs["Vector_10X_norm"] = adata.obs["Vector_10X"] * adata.obs["norm"]

    if tf[0] in adata.var["symbol"].tolist():
        gene_id = adata.var.index[adata.var.symbol == tf[0]][0]
        adata.obs["Vector+endo"] = (
            adata.obs["Vector_10X"] + sc.get.obs_df(adata, gene_id).values[:, 0]
        )
        adata.obs["Vector+endo_norm"] = adata.obs["Vector+endo"] * adata.obs["norm"]

    return adata


def plot(
    mtx,
    obs_alltfs,
    var_alltfs,
    tf,
    gene,
    housekeeping,
    adata,
    datasets,
    adatamyo,
    adataadipo,
    n_top=20,
    width=2,
    citations=None,
    add_citations=True,
):
    if citations is None:
        citations = {}
    adatatf = extract_tf_dataset(mtx, obs_alltfs, var_alltfs, tf, housekeeping)
    adatatf = adatatf[adatatf.obs["Phase_corrected"] == "G1"]
    # adatatf_norm = adatatf.copy()

    # sc.pp.normalize_total(adatatf_norm)
    # sc.pp.log1p(adatatf_norm)
    # sc.pp.pca(adatatf_norm)
    # ref = adatatf_norm.obsm["X_pca"][adatatf_norm.obs["TF"] != tf].mean(0)
    # diff = adatatf_norm.obsm["X_pca"] - ref
    # adatatf_norm.obs["diff"] = np.linalg.norm(diff, axis=1)

    # targets = eyck.m.t.compare_two_groups(
    #     adatatf_norm, adatatf.obs["TF"] == tf, adatatf.obs["TF"] != tf
    # )
    # targets = targets["symbol"].head(50)
    # sc.tl.score_genes(
    #     adatatf_norm, eyck.m.t.gene_id(adatatf.var, targets), score_name="target"
    # )
    # adatatf.obs["target"] = adatatf_norm.obs["target"]
    # adatatf.obs["target"] = adatatf_norm.obs["target"] = adatatf_norm.obs["diff"]

    adatatf.obs["target"] = adatatf.obs["activity_to_D0"]

    fig = pp.grid.Figure(pp.grid.Grid(padding_height=0))

    ##
    ax = fig.main.add_under(pp.grid.Panel((width, 0.5)))
    ax.set_title(gene, fontstyle="italic")

    plotdata = pd.DataFrame(
        {
            "expression": np.log(adatatf.obs["Vector_10X_norm"] + 1),
            "dataset_cell_type": "scTF-seq",
            "target": adatatf.obs["target"],
        }
    )
    plotdata_non0 = plotdata.loc[plotdata["expression"] != 0]
    ax.scatter(
        x=plotdata_non0["expression"],
        y=plotdata_non0["target"],
        s=1 if len(plotdata_non0) > 1000 else 2,
        color="#FF6347",
        lw=0,
        rasterized=True,
    )

    extent = plotdata_non0["expression"].max()
    if extent < np.log1p(4):
        extent = np.log1p(4)
    extent_left = -extent * 0.02
    plotdata_smooth = pd.DataFrame(
        {
            "expression": np.linspace(0, extent, 100),
        }
    )
    plotdata_smooth["target"], plotdata_smooth["target_se"] = smooth_spline_fit_se(
        plotdata["expression"], plotdata["target"], plotdata_smooth["expression"]
    ).T
    ax.plot(plotdata_smooth["expression"], plotdata_smooth["target"], color="#333")
    ax.fill_between(
        plotdata_smooth["expression"],
        plotdata_smooth["target"] - plotdata_smooth["target_se"] * 2,
        plotdata_smooth["target"] + plotdata_smooth["target_se"] * 2,
        color="#333",
        alpha=0.3,
        lw=0.0,
    )
    ax.set_ylim(0, 0.7)
    ax.set_yticks([0, 0.35, 0.7])
    ax.set_yticklabels([0, 0.35, 0.7], fontsize=8)
    ax.get_yticklabels()[0].set_verticalalignment("bottom")

    ax.set_xlim(extent_left, extent)
    ax.set_ylabel(
        "Overall\ntranscriptomic\nchange",
        rotation=0,
        ha="right",
        va="center",
        fontsize=8,
    )
    ax.set_xticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if plotdata_smooth["target"].max() > 0.23:
        ix = plotdata_smooth.index[plotdata_smooth["target"] > 0.23][0]
        cut_x = plotdata_smooth["expression"][ix]
        cut_y = plotdata_smooth["target"][ix]
        ax.plot([cut_x, cut_x], [0, cut_y], color="#333", linestyle="dotted")
        ax.scatter(cut_x, cut_y, color="#333", marker="o", zorder=10)
    else:
        cut_x = 99999

    ##
    # make cellxgene data
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

    # add exogenous and endogenous
    adatatf_oi = adatatf[adatatf.obs["TF"] == gene]
    plotdata = pd.concat(
        [
            plotdata,
            pd.DataFrame(
                {
                    "expression": adatatf_oi.obs["Vector_10X_norm"],
                    "dataset_cell_type": "scTF-seq exogenous",
                },
            )
            .sample(2000, random_state=42, replace=True)
            .drop_duplicates(),
            pd.DataFrame(
                {
                    "expression": sc.get.obs_df(adata, gene, layer="norm").values[:, 0],
                    "dataset_cell_type": "scTF-seq endogenous",
                },
            )
            .sample(2000, random_state=42, replace=True)
            .drop_duplicates(),
            pd.DataFrame(
                {
                    "expression": adatatf_oi.obs["Vector+endo_norm"],
                    "dataset_cell_type": "scTF-seq exo+endo",
                },
            )
            .sample(2000, random_state=42, replace=True)
            .drop_duplicates(),
            pd.DataFrame(
                {
                    "expression": sc.get.obs_df(
                        adatamyo, eyck.m.t.gene_id(adatamyo.var, gene), layer="norm"
                    ).values[:, 0],
                    "dataset_cell_type": "Induced myogenesis",
                },
            ),
            pd.DataFrame(
                {
                    "expression": sc.get.obs_df(
                        adataadipo, eyck.m.t.gene_id(adataadipo.var, gene), layer="norm"
                    ).values[:, 0],
                    "dataset_cell_type": "Induced adipogenesis",
                },
            ),
        ]
    )

    # select top
    plotdata_mean = plotdata.groupby("dataset_cell_type")["expression"].mean()
    top = plotdata_mean.sort_values().tail(n_top).index

    # add extra if not in top
    if "scTF-seq exo+endo" not in top:
        top = ["scTF-seq exo+endo", *top]
    if "scTF-seq exogenous" not in top:
        top = ["scTF-seq exogenous", *top]
    if "scTF-seq endogenous" not in top:
        top = ["scTF-seq endogenous", *top]
    if "Induced myogenesis" not in top:
        top = ["Induced myogenesis", *top]
    if "Induced adipogenesis" not in top:
        top = ["Induced adipogenesis", *top]

    # resort
    top = plotdata_mean.loc[top].sort_values().index

    # make dataset labels
    dataset_ixs = np.unique(
        [int(x.split(" ")[-1]) for x in top if x.split(" ")[-1].isdigit()]
    )
    footnote_symbols = [*"abcdefghijklnopqrstuvwxyz"]
    footnote_symbols = [x for x in footnote_symbols if x not in citations.values()]
    datasets_oi = datasets.loc[dataset_ixs].copy()
    datasets_oi["footnote_symbol"] = [
        footnote_symbols.pop(0)
        if dataset_ix not in citations
        else citations[dataset_ix]
        for dataset_ix in dataset_ixs
    ]
    citations.update(datasets_oi["footnote_symbol"].to_dict())
    dataset_cell_type_labels = {}
    for dataset_cell_type in top:
        if dataset_cell_type.split(" ")[-1].isdigit():
            dataset_ix = int(dataset_cell_type.split(" ")[-1])

            label = dataset_cell_type.replace(
                f" {dataset_ix}",
                f"$^{datasets_oi.loc[dataset_ix]['footnote_symbol']}$",
            )
            label = label[0].lower() + label[1:]

            # if len(label) > 30:
            #     label = "…" + label[-30:]
            dataset_cell_type_labels[dataset_cell_type] = label
        else:
            dataset_cell_type_labels[dataset_cell_type] = dataset_cell_type

    # make plot data
    plotdata = plotdata.loc[plotdata["dataset_cell_type"].isin(top)]
    plotdata["dataset_cell_type"] = (
        plotdata["dataset_cell_type"]
        .astype("category")
        .cat.set_categories(top, ordered=True)
    )
    plotdata_mean = (
        plotdata.groupby("dataset_cell_type", observed=True)["expression"]
        .mean()
        .to_frame()
    )
    plotdata_mean["ix"] = range(plotdata_mean.shape[0])

    fc_palette = {
        "scTF-seq exo+endo": "#85144b",
        "scTF-seq exogenous": "#FF4136",
        "scTF-seq endogenous": "#0074D9",
        "Induced myogenesis": "#2ECC40",
        "Induced adipogenesis": "#2ECC40",
    }

    ax = fig.main.add_under(pp.grid.Panel((width, 2 if n_top == 20 else 0.12 * n_top)))

    for ix, (dataset_cell_type, data) in enumerate(
        plotdata.groupby("dataset_cell_type", observed=True)
    ):
        fc = (
            fc_palette[dataset_cell_type] + "FF"
            if dataset_cell_type in fc_palette
            else "#555"
        )
        ec = (
            fc_palette[dataset_cell_type]
            if dataset_cell_type in fc_palette
            else "#555555"
        )
        # median_color = "#333" if dataset_cell_type == "scTF-seq" else "#333"
        median_color = "#00000000"
        ax.boxplot(
            data["expression"],
            positions=[ix],
            showfliers=True,
            widths=0.95,
            patch_artist=True,
            boxprops=dict(facecolor=fc, edgecolor=ec, lw=0),
            vert=False,
            medianprops=dict(color=median_color),
            flierprops=dict(marker=".", mec="#00000000", markersize=3, mfc=ec),
            whiskerprops=dict(color=ec),
            capprops=dict(color=ec),
            whis=[5, 95],
        )
        mean = plotdata_mean.loc[dataset_cell_type, "expression"]
        ax.scatter(mean, ix, s=50, marker=".", zorder=10, color="#FFF", lw=0)
        ax.scatter(mean, ix, s=100, marker=".", zorder=9, color="#333333", lw=0)
    ax.set_yticks(plotdata_mean["ix"])
    ax.set_yticklabels(
        [dataset_cell_type_labels[x] for x in plotdata_mean.index], fontsize=6
    )
    for ticklabel in ax.get_yticklabels():
        if ticklabel.get_text() in fc_palette:
            ticklabel.set_color(fc_palette[ticklabel.get_text()])
            ticklabel.set_fontweight("bold")

    ax.tick_params(axis="y", which="major", pad=1, length=0)

    ax.set_xticks(np.log1p([0, 0.5, 1, 2, 5, 10, 20, 50]))
    ax.set_xticklabels([0, "½", 1, 2, 5, 10, 20, 50])
    ax.set_xlim(extent_left, extent)

    ax.set_xlabel("log1p TF / HK expression", fontsize=8)

    ax.axvline(cut_x, color="#000", linestyle="dotted", zorder=10, lw=1)

    plotdata = datasets_oi.groupby("collection_doi_label").agg(
        {"footnote_symbol": tuple}
    )
    plotdata = plotdata.sort_values("footnote_symbol", ascending=False)
    if add_citations:
        ## labels
        plotdata["ix"] = range(plotdata.shape[0])
        ax = fig.main.add_under(pp.grid.Panel((2, 0.08 * len(plotdata))))
        fig.main.paddings_height[1] = 0.5
        ax.axis("off")

        for label, dataset in plotdata.iterrows():
            ax.text(
                0,
                dataset["ix"],
                f"$^{{{','.join(dataset['footnote_symbol'])}}}$ {label}",
                ha="left",
                va="center",
                fontsize=5,
            )
        ax.set_ylim(-0.5, len(plotdata) - 0.5)

    return fig
