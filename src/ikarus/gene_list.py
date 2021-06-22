import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path


def select_genes(
    adata,
    obs_name,
    label_upreg,
    label_downreg=None,
    lfc_threshold=3,
    pval_threshold=0.1,
    DE_method="t-test_overestim_var",
    sort_by="logfoldchanges",
    sort_ascending=False,
):
    if label_downreg is None:
        adata = adata.copy()
        adata.obs["1_vs_all"] = [
            label_upreg if label == label_upreg else "Other"
            for label in adata.obs[obs_name]
        ]
        obs_name = "1_vs_all"
        label_downreg = "Other"
    else:
        unique_labels = np.unique(adata.obs[obs_name])
        if label_upreg not in unique_labels or label_downreg not in unique_labels:
            raise IndexError(
                f"Either {label_upreg} or {label_downreg} is not available in given obs_name."
            )

    adata = sc.tl.rank_genes_groups(
        adata,
        groupby=obs_name,
        groups=[label_upreg, label_downreg],
        key_added=f"{DE_method}_results",
        n_genes=15000,
        copy=True,
        method=DE_method,
    )

    DE_df = pd.DataFrame()
    DE_df["gene_symbol"] = adata.uns[f"{DE_method}_results"]["names"][label_upreg]
    DE_df["logfoldchanges"] = adata.uns[f"{DE_method}_results"]["logfoldchanges"][
        label_upreg
    ]
    DE_df["p"] = adata.uns[f"{DE_method}_results"]["pvals"][label_upreg]
    DE_df["padj"] = adata.uns[f"{DE_method}_results"]["pvals_adj"][label_upreg]
    DE_df["scores"] = adata.uns[f"{DE_method}_results"]["scores"][label_upreg]
    DE_df = DE_df.loc[
        ((DE_df["padj"] < pval_threshold) & (DE_df["logfoldchanges"] > lfc_threshold))
    ]
    DE_df.sort_values(by=sort_by, ascending=sort_ascending, inplace=True)
    return DE_df


def integrate(
    DE_dfs_list,
    integration_fun,
    integrate_by="logfoldchanges",
    sort_ascending=False,
    top_x=100,
):
    dfs = [df.copy() for df in DE_dfs_list]
    if len(dfs) == 0:
        raise RuntimeError(
            "Neither input dataset contains either upregulated or downregulated labels."
        )
    for i, df in enumerate(dfs):
        dfs[i].set_index(df["gene_symbol"], inplace=True)
        dfs[i] = df[integrate_by]
        dfs[i].name = f"{integrate_by}{i}"
    DE_df = integration_fun(dfs)

    for i in range(len(dfs)):
        DE_df[f"{integrate_by}{i}"] /= DE_df[f"{integrate_by}{i}"].max()
    DE_df["weighted_avg"] = DE_df[[f"{integrate_by}{i}" for i in range(len(dfs))]].mean(
        axis=1
    )
    DE_df.sort_values(by="weighted_avg", ascending=sort_ascending, inplace=True)
    DE_df[integrate_by] = DE_df["weighted_avg"]
    DE_df["gene_symbol"] = DE_df.index.values

    gene_list = list(DE_df.index.values)
    gene_list = gene_list[: int(top_x)] if len(gene_list) >= top_x else gene_list
    return gene_list, DE_df


def save_gmt(gene_lists, gene_list_names, out_dir):
    gmt = pd.DataFrame(gene_lists, index=gene_list_names)
    gmt.insert(0, "00", "ikarus")
    path = Path.cwd() / out_dir
    path.mkdir(parents=True, exist_ok=True)
    gmt.to_csv(path / "signatures.gmt", header=None, sep="\t")
