import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from ikarus import utils


def select_genes(
    adata,
    obs_name,
    label_upreg,
    label_downreg=None,
    lfc_threshold=1,
    pval_threshold=0.1,
    DE_method="t-test_overestim_var",
    sort_by="logfoldchanges",
    sort_ascending=False,
):
    """Method selects genes for further consideration.
    
    t-test with overestimated variance is used to compute an 
    approximation of log 2 fold changes between two cell groups.
    Select genes coping with thresholds for log 2 fold changes 
    and p-val. 
    """
    
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
            print(
                f"Either {label_upreg} or {label_downreg} is not available in given obs_name. "
                f"None is returned."
            )
            return None

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
    """Method integrates pair-wise signatures of different datasets.
    
    For each pair-wise signature compute the weighted average of 
    log 2 fold changes of genes across different datasets. Sort 
    genes accordingly and pick top_x genes.
    """
    
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


def create_all(
    label_upregs_list,
    label_downregs_list,
    adatas_dict,
    names_list,
    obs_names_list,
    integration_fun=utils.intersection_fun,
    top_x=300
):
    """Method creates multiple pair-wise gene list signatures.
    
    Loop over all pair-wise comparisons and over given input datasets,
    and perform t-tests with overestimated variance so that one obtains 
    for each gene n (one for each dataset) estimates of log2fold changes.
    For integrating the outcome of all datasets, it is considered the 
    intersection or union of genes available in different datasets. Then 
    for each gene the weighted average of log2fold changes is computed. 
    For each pair-wise comparison it is taken the top_x of upregulated 
    (highest averaged log2fold changes) genes.
    
    Parameters
    ----------
    label_upregs_list : list of str
        label which should be considered up-regulated. Should be available 
        in the given obs column of the anndata object.
    label_downregs_list : list of str
        label which should be considered down-regulated. Should be available 
        in the given obs column of the anndata object.
    adatas_dict : dict {"name": anndata_object}
        given anndata objects (dict value) with the corresponding names
        (dict key).
    names_list : list of str
        names of anndata objects in accordance to adatas_dict keys.
    obs_names_list : list of str
        for each anndata object in adatas_dict the corresponding obs column 
        where to search for up- or down-regulated labels.
    integration_fun : function
        function to integrate target genes from different datasets. Supported
        "utils.intersection_fun" and "utils.union_fun".
    top_x : int
        for each pair-wise comparison take the top_x genes with highest 
        averaged (across datasets) log2fold changes.
        
    Returns
    -------
    dict : {"{label_upreg}_vs_{label_downreg}": signature}
        Returns dict of all pair-wise signatures.
    """
    
    signatures = {}
    for label_upreg, label_downreg in zip(label_upregs_list, label_downregs_list):
        DE_dfs_list = []
        for name, obs_name in zip(names_list, obs_names_list):
            DE_dfs_list.append(
                select_genes(
                    adatas_dict[name], 
                    obs_name,
                    label_upreg, 
                    label_downreg
                )
            )
        # If either label_upreg or label_downreg is not available in the provided data set
        # 'select_genes' returns None 
        # before continuing we have to remove all the Nones from the DE dataframe list
        DE_dfs_list = [d for d in DE_dfs_list if d is not None]
        signatures[f"{label_upreg}_vs_{label_downreg}"], _ = integrate(
            DE_dfs_list, 
            integration_fun=integration_fun,
            top_x=top_x
        )
    return signatures


def save_gmt(gene_lists, gene_list_names, out_dir):
    gmt = pd.DataFrame(gene_lists, index=gene_list_names)
    gmt.insert(0, "00", "ikarus")
    path = Path.cwd() / out_dir
    path.mkdir(parents=True, exist_ok=True)
    gmt.to_csv(path / "signatures.gmt", header=None, sep="\t")
