import scipy
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData, read_h5ad
from sklearn.linear_model import LogisticRegression
from PySingscore.singscore import singscore


def preprocess_adata(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def load_adata(path, adata_is_given=False):
    if adata_is_given:
        adata = read_h5ad(f'{path}adata.h5ad')
    else:
        cl = pd.read_csv(f'{path}cell_labels.csv')
        genes = pd.read_csv(f'{path}genes_symbol.csv', header=None, names=['gene_symbol'])
        genes.index = genes['gene_symbol'].values
        sparse = scipy.sparse.load_npz(f'{path}matrix_sparse.npz')
        adata = AnnData(sparse, var=genes, obs=cl)
    adata.var_names_make_unique()
    return adata


def gene_selector(
    adata,
    obs_name,
    label_upreg,
    label_downreg=None,
    lfc_threshold=3,
    pval_threshold=0.1,
    DE_method='t-test_overestim_var',
    sort_by='logfoldchanges',
    sort_ascending=False
):
    if label_downreg is None:
        adata = adata.copy()
        adata.obs['1_vs_all'] = [label_upreg if label == label_upreg
                                 else 'Other' for label in adata.obs[obs_name]]
        obs_name = '1_vs_all'
        label_downreg = 'Other'

    if label_downreg is not None:
        unique_labels = np.unique(adata.obs[obs_name])
        if label_upreg not in unique_labels or label_downreg not in unique_labels:
            return None

    adata = sc.tl.rank_genes_groups(
        adata,
        groupby=obs_name,
        groups=[label_upreg, label_downreg],
        key_added=f'{DE_method}_results',
        n_genes=15000,
        copy=True,
        method=DE_method
    )

    DE_results_df = pd.DataFrame()
    DE_results_df['gene_symbol'] = adata.uns[f'{DE_method}_results']['names'][label_upreg]
    DE_results_df['logfoldchanges'] = adata.uns[f'{DE_method}_results']['logfoldchanges'][label_upreg]
    DE_results_df['p'] = adata.uns[f'{DE_method}_results']['pvals'][label_upreg]
    DE_results_df['padj'] = adata.uns[f'{DE_method}_results']['pvals_adj'][label_upreg]
    DE_results_df['scores'] = adata.uns[f'{DE_method}_results']['scores'][label_upreg]
    DE_results_df = DE_results_df.loc[((DE_results_df['padj'] < pval_threshold)
                                       & (DE_results_df['logfoldchanges'] > lfc_threshold))]
    DE_results_df.sort_values(by=sort_by, ascending=sort_ascending, inplace=True)
    # gene_list = list(DE_results_df['gene_symbol'].values)
    return DE_results_df


def gene_list_integrator(
    list_of_DE_results_df,
    integration_fun,
    integrate_by='logfoldchanges',
    sort_ascending=False,
    top_x=100
):
    dfs = [df.copy() for df in list_of_DE_results_df]
    if len(dfs) == 0:
        raise('Error: Neither input dataset contains either upregulated or downregulated labels.')
    for i, df in enumerate(dfs):
        dfs[i].set_index(df['gene_symbol'], inplace=True)
        dfs[i] = df[integrate_by]
        dfs[i].name = f'{integrate_by}{i}'
    DE_results_df = integration_fun(dfs)

    for i in range(len(dfs)):
        DE_results_df[f'{integrate_by}{i}'] /= DE_results_df[f'{integrate_by}{i}'].max()
    DE_results_df['weighted_avg'] = (
        DE_results_df[[
            f'{integrate_by}{i}' for i in range(len(dfs))
        ]].mean(axis=1)
    )
    DE_results_df.sort_values(by='weighted_avg', ascending=sort_ascending, inplace=True)
    DE_results_df[integrate_by] = DE_results_df['weighted_avg']
    DE_results_df['gene_symbol'] = DE_results_df.index.values

    gene_list = list(DE_results_df.index.values)
    gene_list = gene_list[:int(top_x)] if len(gene_list) >= top_x else gene_list
    return gene_list, DE_results_df


def intersection_fun(x): return pd.concat(x, axis=1, join='inner')
def union_fun(x): return pd.concat(x, axis=1, join='outer')


def cell_scorer(
    adata,
    name,
    gene_list_dict,
    scoring_fun,
    obs_name
):
    df = adata.to_df()
    df = df.T
    for label_upreg, gene_list in gene_list_dict.items():
        all_scores = pd.DataFrame(index=adata.obs[obs_name].values)
        all_scaled_scores = pd.DataFrame(index=adata.obs[obs_name].values)
        scores = scoring_fun(gene_list, df)
        scores = scores['total_score'].values
        all_scores[label_upreg] = scores
        all_scaled_scores[label_upreg] = ((scores - scores.min())
                                          / (scores.max() - scores.min()))
        all_scores.to_csv(
            f'out/{name}/{label_upreg}_score.csv',
            index_label='cells'
        )
        all_scaled_scores.to_csv(
            f'out/{name}/{label_upreg}_scaled_score.csv',
            index_label='cells'
        )
        adata.obs[f'{label_upreg}_score'] = all_scores[label_upreg].values
        adata.obs[f'{label_upreg}_scaled_score'] = all_scaled_scores[label_upreg].values

    # adata_out.write(f"out/scored_adata.h5ad")

    # return all_scores, all_scaled_scores, adata


def singscore_fun(gene_list, df):
    return singscore.score(
        up_gene=gene_list,
        sample=df,
        down_gene=False,
        norm_method='standard',
        full_data=False
    )


def cell_annotator(
    connectivities_dict,
    results_dict,
    names_list,
    training_names_list,
    test_name,
    certainty_threshold,
    n_iter
):
    #certain?
    for name in names_list:
        results_dict[name]['certain'] = False
        results_dict[name].loc[
            (abs(results_dict[name]['Normal'] - results_dict[name]['Tumor']) > certainty_threshold),
            'certain'
            ] = True

    # first training, pre label propagation
    # Tumor vs Normal classification
    X_features = ['Normal', 'Tumor', 'Epithelial']
    X_dict = {}
    y_dict = {}
    Model = LogisticRegression()
    for name in names_list:
        X_dict[name] = results_dict[name].loc[:, X_features]
        y_dict[name] = results_dict[name].loc[:, 'tier_0']
    
    y_train = pd.concat(
        [y_dict[name] for name in training_names_list], 
        axis=0, ignore_index=True)
    X_train = pd.concat(
        [X_dict[name] for name in training_names_list], 
        axis=0, ignore_index=True)
    X_test = X_dict[test_name]
    
    _ = Model.fit(X_train, y_train)
    y_pred_lr = Model.predict(X_test)
    y_pred_proba_lr = Model.predict_proba(X_test)

    results_dict[test_name]['LR_proba_Normal'] = y_pred_proba_lr[:, 0]
    results_dict[test_name]['LR_proba_Tumor'] = y_pred_proba_lr[:, 1]
    results_dict[test_name]['tier_0_prediction_LR'] = y_pred_lr

    proba_tier_0 = results_dict[test_name].loc[:, ['LR_proba_Normal', 'LR_proba_Tumor']].copy()
    proba_tier_0.columns = ['Normal', 'Tumor']
    proba_tier_0.loc[results_dict[test_name]['certain'] == False] = 0

    for _ in range(n_iter):
        proba_tier_0 = proba_tier_0.apply(
            lambda x: np.dot(x, connectivities_dict[test_name]
            ), axis=0, raw=True)

    # get prediction and arrange output file
    # q&d: find better solution to this mess
    results_dict[test_name][
        'tier_0_prediction_LR_with_label_propagation'
        ] = results_dict[test_name]['tier_0_prediction_LR'].copy()
    results_dict[test_name][
        'tier_0_prediction_LR_with_label_propagation'
        ] = proba_tier_0.idxmax(axis=1)

    results_dict[test_name]['Normal_lr_proba_ns'] = proba_tier_0['Normal']
    results_dict[test_name]['Tumor_lr_proba_ns'] = proba_tier_0['Tumor']

    lr_proba_smoothed = results_dict[test_name].loc[:, [
        'Normal_lr_proba_ns', 
        'Tumor_lr_proba_ns'
        ]]
    lr_proba_smoothed = lr_proba_smoothed.divide(lr_proba_smoothed.sum(axis=1), axis=0)
    
    results_dict[test_name] = results_dict[test_name].loc[:, [
        'raw', 
        'tier_0', 
        'Normal', 
        'Tumor', 
        'Epithelial', 
        'tier_0_prediction_LR', 
        'LR_proba_Normal', 
        'LR_proba_Tumor', 
        'tier_0_prediction_LR_with_label_propagation'
        ]]
    results_dict[test_name].columns = [
        'raw', 
        'tier_0', 
        'Singscore_Normal', 
        'Singscore_Tumor', 
        'Singscore_Epithelial', 
        'tier_0_prediction_LR', 
        'LR_proba_Normal', 
        'LR_proba_Tumor', 
        'tier_0_prediction_LR_with_label_propagation'
        ]
    results_dict[test_name][
        'LR_with_label_propagation_proba_Normal'
        ] = lr_proba_smoothed.iloc[:, 0]
    results_dict[test_name][
        'LR_with_label_propagation_proba_Tumor'
        ] = lr_proba_smoothed.iloc[:, 1]
    
    results_dict[test_name].to_csv(f"/out/{test_name}/results_final.csv")



def load_scores(
    name,
    adata
):
    if "tier_0_hallmark_corrected" in adata.obs.columns:
        adata.obs["tier_0_raw"] = adata.obs["tier_0"]
        adata.obs["tier_0"] = adata.obs["tier_0_hallmark_corrected"]

    scores = pd.DataFrame()
    fnames = [
        f"/out/{name}/Normal",
        f"/out/{name}/Tumor",
        f"/out/{name}/Epithelial",
    ]
    for fname in fnames:
        s = pd.read_csv(f"{fname}_scaled_score.csv", index_col=False)
        s.drop(["cells"], axis=1, inplace=True)
        scores = pd.concat([scores, s], axis=1)

    result_df = pd.concat(
        [adata.obs.reset_index(drop=True), scores], 
        axis=1
    )
    tier_0 = result_df[["Normal", "Tumor"]]
    result_df["tier_0_pred"] = tier_0.idxmax(axis=1)

    return result_df


def calculate_connectivities(
    adata,
    n_neighbors,
    use_highly_variable
):
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata, use_highly_variable=use_highly_variable)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, method='umap')
    connectivities= adata.uns['neighbors']['connectivities'].todense()
    return connectivities


def load_connectivities(
    name
):
    connectivities = scipy.sparse.load_npz(
        f'/out/{name}/connectivities_sparse.npz'
        ).todense()
    return connectivities