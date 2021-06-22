import scipy
import pandas as pd
import scanpy as sc
from anndata import AnnData, read_h5ad
from scipy.sparse import issparse, csr_matrix


def load_h5ad(adata_h5ad):
    adata = read_h5ad(adata_h5ad)
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    adata.var_names_make_unique()
    return adata


def load_sparse(gex_npz, gene_symbols, cell_labels):
    cl = pd.read_csv(cell_labels)
    genes = pd.read_csv(gene_symbols, header=None, names=["gene_symbol"])
    genes.index = genes["gene_symbol"].values
    sparse = scipy.sparse.load_npz(gex_npz)
    adata = AnnData(sparse, var=genes, obs=cl)
    adata.var_names_make_unique()
    return adata


def load_dense(gex_csv, gene_symbols, cell_labels):
    cl = pd.read_csv(cell_labels)
    genes = pd.read_csv(gene_symbols, header=None, names=["gene_symbol"])
    genes.index = genes["gene_symbol"].values
    dense = pd.read_csv(gex_csv, index_col=0)
    sparse = csr_matrix(dense)
    adata = AnnData(sparse, var=genes, obs=cl)
    adata.var_names_make_unique()
    return adata


def preprocess_adata(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata