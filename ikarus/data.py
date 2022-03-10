import scipy
import pandas as pd
import scanpy as sc
from anndata import AnnData, read_h5ad
from scipy.sparse import issparse, csr_matrix


def load_h5ad(adata_h5ad):
    """Helper function to load an anndata object.
    
    Please see 
    "https://anndata.readthedocs.io/en/latest/api.html#reading" 
    for further information.
    """
    
    adata = read_h5ad(adata_h5ad)
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    adata.var_names_make_unique()
    return adata


def load_sparse(gex_npz, gene_symbols, cell_labels):
    """Helper function to create an anndata object.
    
    Load sparse expression matrix, genes (symbols) and 
    annotation and create an anndata object. The input 
    format is very specific. Please see
    "https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html" 
    for further information.
    """
        
    cl = pd.read_csv(cell_labels)
    genes = pd.read_csv(gene_symbols, header=None, names=["gene_symbol"])
    genes.index = genes["gene_symbol"].values
    sparse = scipy.sparse.load_npz(gex_npz)
    adata = AnnData(sparse, var=genes, obs=cl)
    adata.var_names_make_unique()
    return adata


def load_dense(gex_csv, gene_symbols, cell_labels):
    """Helper function to create an anndata object.
    
    Load dense expression matrix, genes (symbols) and 
    annotation and create an anndata object with sparse 
    values. The input format is very specific. Please see
    "https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html" 
    for further information.
    """    
    
    cl = pd.read_csv(cell_labels)
    genes = pd.read_csv(gene_symbols, header=None, names=["gene_symbol"])
    genes.index = genes["gene_symbol"].values
    dense = pd.read_csv(gex_csv, index_col=0)
    sparse = csr_matrix(dense)
    adata = AnnData(sparse, var=genes, obs=cl)
    adata.var_names_make_unique()
    return adata


def preprocess_adata(adata):
    """Helper function to preprocess an anndata object.
    
    Basic preprocessing pipeline based on scanpy. Perform 
    normalization of counts per cell and logarithmize the 
    expression. Feel free to adapt the preprocessing. Please 
    see
    "https://scanpy.readthedocs.io/en/stable/api.html#module-scanpy.pp"
    for further information.
    """    
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata