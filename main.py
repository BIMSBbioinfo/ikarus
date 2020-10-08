import yaml
import argparse
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
from PySingscore.singscore import singscore
from utils import *


parser = argparse.ArgumentParser(
    description='Ikarus pipeline',
    formatter_class=argparse.RawTextHelpFormatter
)
parser.format_help()
parser.add_argument(
    '--config_fname',
    '-c',
    metavar='',
    type=str,
    default='config.yaml',
    help=
    'config_fname corresponds to the file name \n'
    'of the used .yaml configuration file. \n'
    '(default: config.yaml)'
)
args = parser.parse_args()

with open(args.config_fname) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


# gene selection
if config['run']['gene_selector']:
    adatas = {}
    for p, n, a, preproc in zip(
        config['gene_selector']['paths'],
        config['gene_selector']['names'],
        config['gene_selector']['adata_is_given'],
        config['gene_selector']['is_preprocessed']
    ):
        adata = load_adata(p, a)
        adatas[n] = adata if preproc else preprocess_adata(adata)
        if not a and not preproc:
            adatas[n].write(p + 'adata.h5ad')

    list_of_DE_results_df = [
        gene_selector(
            adatas[name],
            config['gene_selector']['obs_name'],
            config['gene_selector']['label_upreg'],
            config['gene_selector']['label_downreg'],
            config['gene_selector']['lfc_threshold'],
            config['gene_selector']['pval_threshold'],
            config['gene_selector']['DE_method'],
            config['gene_selector']['sort_by'],
            config['gene_selector']['sort_ascending']
        ) for name in config['gene_selector']['names']
    ]
    list_of_DE_results_df = [d for d in list_of_DE_results_df if d is not None]


# gene list integration
if config['run']['gene_list_integrator']:
    if config['gene_list_integrator']['integration_fun'] == 'intersection':
        integration_fun = intersection_fun
    elif config['gene_list_integrator']['integration_fun'] == 'union':
        integration_fun = union_fun
    else:
        raise("Error: given integration currently not provided. Use 'intersection' or 'union' instead.")


    if not config['run']['gene_selector']:
        list_of_DE_results_df = [
            pd.read_csv(fname) for fname in config['gene_list_integrator']['DE_results_dfs']
        ]
    gene_list, DE_results_df = gene_list_integrator(
        list_of_DE_results_df,
        integration_fun,
        config['gene_list_integrator']['integrate_by'],
        config['gene_list_integrator']['sort_ascending'],
        config['gene_list_integrator']['top_x']
    )
    if config['gene_selector']['label_downreg'] == None or not config['run']['gene_selector']:
        write_label_downreg = 'all'
    else:
        write_label_downreg = config['gene_selector']['label_downreg']
    pd.DataFrame(gene_list).to_csv(
        (
            f"out/{config['gene_selector']['label_upreg']}"
            f"_vs_{write_label_downreg}_gene_list.csv"
        ),
        header=None,
        index=None
    )
    DE_results_df.to_csv(
        (
            f"out/{config['gene_selector']['label_upreg']}"
            f"_vs_{write_label_downreg}_DE_results.csv"
        ),
        index=None
    )


# cell annotation
if config['run']['cell_annotator']:
    adatas = {}
    adata = load_adata(
        config['cell_annotator']['path'],
        config['cell_annotator']['adata_is_given']
    )
    adatas[config['cell_annotator']['name']] = (
        adata if config['cell_annotator']['is_preprocessed'] else preprocess_adata(adata)
    )

    if not config['cell_annotator']['adata_is_given'] and not config['cell_annotator']['is_preprocessed']:
        adatas[config['cell_annotator']['name']].write(config['cell_annotator']['path'] + 'adata.h5ad')
    # use first n cells for test reason
    adatas[config['cell_annotator']['name']] = (
        adatas[config['cell_annotator']['name']][:config['cell_annotator']['n_cells']]
    )

    if config['run']['gene_list_integrator']:
        gene_list_dict = {config['gene_selector']['label_upreg']: gene_list}
    else:
        gene_lists = [
            pd.read_csv(gl, header=None).values.ravel().tolist() for gl in config['cell_annotator']['gene_lists']
        ]
        label_upregs = config['cell_annotator']['label_upregs']
        gene_list_dict = dict(zip(label_upregs, gene_lists))

    cell_annotator(
        adatas[config['cell_annotator']['name']],
        gene_list_dict,
        singscore_fun,
        config['cell_annotator']['obs_name'])

