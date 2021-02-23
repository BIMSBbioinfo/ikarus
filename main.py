import yaml
import argparse
from pathlib import Path
import pandas as pd
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

    label_upregs_list = config['gene_selector']['label_upregs']
    label_downregs_list = config['gene_selector']['label_downregs']
    for label_upreg, label_downreg in zip(label_upregs_list, label_downregs_list):
        list_of_DE_results_df = [
            gene_selector(
                adatas[name],
                config['gene_selector']['obs_name'],
                label_upreg,
                label_downreg,
                config['gene_selector']['lfc_threshold'],
                config['gene_selector']['pval_threshold'],
                config['gene_selector']['DE_method'],
                config['gene_selector']['sort_by'],
                config['gene_selector']['sort_ascending']
            ) for name in config['gene_selector']['names']
        ]
        list_of_DE_results_df = [d for d in list_of_DE_results_df if d is not None]


        # gene list integration
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
        if label_downreg == None or not config['run']['gene_selector']:
            write_label_downreg = 'all'
        else:
            write_label_downreg = label_downreg
        (Path.cwd() / 'out').mkdir(parents=True, exist_ok=True)
        pd.DataFrame(gene_list).to_csv(
            (
                f"out/{label_upreg}"
                f"_vs_{write_label_downreg}_gene_list.csv"
            ),
            header=None,
            index=None
        )
        DE_results_df.to_csv(
            (
                f"out/{label_upreg}"
                f"_vs_{write_label_downreg}_DE_results.csv"
            ),
            index=None
        )


# cell scoring
if config['run']['cell_scorer']:
    paths = config['cell_scorer']['paths']
    names = config['cell_scorer']['names']
    obs_names = config['cell_scorer']['obs_names']
    for path, name, obs_name in zip(paths, names, obs_names):
        adata = load_adata(
            path,
            config['cell_scorer']['adata_is_given'],
            config['cell_scorer']['sparse_is_given']

        )
        adata = (
            adata if config['cell_scorer']['is_preprocessed'] else preprocess_adata(adata)
        )

        if not config['cell_scorer']['adata_is_given'] and not config['cell_scorer']['is_preprocessed']:
            adata.write(path + 'adata.h5ad')
        # use first n cells for test reason
        adata = (
            adata[:config['cell_scorer']['n_cells']]
        )

        if config['run']['gene_selector']:
            gene_list_dict = {}
            for label_upreg, label_downreg in zip(config['gene_selector']['label_upregs'], config['gene_selector']['label_downregs']):
                if label_downreg == None:
                    write_label_downreg = 'all'
                else:
                    write_label_downreg = label_downreg
                gene_list_dict[label_upreg] = pd.read_csv(f"out/{label_upreg}_vs_{write_label_downreg}_gene_list.csv", header=None).values.ravel().tolist()

        else:
            gene_lists = [
                pd.read_csv(gl, header=None).values.ravel().tolist() for gl in config['cell_scorer']['gene_lists']
            ]
            label_upregs = config['cell_scorer']['label_upregs']
            gene_list_dict = dict(zip(label_upregs, gene_lists))

        (Path.cwd() / 'out' / f"{name}").mkdir(parents=True, exist_ok=True)
        cell_scorer(
            adata,
            name,
            gene_list_dict,
            singscore_fun,
            obs_name
            )


# cell annotation
if config['run']['cell_annotator']:
    paths = config['cell_annotator']['paths']
    names = config['cell_annotator']['names']
    obs_names = config['cell_annotator']['obs_names']
    test_name = config['cell_annotator']['test_name']
    training_names = config['cell_annotator']['training_names']

    results = {}
    adatas = {}
    connectivities = {}

    for name, path in zip(names, paths):
        # adata
        adatas[name] = load_adata(path, adata_is_given=True)

        # scores
        results[name] = load_scores(name, adatas[name])
        
    # connectivities
    if config['cell_annotator']['connectivities_given']:
        connectivities[test_name] = load_connectivities(test_name)
    else:
        connectivities[test_name] = calculate_connectivities(
            adatas[test_name], 
            n_neighbors=100,
            use_highly_variable=False
            )

    # label propagation
    cell_annotator(
        connectivities, 
        results,
        names,
        obs_names,
        training_names,
        test_name,
        config['cell_annotator']['certainty_threshold'],
        config['cell_annotator']['n_iter']
        )