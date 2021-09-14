# Ikarus
Ikarus is a stepwise machine learning pipeline that tries to cope with a task of distinguishing tumor cells from normal cells. Leveraging multiple annotated single cell datasets it can be used to define a gene set specific to tumor cells. First, the latter gene set is used to rank cells and then to train a logistic classifier for the robust classification of tumor and normal cells. Finally, sensitivity is increased by propagating the cell labels based on a custom cell-cell network. Ikarus is tested on multiple single cell datasets to ascertain that it achieves high sensitivity and specificity in multiple experimental contexts. 


## Installation
Ikarus can be installed by cloning the repo, navigating into the directory and running:
```
pip install -e .
```


## Usage
The easiest option to get started is to use the provided Tumor/Normal gene lists and the pretrained model.
```
from ikarus import classifier

model = classifier.Ikarus(signatures_gmt=signatures_path)
model.load_core_model(model_path)
predictions = model.predict(test_adata, 'test_name')
```

More information on how to train a model or how to create own gene lists is provided in the [tutorial notebook](tutorial.ipynb).
