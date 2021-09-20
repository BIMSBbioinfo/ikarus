<!--Get the samples from https://www.adobe.com/go/pdfembedapi_samples-->
<!DOCTYPE html>
<html>
<head>
 <title>Adobe Document Services PDF Embed API Sample</title>
 <meta charset="utf-8"/>
 <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1"/>
 <meta id="viewport" name="viewport" content="width=device-width, initial-scale=1"/>
</head>
<body style="margin: 0px">
 <div id="adobe-dc-view"></div>
 <script src="https://documentcloud.adobe.com/view-sdk/main.js"></script>
 <script type="text/javascript">
    document.addEventListener("adobe_dc_view_sdk.ready", function()
    {
        var adobeDCView = new AdobeDC.View({clientId: "<YOUR_CLIENT_ID>", divId: "adobe-dc-view"});
        adobeDCView.previewFile(
       {
          content:   {location: {url: "https://documentcloud.adobe.com/view-sdk-demo/PDFs/Bodea Brochure.pdf"}},
          metaData: {fileName: "Bodea Brochure.pdf"}
       });
    });
 </script>
</body>
</html>


# Ikarus
Ikarus is a stepwise machine learning pipeline that tries to cope with a task of distinguishing tumor cells from normal cells. Leveraging multiple annotated single cell datasets it can be used to define a gene set specific to tumor cells. First, the latter gene set is used to rank cells and then to train a logistic classifier for the robust classification of tumor and normal cells. Finally, sensitivity is increased by propagating the cell labels based on a custom cell-cell network. Ikarus is tested on multiple single cell datasets to ascertain that it achieves high sensitivity and specificity in multiple experimental contexts. 

![Chema](ikarus_schema.pdf)

## Installation
Make sure you are using python >= 3.8 before installing ikarus. If that requirement is fulfilled, ikarus can be installed from a gitthub repo:
```
git clone https://github.com/BIMSBbioinfo/ikarus.git
cd ikarus
pip install -e .
```
Alterantively, one can install ikarus' master branch directly from github:
```
python -m pip install git+https://github.com/BIMSBbioinfo/ikarus.git
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
