# GeneSigNet

A computational method for constructing Gene - Signature causal network from gene expression and exposure of mutational signatures. 
#### To install the GeneSigNet package

```
pip install GeneSigNet
```
Installation source: [GeneSigNet 0.1.0](https://pypi.org/project/GeneSigNet/)

#### Pre-required installations

Python libraries pandas, numpy, scipy, sklearn and warnings are required to be installed before running GeneSigNet, and pandas and pyvis are required to be installed before running the network visualization module VisualizeNetwork.

#### Sample python script to load data and run GeneSigNet

```
import pandas as pd
import GeneSigNet as GSN

# A gene expression data matrix (columns are genes and rows are samples) 
ExpData = pd.read_csv('Gene_Expression_Data.csv', index_col=0) 

# A signature exposure datamatrix (columns are signatures and rows are samples)
SigData = pd.read_csv('Signature_Exposure_Data.csv', index_col=0)  

maxit=10000     # maximum number of iterations for selecting sparce partial correlation
tolerance=1e-12 # tolerance of iteration for selecting sparce partial correlation

D=pd.concat([ExpData, SigData], axis=1)
Net=GSN.WeightMatrix(D, maxit, tolerance)  
Weight_Matrix=Net.ConstructNet()
Weight_Matrix.to_csv('Weight_Matrix.csv')
```
#### Sample python script to visualize the resulting network
```
import pandas as pd
from pyvis.network import Network

Signatures=SigData.columns # or Signatures=['SBS1', 'SBS2',...,] (names of signature nodes)
th=0.05                    # The threshold parameter for selecting edges to be included in the visualizing network
VisualizeNetwork(Weight_Matrix, Signatures, th)
```

#### Python [Package](package) 
* The GeneSigNet method is implemented in python and the codes are available as [Python module](package/GeneSigNet.py) and [Jupyter Notebook module](package/GeneSigNet.ipynb).

* The python scripts for running GeneSigNet are available as [Python script](package/Call_GeneSigNet.py) and [Jupyter Notebook script](package/Call_GeneSigNet.ipynb). The scripts are recommended to load data and run GeneSigNet.

* Simalated data matrices, [Gene Expression](package/Gene_Expression_Simulated_Data.cvs) and [Exposure of Mutational Signatures](package/Signature_Exposure_Simulated_Data.cvs), are provided as example inputs to run the package. 

* Python module [VisualizeNetwork](package/VisualizeNetwork.ipynb) visualize the subnetwork covering hub-nodes (Signatures) and their up and downstream nodes (causal and affected genes). The weight matrix inferred by GeneSigNet method is as input of VisualizeNetwork module.
  
#### Results of the analysis on breast and lung cancer data.
The following files proves the results for the analysis on breast and lung cancer data 
- [Results](results): As result of the analysis in the two cancer data sets, the weight matrices inferred to represent the directed interactions among genes and signatures, and subnetwork figures representing the interaction between signatures and their up and down stream genes.   

#### Cancer Data (.cvs format)

The following files provides the gene expression and exposure of mutational signatures for cancer patients.

* [Gene_Expression_BRCA.csv](data/Gene_Expression_BRCA.csv): The normalized gene expression data (ICGC data portal) for 266 breast cancer (BRCA) patients. Gene expression profiles for 2,204 genes involved in either DNA metabolic or immune response processes of the Gene Ontology (GO) database were
selected for the analysis.
* [Signature_Exposure_BRCA.csv](data/Signature_Exposure_BRCA.csv): Exposure of 6 mutational signatures and HRD mutational status in breast cancer are used for the analysis
* [Gene_Expression_LUAD.csv](data/Gene_Expression_LUAD.csv): The normalized gene expression data (TCGA data portal) for 466 lung cancer (LUAD) patients. Gene expression profiles for 2,433 genes involved in either DNA metabolic or immune response processes of the Gene Ontology (GO) database were
selected for the analysis.
* [Signature_Exposure_LUAD.csv](data/Signature_Exposure_LUAD.csv): Exposure of 6 mutational signatures known to be operative in lung cancer are used for the analysis.
