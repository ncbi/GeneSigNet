# GeneSigNet
A computational method for constructing Gene - Signature causal network from gene expression and exposure of mutational signatures. 
#### Source install

```
pip install GeneSigNet
```

#### Sample python script to load data and run GeneSigNet

```
import pandas as pd
import GeneSigNet as GSN

ExpData = pd.read_csv('Gene_Expression_Simulated_Data.csv', index_col=0)  
SigData = pd.read_csv('Signature_Exposure_Simulated_Data.csv', index_col=0)

maxit=10000; tolerance=1e-12; 
D=pd.concat([ExpData, SigData], axis=1)
Net=GSN.WeightMatrix(D, maxit, tolerance)  
Weight_Matrix=Net.ConstructNet()
Weight_Matrix.to_csv('Weight_Matrix_Simulated.csv')
```
#### Sample python script to visualize the result as network
```
Signatures=SigData.columns # or Signatures=['SBS1', 'SBS2',...,] (names of signature nodes)
th=0.05
VisualizeNetwork(Weight_Matrix, Signatures, th)
```
* [Weight_Matrix_Simulated.cvs](Package/Weight_Matrix_Simulated.cvs): The

#### Results (.cvs format)
The following files proves the results for the analysis on breast and lung cancer data 
* [Weight_Matrix_BRCA.cvs](Results/Weight_Matrix_BRCA.cvs): The
* [Weight_Matrix_LUAD.cvs](Results/Weight_Matrix_LUAD.cvs):

- [SignatureNetwork_BRCA.html](Results/SignatureNetwork_BRCA.html):
- [SignatureNetwork_LUAD.html](Results/SignatureNetwork_LUAD.html):
#### Cancer Data (.cvs format)

The following files provides the gene expression and exposure of mutational signatures for cancer patients.

* [Gene_Expression_BRCA.csv](data/Gene_Expression_BRCA.csv): The normalized gene expression data (ICGC data portal) for 266 breast cancer (BRCA) patients. Gene expression profiles for 2,204 genes involved in either DNA metabolic or immune response processes of the Gene Ontology (GO) database were
selected for the analysis.
* [Signature_Exposure_BRCA.csv](data/Signature_Exposure_BRCA.csv): Exposure of 6 mutational signatures and HRD mutational status in breast cancer are used for the analysis
* [Gene_Expression_LUAD.csv](data/Gene_Expression_LUAD.csv): The normalized gene expression data (TCGA data portal) for 466 lung cancer (LUAD) patients. Gene expression profiles for 2,433 genes involved in either DNA metabolic or immune response processes of the Gene Ontology (GO) database were
selected for the analysis.
* [Signature_Exposure_LUAD.csv](data/Signature_Exposure_LUAD.csv): Exposure of 6 mutational signatures known to be operative in lung cancer are used for the analysis.

