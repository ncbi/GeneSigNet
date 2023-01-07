# GeneSigNet
A computational method for constructing Gene - Signature causal network from gene expression and exposure of mutational signatures. 
#### Source install

```
pip install GeneSigNet
```

#### Vignettes

```
To use GeneSigNet:
```
import pandas as pd
import GeneSigNet as GSN

maxit=10000; tolerance=1e-12; 
D=pd.concat([ExpData, SigData], axis=1)
Net=GSN.WeightMatrix(D, maxit, tolerance)  
Weight_Matrix=Net.ConstructNet()
Weight_Matrix.to_csv('Weight_Matrix_Simulated.csv')

#### Visualize the result as network

Signatures=SigData.columns
th=0.05
VisualizeNetwork(Weight_Matrix, Signatures, th)
