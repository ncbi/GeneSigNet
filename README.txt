--------------------------------------
Introduction
--------------------------------------

This README file includes the description of GeneSig package,
which take gene expressions and mutational signature exposures from cancer patients as an input
and construct a weighted-directed network representing casual the relationships between genes and signatures as output  

References:

GeneSig: A network-based approach for inferring causality flows over gene-signature network
Bayarbaatar Amgalan, Yoo-Ah Kim and Teresa M. Przytycka


Questions? bayarbaatar.amgalan@nih.gov

--------------------------------------
How to run the codes
--------------------------------------

1. Construct a weighted-directed network representing casual the relationship genes and signatures
Input Data:   "DATA_LUAD.csv" or "DATA_BRCA": (N_patients, N_variables) a dataframe representing gene expressions and 
mutational exposures. Its rows denote samples while the columns denote variables.

Output:  "Weight_Matrix_LUAD.csv" or "Weight_Matrix_BRCA.csv": (N_variables, N_variables) a dataframe representing relationship between 
variables and its non-zero entries in i-th row represents the weights of incoming effects on i-th variables 
from the other variables. Thus, The nonzero value in i-th row and j-th column represents the weight of effect on i-th variable from j-th variable (j-->i)

Go to the directory where Construct_Network.py and Input Data are located. Then run Construct_Network.py as following

>> python Construct_Network.py


2. Obtain the list of causals (regulators) and affecteds (regulators) of focused nodes (The nodes can be genes and signatures) 

Input:   "Weight_Matrix_LUAD.csv" or "Weight_Matrix_BRCA.csv": (N_variables, N_variables) a dataframe representing relationship between 
variables and its non-zero entries in i-th row represents the weights of incoming effects on i-th variables 
from the other variables. Thus, w_ij represents the weight of effect on i-th variable from j-th variable (j-->i)

Input:   "FocusedNodes" must be provided in python code (Find_Causals_and_Affecteds.py), not as a separate file. A set of focus genes or signatures that we are interested in identifying the regulators and targets of them. For our settings, the mutation signatures (SigStates) are considered as the focused nodes

Output:  "Causals_and_Affecteds_LUAD.csv" or "Causals_and_Affecteds_BRCA.csv" depending on which cancer is chosen for the analyis. Regulators and targets of focused genes and mutational signatures with their interaction weights

>> python Find_Causals_and_Affecteds.py

--------------------------------------
Source Files
--------------------------------------

* Construct_Network.py: construct a weight matrix representing a weighted-directed network
* Find_Causals_and_Affecteds.py: obtain the causal and affected nodes for the focused gene and signatures

--------------------------------------
Data/Results
--------------------------------------
data/
      "DATA_LUAD.csv": a data matrix representing gene expressions and mutational exposures. 
      Its rows denote samples while the columns denote variables.

      "FocusedNodes" must be provided (specified) in the python code (Find_Causals_and_Affecteds.py) A set of focus genes or
      signatures that we are interested in identifying the causal (regulators) and affected (targets) nodes of them

results/
	"Weight_Matrix_LUAD.csv" or "Weight_Matrix_BRCA.csv": a matrix representing relationship between 
      variables and its non-zero entries of i-th row represent the weights of incoming effects on i-th variables 
      from the other variables (nodes in the Gene Signature Network)

      "Causals_and_Affecteds_LUAD.csv" or "Causals_and_Affecteds_LUAD.csv": Causal and affected set for the focused genes and mutational signatures with their
      interaction weights.
