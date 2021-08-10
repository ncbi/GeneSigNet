""" 
Parameters
----------
W: (N_variables, N_variables) The weight matrix representing weights of directional influence within nodes 
    
Returns
-------
Regulators and targets of focused genes and mutational signatures with their interaction weights """
import numpy as np
import pandas as pd

W=pd.read_csv('Weight_Matrix_LUAD.csv',index_col=0) 

#-----------------------------------------------------------------------------------------------
# The focus genes or signatures which we are interested in indentifing their causal and affected gene or signatures
# Here the focus is to obtain upstream and downstream of the 6 signatures
FocusedNodes = ['SBS1','SBS2','SBS4','SBS5','SBS13','SBS40']
#For each focused node, upstream and downstream nodes are clustered as Positive causal /PC/, Negative causal /NC/, Positive affected /PA/ and Negartive affected /NA/, and
#are included in "Causals_and_Affecteds.csv" file as output
#-----------------------------------------------------------------------------------------------

GeneSets=[]
ColName=[]
df_all=pd.DataFrame(np.zeros(20), columns=['~'])
for node in FocusedNodes:
    x=W.loc[node,:].sort_values(ascending=False)
    df1=pd.DataFrame(np.column_stack((x[x>0].index,np.array(x[x>0]))), columns=[node+'/PC/',node+'/PCW/'])
    GeneSets.append(list(x[x>0].index)); ColName.append(node+'/PC/')
    x=W.loc[node,:].sort_values()
    df2=pd.DataFrame(np.column_stack((x[x<0].index,np.array(x[x<0]))), columns=[node+'/NC/',node+'/NCW/'])
    GeneSets.append(list(x[x<0].index)); ColName.append(node+'/NC/')
    Regulators=pd.concat([df1, df2], axis=1, sort=False)
    x=W.loc[:,node].sort_values(ascending=False)
    df1=pd.DataFrame(np.column_stack((x[x>0].index,np.array(x[x>0]))), columns=[node+'/PA/',node+'/PAW/'])
    GeneSets.append(list(x[x>0].index))
    ColName.append(node+'/PT/')
    x=W.loc[:,node].sort_values()
    df2=pd.DataFrame(np.column_stack((x[x<0].index,np.array(x[x<0]))), columns=[node+'/NA/',node+'/NAW/'])
    GeneSets.append(list(x[x<0].index))
    ColName.append(node+'/NT/')
    Targets=pd.concat([df1, df2], axis=1, sort=False)
    df=pd.concat([Regulators, Targets], axis=1, sort=False)
    df_all=pd.concat([df_all, df], axis=1, sort=False)
df_all=df_all.drop(columns=['~'])
df_all = df_all.dropna(axis='columns', how='all')
df_all.to_csv('Causals_and_Affecteds_LUAD.csv',index=False)
