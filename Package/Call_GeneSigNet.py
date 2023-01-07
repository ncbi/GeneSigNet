import pandas as pd
import GeneSigNet as GSN
from pyvis.network import Network

def VisualizeNetwork(W, Signatures, th):
    net=Network(height='1000px', width='100%', bgcolor='white', font_color='black', directed = True)
    net.barnes_hut()

    sources=[]; targets=[]; weights=[]
    for sig in Signatures:
        t=W.loc[:,sig]
        t=t[abs(t)>th]
        r=W.loc[sig,:]
        r=r[abs(r)>th]
        sources =sources+[sig] * len(t)+list(r.index)
        targets=targets+list(t.index)+[sig] * (len(r))
        weights=weights+list(t)+list(r)
    edges=zip(sources,targets, weights)

    net.toggle_physics(False)
    for e in edges:
        reg = e[0]
        tar=e[1]
        w=e[2]
        if reg in Signatures:
            net.add_node(reg, reg, title=reg, color='orange', shape = 'owl')
        else:
            net.add_node(reg, reg, title=reg, color='green')
        if tar in Signatures:
            net.add_node(tar, tar, title=tar, color='orange', shape = 'owl')
        else:
            net.add_node(tar, tar, title=reg, color='green')    
        if w>0:
            net.add_edge(reg, tar, value=w, color='red')
        else:
            net.add_edge(reg, tar, value=w, color='blue')
    
    neighbor_map = net.get_adj_list()
    for node in net.nodes:
        node['title'] += ' Neighbor:<>'+'<br>'.join(neighbor_map[node['id']])
        node['value'] = len(neighbor_map[node['id']])
    net.show('SignatureNetwork_Simulated.html')

    #from IPython.core.display import display, HTML
    #display(HTML('SignatureNetwork.html'))
    #net.toggle_physics(False)

maxit=10000; tolerance=1e-12; 
ExpData = pd.read_csv('Gene_Expression_Simulated_Data.csv', index_col=0)  
SigData = pd.read_csv('Signature_Exposure_Simulated_Data.csv', index_col=0)
D=pd.concat([ExpData, SigData], axis=1)
Net=GSN.WeightMatrix(D, maxit, tolerance)  
Weight_Matrix=Net.ConstructNet()
Weight_Matrix.to_csv('Weight_Matrix_Simulated.csv')
Signatures=SigData.columns
th=0.05
VisualizeNetwork(Weight_Matrix, Signatures, th)