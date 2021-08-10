import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import moment
from numpy import linalg as LA
from sklearn.linear_model import Lasso, LassoLarsIC
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
simplefilter("ignore", category=ConvergenceWarning)

#Estimation of a Sparce Partial correlation matrix for a given data
def WeightMatrix(D, maxit, tolerance):
    """Compute a weight matrix using "sklearn.linear_model.Lasso" ( Sparce Partial correlation matrix )
       The objective function to be minimized: (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1.
       For each response y, weight vector w is the incoming effect on y from columns of X and it is a row of weight matrix W 
    
    Parameters (Input)
    ----------
    D: (N_patients, N_variables) a dataframe representing gene expressions and mutational exposures.
    Its rows denote samples while the columns denote variables,
    
    maxit: The maximum number of iterations  
    
    tolerance: The tolerance for the optimization
    
    Returns (Output)
    -------
    W: (N_variables, N_variables) a dataframe representing relationship between variables. 
    its non-zero entries in i-th row represents the weights of incoming effects on i-th variables from the other variables """
       
    print('Running for estimating the SPC coefficients ...')
    [N_patients,N_variables]=D.shape
    w=np.zeros((N_variables, N_variables), dtype=float)
    Data=np.array(D)
    for k in range(N_variables):
        X=np.delete(Data,k,1)
        y=Data[:,k]
        model = LassoLarsIC(criterion='aic')
        model.fit(X, y)
        dense_lasso = Lasso(alpha=model.alpha_, fit_intercept=True, max_iter=maxit, tol=tolerance)
        dense_lasso.fit(X, y)
        w[k,:]=np.insert(dense_lasso.coef_, k, 0)
    W = pd.DataFrame(w, columns=D.columns, index=D.columns)
    return W
def Directionality(W, Data):
    """Refine the weight matrix W by removing edges with low interaction weight using a compounder-free for directionality
    For reference of the compounder-free model, see https://www.tandfonline.com/doi/abs/10.1080/00273171.2018.1528542?journalCode=hmbr20
    
    Parameters (Input)
    ----------
    W: (N_variables, N_variables) the Sparse Partial Correlation matrix obtained using WeightMatrix(D, maxit, tolerance)
    Data: (N_patients, N_variables) a dataframe representing gene expressions and mutational exposures.This data matrix is the same as D used in WeightMatrix
    Its rows denote samples while the columns denote variables,
    
    Returns (Output)
    -------
    W: (N_variables, N_variables) a revised version of the Sparse Partial Correlation matrix W"""
    eps=0.0521 #Threshold for using confounder free model (In the analysis on breast cancer date, we used eps=0.0391. However, users can use another value for threshhold)
    w=np.array(W); A=np.array(W)
    A[A!=0]=1; A=A+A.T; A[A!=2]=0
    for i in range(len(W.iloc[:,1])):
        A[i,0:i]=0
    index=np.transpose(np.nonzero(A))
    Genes=W.columns;
    for k in range(len(index)): 
        [x,y]=np.array(index[k])
        if abs(W.iloc[x,y])>eps and abs(W.iloc[y,x])>eps:
            Wx=w[x,:]; Wx[y]=0; X=Data.iloc[:,x]-Data.dot(Wx)
            Wy=w[y,:]; Wy[x]=0; Y=D.iloc[:,y]-D.dot(Wy)
            GamX=abs(moment(X, moment=3)/pow(stats.tstd(X),3))
            GamY=abs(moment(Y, moment=3)/pow(stats.tstd(Y),3))
            DelX=abs(moment(X, moment=4)/pow(stats.tstd(X),4)-3)
            DelY=abs(moment(Y, moment=4)/pow(stats.tstd(Y),4)-3)
            if abs(W.iloc[y,x])>=abs(W.iloc[x,y]) and GamX>GamY and DelX>DelY:
                W.iloc[x,y]=0
            elif abs(W.iloc[y,x])<=abs(W.iloc[x,y]) and GamX<GamY and DelX<DelY:
                W.iloc[y,x]=0
        elif abs(W.iloc[x,y])>eps and abs(W.iloc[y,x])<=eps:  
            W.iloc[y,x]=0
        elif abs(W.iloc[x,y])<=eps and abs(W.iloc[y,x])>eps: 
            W.iloc[x,y]=0
    return W

def Inside_L1ball(v): 
    """ Compute the Euclidean projection a vector (columns and rows of W) into the unit L1 norm ball
    Solves the optimization problem:
        minimize: (1/2) * || w - v ||_2^2 , s.t. || w ||_1 <= 1
    
    For reference, see https://stanford.edu/~jduchi/projects/DuchiShSiCh08.pdf
    Parameters (Input)
    ----------
    v: (n,) numpy array,
       n-dimensional vector to project
    Returns (Output)
    -------
    x: (n,) numpy array,
       Euclidean projection of v into the unit L1 norm ball"""
    n, = v.shape  
    s = np.abs(v)
    if s.sum() <= 1:
        return v
    else:
        u = np.sort(s)[::-1]
        cSum = np.cumsum(u)
        rho = np.nonzero(u * np.arange(1, n+1) > (cSum - 1))[0][-1]
        theta = float(cSum[rho] - 1) / rho
        x = (s - theta).clip(min=0)
        x *= np.sign(v)
    return x

def AM_Normalization(W):
    """Normalize the weight matrix W based on Alternate Minimization strategy
    Parameters (Input)
    ----------
    W: (N_variables, N_variables) the matrix after refinement on Sparse Partial Correlation matrix 
    
    Returns (Output)
    -------
    W: (N_variables, N_variables) a L1 norm ball normalized weight matrix consisting of directional interaction weights """
    NormW=0; k=0
    while NormW<LA.norm(W) or k<1:
        for i in range(W.shape[0]):
            W.iloc[i,:]=Inside_L1ball(np.array(W.iloc[i,:]))
        for i in range(W.shape[1]):
            W.iloc[:,i]=Inside_L1ball(np.array(W.iloc[:,i]))
        NormW=LA.norm(W)
        k=k+1
    return W

#Main body of code 

maxit=10000; tolerance=1e-12; 
D = pd.read_csv('DATA_LUAD.csv', index_col=0)   # Gene expression and mutational signature exposures as input
W=WeightMatrix(D, maxit, tolerance)               # Compute weight matrix W (a sparse partial correlation matrix) 
W=Directionality(W, D)                            # Decide directionality for bidirectional edges in W if solvable
W= AM_Normalization(W)                            # Project W onto L1 norm ball to extract more significant interactions

#Save the weight matrix W as a directed-weighted graph representing weights of directional effects between each nodes
W.to_csv('Weight_Matrix_LUAD.csv')       
