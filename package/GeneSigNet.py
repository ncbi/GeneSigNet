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
class WeightMatrix:
    def __init__(self, Data, maxit, tolerance):

        self.D = Data
        self.maxit =maxit
        self.tol=tolerance
        self.var=Data.columns
        self.W=0

    def ConstructNet(self):
        """Compute a weight matrix using "sklearn.linear_model.Lasso" ( Sparce Partial Correlation Matrix )
    
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

        N_variables=len(self.var)  # total number of nodes in the network (the number of variables in the multivatiate analysis) 
        w=np.zeros((N_variables, N_variables), dtype=float)
        Data=np.array(self.D)
        for k in range(N_variables):
            X=np.delete(Data,k,1)
            y=Data[:,k]
            model = LassoLarsIC(criterion='aic')
            model.fit(X, y)
            dense_lasso = Lasso(alpha=model.alpha_, fit_intercept=True, max_iter=self.maxit, tol=self.tol)
            dense_lasso.fit(X, y)
            w[k,:]=np.insert(dense_lasso.coef_, k, 0)
        self.W = pd.DataFrame(w, columns=self.var, index=self.var)
        self.W=Directionality(self.W, self.D)
        self.W=AM_Normalization(self.W)
        return self.W
def Directionality(W, D):
    """Resolve directions for bidirected edges in the weight matrix W by using the proposed partial higher-moment strategy
    
    Parameters (Input)
    ----------
    W: (N_variables, N_variables) the Sparse Partial Correlation matrix obtained using WeightMatrix(D, maxit, tolerance)
    Data: (N_patients, N_variables) a dataframe representing gene expressions and mutational exposures.This data matrix is the same as D used in WeightMatrix
    Its rows denote samples while the columns denote variables,
    
    Returns (Output)
    -------
    M: (N_variables, N_variables) a revised version of the Sparse Partial Correlation Matrix W"""

    tau=0.001 # A threshold parameter for edge weights (described in Section S1.3 of the paper: We used tau=0.0391 for BRCA and 0.0521 for LUAD). tau=0.01 is default 
    w=np.array(W); A=np.array(W)
    A[A!=0]=1; A=A+A.T; A[A!=2]=0
    for i in range(len(W.iloc[:,1])):
        A[i,0:i]=0
    index=np.transpose(np.nonzero(A))
    var=W.columns #list of variables (Genes and Signatures)
    M=pd.DataFrame(w, columns=var, index=var)
    for k in range(len(index)): 
        [x,y]=np.array(index[k])
        """For each bi-directed edge in W, use the partial higher moments to decide direction if the weights in both direction are strong"""
        if abs(W.iloc[x,y])>tau and abs(W.iloc[y,x])>tau: # For each bi-directed edge in W, use the partial higher moments to decide direction if the weights in both direction are strong
            Wx=w[x,:]; Wx[y]=0; rX=D.iloc[:,x]-D.dot(Wx)   # Calculate residual rX by removing confounding effects from the observation data of x-th variable
            Wy=w[y,:]; Wy[x]=0; rY=D.iloc[:,y]-D.dot(Wy)   # Calculate residual rY by removing confounding effects from the observation data of y-th variable
            GamX=abs(moment(rX, moment=3)/pow(stats.tstd(rX),3))    # Skewnesses of rX
            GamY=abs(moment(rY, moment=3)/pow(stats.tstd(rY),3))    # Skewnesses of rY
            DelX=abs(moment(rX, moment=4)/pow(stats.tstd(rX),4)-3)  # Kurtosises of rX
            DelY=abs(moment(rY, moment=4)/pow(stats.tstd(rY),4)-3)  # Kurtosises of rX
            if abs(W.iloc[y,x])>=abs(W.iloc[x,y]) and GamX>GamY and DelX>DelY:    # Remove edge from y to x if the 3 conditions satisfied
                M.iloc[x,y]=0    
            elif abs(W.iloc[y,x])<=abs(W.iloc[x,y]) and GamX<GamY and DelX<DelY:  # Remove edge from x to y if the 3 conditions satisfied
                M.iloc[y,x]=0                                                     
        elif abs(W.iloc[x,y])>tau and abs(W.iloc[y,x])<=tau:                      # Remove edge from x to y if the 2 conditions satisfied
            M.iloc[y,x]=0
        elif abs(W.iloc[x,y])<=tau and abs(W.iloc[y,x])>tau:                      # Remove edge from y to x if the 2 conditions satisfied
            M.iloc[x,y]=0
    return M

def Inside_L1ball(v): 
    """ Compute the Euclidean projection a vector (columns and rows of W) into the unit L1 norm ball
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
        theta = float(cSum[rho] - 1) / (rho+1)
        x = (s - theta).clip(min=0)
        x *= np.sign(v)
    return x

def AM_Normalization(W):
    """Normalize the weight matrix W based on Alternate Minimization strategy
    Parameters (Input)
    ----------
    W: (N_variables, N_variables) the matrix after refinement on Sparse Partial Correlation Matrix 
    
    Returns (Output)
    -------
    W: (N_variables, N_variables) sparse weight matrix whose columns and rows are normalized onto L1 norm ball """
    NormW=0; k=0
    while NormW<LA.norm(W) or k<1:
        for i in range(W.shape[0]): # Normalize each row onto L1 norm ball 
            W.iloc[i,:]=Inside_L1ball(np.array(W.iloc[i,:]))
        for i in range(W.shape[1]): # Normalize each column onto L1 norm ball 
            W.iloc[:,i]=Inside_L1ball(np.array(W.iloc[:,i]))
        NormW=LA.norm(W)
        k=k+1
    return W