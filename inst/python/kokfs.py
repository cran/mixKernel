#############################################################################################################
# Author :
#   Céline Brouard, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#   Jerome Mariette, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#   Nathalie Villa-Vialaneix, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#
# Copyright (C) 2017
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


import scipy as sp
import autograd
import autograd.numpy as np
from autograd import grad
from autograd.numpy import linalg as la
from autograd.numpy import ndarray
from sklearn.model_selection import KFold

#### Kernel functions
def linear(x1, x2):
    """Computation of the kernel matrix for the linear kernel"""
    return np.dot(x1.T,x2)

def linear_sel(wl,wr,x1,x2):
    """Computation of the linear kernel where the variables are weighted"""
    x1_weighted = np.diag(wl)@x1
    x2_weighted = np.diag(wr)@x2
    return linear(x1_weighted,x2_weighted)

def gaussian_radial_basis(x1,x2,gamma=1):
    """Computation of the kernel matrix for the Gaussian kernel"""
    return np.exp(-gamma*euclidean_distkokfs(x1,x2))

def gaussian_radial_basis_sel(wl,wr,x1,x2,gamma=1):
    """Computation of the Gaussian kernel where the variables are weighted"""
    x1_weighted = np.diag(wl)@x1
    x2_weighted = np.diag(wr)@x2
    return gaussian_radial_basis(x1_weighted,x2_weighted,gamma)

def kernel(kx_func):
    if kx_func == "linear":
        kernel_k = linear
    elif kx_func == "gaussian.radial.basis":
        kernel_k = gaussian_radial_basis
    return kernel_k

def kernel_sel(kx_func):
    if kx_func == "linear":
        kernel_k = linear_sel
    elif kx_func == "gaussian.radial.basis":
        kernel_k = gaussian_radial_basis_sel
    return kernel_k
####

#### optim functions
def compute_f(wl,wr,x,kernel_sel,KY,lambda1,**kwargs):
    """
    Compute the loss f(h,w) + lambda1 ||h||^2 where h is replaced by the solution obtained when w is fixed
    This function is called by the select_variables function
    """

    KX_l_r = kernel_sel(wl,wr,x,x,**kwargs)
    KX_l_l = kernel_sel(wl,wl,x,x,**kwargs)
    n = KX_l_l.shape[0]
    I_n = np.identity(n)
    A = la.solve(lambda1*I_n + KX_l_l, KX_l_r)
    A1 = la.solve(lambda1*I_n + KX_l_l, KX_l_l)
    A2 = la.solve(lambda1*I_n + KX_l_l, KY)
    return np.trace(np.transpose(I_n - A) @ KY @ (I_n - A)) + lambda1*np.trace(A1@A2)

def learn_w(x, kx_func, kx_param, KY, w0, lambda1, lambda2, backtrack=True, nbitermax=1000,
                     stopvarw=1e-9, stopvarj=1e-9, t0=10., m_back=1,sigma=1e-9, eta=2, nbitermax_back=100, bbrule=False, verbose=False):
    """
    Main function for the proposed approach that outputs the learnt weight vector w
    
    INPUTS:
    ----------
    x: input data
    kx_func: input kernel ('linear' or 'gaussian.radial.basis')
    kx_param: dict containing the potential parameter(s) of the input kernel
    KY: output kernel matrix
    w0: initial value for the weight vector w
    lambda1: value of the lambda1 regularization parameter (>0)
    lambda2: value of the lambda2 regularization parameter (>0)
    other options: see the description given in the solver function
    
    OUTPUTS:
    ----------
    w: optimal solution w
    loss: value of the cost function at each iteration
        """
    
    # Input kernel
    kx_kernel = kernel(kx_func)
    kx_kernel_sel = kernel_sel(kx_func)
    
    grad_f_w = grad(compute_f,1)
    
    f = lambda wl,wr:compute_f(wl,wr,x,kx_kernel_sel,KY,lambda1,**kx_param)
    df = lambda wl,wr:grad_f_w(wl,wr,x,kx_kernel_sel,KY,lambda1,**kx_param)
    g = lambda x, lambd : lambd * np.sum(np.abs(x))
    prox_g = prox_l1_pos
    
    w,loss = solver(f,df,g,prox_g,w0,lambda2,backtrack,nbitermax,
                    stopvarw,stopvarj,t0,m_back,sigma,eta,nbitermax_back,bbrule,verbose)
    return w,loss


def solver(f,df,g,prox_g,w0,lambda2,backtrack=True, nbitermax=1000,
           stopvarw=1e-9, stopvarj=1e-9, t0=10., m_back=1,
           sigma=1e-9, eta=2, nbitermax_back=100, bbrule=False,verbose=False,**kwargs):
    """
    Solver for the alternating proposed algorithm
    This solver has been adapted from the solver for Forward Backward Splitting implemented in the PyOptim library
    
    INPUTS:
    ----------
    f: function computed using the function compute_f
    df: gradient of f
    g: nonsmooth function (for example l1 norm)
    prox_g: proximal of g
    w0: initial weight vector
    lambda2: regularization parameter (>0)
    backtrack: boolean indicating whether backtracking should be performed (default: True)
    nbitermax: max number of iterations
    stopvarw: stopping criterion for relative variation of the norm of w
    stopvarj: stopping criterion for relative variation of the cost
    t0: initial step size (optional)
    m_back: window size for backtrack
    sigma: descent parameter for backtrack
    eta: value multiplying t during backtrack
    nbitermax_back: max number of backtrack iterations
    bbrule: optional boolean indicating whether bbrule should be used for updating the step size
    verbose: print optimization information (optional)
    
    OUTPUTS:
    ----------
    w: optimal solution w
    loss: value of the cost function at each iteration
        """
    
    # Initialization
    w = w0.copy()

    grad = df(w0,w0,**kwargs)

    loss = list()
    loss.append(f(w0,w0,**kwargs) + g(w0,lambda2))

    t = t0

    # Loop
    loop = True
    it = 1
    while loop:
        w_1 = w.copy()
        grad_1 = grad.copy()
    
        # gradient of the function f
        grad = df(w_1,w_1,**kwargs)
        
        # proximal operator
        w = prox_g(w_1 - grad / t, lambda2 /t)
        
        # cost computation
        loss.append(f(w_1,w,**kwargs) + g(w,lambda2,**kwargs))
    
        # line search backtrack for determining the step size t
        it2 = 0
        thr_back = np.max([loss[-2 - k] - sigma / 2 * t * norm(w - w_1)
                           ** 2 for k in range(min(m_back, len(loss) - 1))])
        while loss[-1] > thr_back and it2 < nbitermax_back and backtrack:
            t = t * eta
            w = prox_g(w_1 - grad / t, lambda2 / t)
            loss[-1] = f(w_1,w, **kwargs) + g(w, lambda2, **kwargs)
            thr_back = np.max([loss[-2 - k] - sigma / 2 * t * norm(w - w_1)
                              ** 2 for k in range(min(m_back, len(loss) - 1))
                              ])
            it2 += 1
        if it2 == nbitermax_back:
            print("Warning: backtrack failed")
    
        # BB rule
        xbb = w - w_1
        ybb = grad - grad_1
        if it >= 1 and norm(xbb) > 1e-12 and norm(ybb) > 1e-12 and bbrule:
            t = abs(np.sum(xbb * ybb) / np.sum(xbb * xbb))

        if abs(loss[-1] - loss[-2]) / abs(loss[-2]) < stopvarj:
            loop = False
            if verbose:
                print("delta loss convergence")
        
        if it >= nbitermax:
            loop = False
            if verbose:
                print("Max number of iteration reached")

        # increment iteration
        it += 1
    
    return w,loss

def select_lambda1(x, kx_func, kx_param, KY, val_lambda1, num_folds=5):
    """
    Selection of the regularization parameter lambda1 that controls the complexity
    of the h function by minimizing the MSE in a cross-validation experiment
    (no variable selection is considered when selecting this hyperparameter)
    
    INPUTS:
    ----------
    x: input data
    kx_func: input kernel ('linear' or 'gaussian.radial.basis')
    kx_param: dict containing the potential parameter(s) of the input kernel
    KY: output kernel matrix
    val_lambda1: grid of positive values among which lambda1 will be selected
    num_folds: number of folds for the cross-validation experiment (default is 5)
    
    OUTPUTS:
    ----------
    lambda1_opt: selected value for the lambda1 regularization parameter
    """
    
    kx_kernel = kernel(kx_func)
    KX = kx_kernel(x,x,**kx_param)
    
    # Cross-validation
    cv = KFold(n_splits=num_folds, shuffle=True)
    mse_test_cv = np.zeros((len(val_lambda1),num_folds))

    for ind_cv, (train_index, test_index) in enumerate(cv.split(x.T)):

        n_train = len(train_index)
        n_test = len(test_index)
        
        KX_train = KX[np.ix_(train_index,train_index)]
        KX_train_test = KX[np.ix_(train_index,test_index)]
        
        KY_train = KY[np.ix_(train_index,train_index)]
        KY_test = KY[np.ix_(test_index,test_index)]
        KY_test_train = KY[np.ix_(test_index,train_index)]
        
        for ind_l1, lambda1 in enumerate(val_lambda1):
            
            # Training step
            I_train = np.identity(n_train)
            A = la.solve(lambda1*I_train + KX_train, KX_train_test)
            
            # Computation of the test MSE for different values of lambda1 and folds
            mse_test_cv[ind_l1,ind_cv] = 1/n_test*np.trace(KY_test -2*KY_test_train @ A + np.transpose(A) @ KY_train @ A)

    # Selection of the value that minimizes the averaged MSE
    mean_mse_cv = np.mean(mse_test_cv,1)
    ind_lambda1_opt = np.argmin(mean_mse_cv)
    lambda1_opt = val_lambda1[ind_lambda1_opt]

    return lambda1_opt
####

#### utility functions
def center_scale(x):
    """center-scaling"""
    num_gene,num_sample = x.shape
    x_center = x - (np.tile(np.mean(x,axis=0),[num_gene,1]))
    x_scale = x_center / (np.tile(np.std(x_center,axis=0),[num_gene,1]))
    return x_scale

def normmat(K):
    """matrix normalization to have diagonal equal to 1"""
    D = np.diag(K)[np.newaxis].T
    Kn = K / np.sqrt(np.dot(D, D.T))
    return Kn

def euclidean_distkokfs(x1, x2):
    """Computes all pairwise euclidean distances between the columns of matrices x1 and x2"""
    [nrow1,ncol1] = x1.shape
    [nrow2,ncol2] = x2.shape
    x1p2 = np.sum(np.square(x1), 0)
    x2p2 = np.sum(np.square(x2), 0)
    return np.tile(x1p2.reshape((-1, 1)),(1,ncol2)) + np.tile(x2p2.reshape((1, -1)),(ncol1,1)) - 2 * np.dot(x1.T, x2)

def norm(x):
    """l2 norm of vector (Frobenius for matrices)"""
    return np.sqrt(np.sum(np.square(x)))

def prox_l1_pos(v, lambd):
    """proximal operator of l1 norm (lasso) + positivity constraint"""
    return np.maximum(v - lambd, 0)

def select_gamma_rbf(x):
    """define the parameter for a Gaussian kernel"""
    n = x.shape[1]
    Dist = euclidean_distkokfs(x, x)
    gamma = (n*(n-1))/(2*np.sum(np.triu(Dist,1)))
    return gamma
####

#### main function
which = lambda lst:list(np.where(lst)[0])
def kokfspy(x, y, kx_func='gaussian.radial.basis', ky_func='gaussian.radial.basis', keepX=40, nstep=50, **kwargs):
    """
    Main function for feature selection
    
    INPUTS:
    ----------
    x: input data
    kx_func: input kernel ('linear' or 'gaussian.radial.basis')
    y: output data
    ky_func: output kernel ('linear' or 'gaussian.radial.basis')
    keepX: number of features to be selected
    
    OUTPUTS:
    ----------
    selected_indices: indices of the selected features
    """
    # force x and y to be a numpy array
    x = np.array(x)
    y = np.array(y)
    p = x.shape[0]
    nstep = int(nstep)
    keepX = int(keepX)
    
    # Setting the kernel parameters
    if kx_func == 'gaussian.radial.basis':
        kx_param = {'gamma':select_gamma_rbf(x)}
    else:
        kx_param = {}
    
    # Computing the Gram matrix of the output kernel
    ky_kernel = kernel(ky_func)
    if ky_func == 'gaussian.radial.basis':
        KY = ky_kernel(y,y,select_gamma_rbf(y))
    else:
        KY = ky_kernel(y,y)
    
    # Selecting the lambda1 regularization parameter
    val_lambda1 = np.logspace(-3, 4, num=25, endpoint=True, base=10.0)
    num_folds = 5
    lambda1_opt = select_lambda1(x, kx_func, kx_param, KY, val_lambda1, num_folds)
    
    # Learning the weights for a grid of lambda2 parameters
    val_lambda2 =  np.logspace(-5, 2, num=nstep, endpoint=True, base=10.0)
    weights = np.zeros((len(val_lambda2),p))
    w0 = np.ones(p)
    for ind_lambda2 in range(len(val_lambda2)):
        w,loss = learn_w(x, kx_func, kx_param, KY, w0, lambda1_opt, val_lambda2[ind_lambda2], **kwargs)
        weights[ind_lambda2,:] = w
        w0 = np.copy(w) # warm restart
        if (len(which(w>0)) <= keepX):
            break
    
    # Extracting the indices of the selected features
    stop = np.zeros(p)
    for i in range(p):
        j = len(val_lambda2) - 1
        while weights[j,i] == 0 :
            j -=1
            if j == 0:
                break
        stop[i] = j
    
    selected_indices = (-stop).argsort()[:int(keepX)]

    return selected_indices
####
