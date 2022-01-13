#############################################################################################################
# Author :
#   Jerome Mariette, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#   Nathalie Vialaneix, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
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


import autograd
import autograd.numpy as np
from scipy.spatial.distance import pdist
from scipy.optimize import minimize
from scipy.linalg import eigh

#### utility functions
def euclidean_dist(x1, x2):
    x1p2 = np.sum(np.square(x1), 1)
    x2p2 = np.sum(np.square(x2), 1)
    return x1p2.reshape((-1, 1)) + x2p2.reshape((1, -1)) - 2 * np.dot(x1, x2.T)
####

#### linear kernel functions
def linear(x1, x2):
    return np.dot(x1,np.transpose(x2))

def linear_sel(w,x1,x2):
    d=w.shape[0]
    return linear(np.reshape(w,(1,d))*x1,np.reshape(w,(1,d))*x2)
####

#### gaussian.radial.basis kernel functions
def gaussian_radial_basis(x1, x2, sigma=1):
    return np.exp(-sigma*euclidean_dist(x1,x2))

def gaussian_radial_basis_sel(w, x1, x2, sigma=1):
    d=w.shape[0]
    return np.exp(-sigma*euclidean_dist(np.reshape(w,(1,d))*x1,np.reshape(w,(1,d))*x2))
####

#### bray-curtis dissimilarity functions
def bray(x1,x2):
    eps=1
    return np.sum(np.abs(x1[:,None,:]-x2[None,:,:]),2)/(np.sum((x1[:,None,:]+x2[None,:,:]),2)+eps) 

def bray_sel(w,x1,x2):
    d=w.shape[0]
    return bray(np.reshape(w,(1,d))*x1,np.reshape(w,(1,d))*x2)
####

#### bray-curtis kernel functions (for kpca)
def brayK(x1,x2):
    eps=0
    diss=np.sum(np.abs(x1[:,None,:]-x2[None,:,:]),2)/(np.sum((x1[:,None,:]+x2[None,:,:]),2)+eps)
    return -0.5 * np.dot(np.dot(np.identity(diss.shape[0]) - 1 / diss.shape[0], diss), np.identity(diss.shape[0]) - 1 / diss.shape[0]) 
    
def brayK_sel(w,x1,x2):
    d=w.shape[0]
    return brayK(np.reshape(w,(1,d))*x1,np.reshape(w,(1,d))*x2)
####

#### optim functions
def reg_l1(x, lambd=1, **kwargs):
    """
    l1 regularization (LASSO)
    x : vector to project
    lambd: regularization term
    """
    return lambd * np.sum(np.abs(x))


def prox_l1(x, lambd=1, **kwargs):
    """
    l1 regularization proximal operator (LASSO)
    x : vector to project
    lambd: regularization term
    """
    return np.sign(x) * np.maximum(np.abs(x) - lambd, 0)

def norm(x):
    """l2 norm of vector (Frobenius for matrices)"""
    return np.sqrt(np.sum(np.square(x)))

# format string for printing in verbose mode
prnt_str_name = "|{it:>5}|{loss:>13}|{dloss:>13}|{step:>13}|\n" \
                "|-----|-------------|-------------|-------------|"
prnt_str_loop = "|{it:5d}|{loss: 10e}|{dloss: 10e}|{step: 10e}|"

def fmin_prox(f, df, g, prox_g, x0, lambd=1., backtrack=True, nbitermax=1000,
              stopvarx=1e-9, stopvarj=1e-9, t0=10., verbose=False, m_back=1,
              sigma=1e-9, eta=2, nbitermax_back=100, bbrule=True, log=False,
              **kwargs):
    r""" Minimize a sum of smooth and nonsmooth function using proximal splitting
    The algorithm is the classical Forward Backward Splitting [1]_
    with BB rule for step estimation [2]_.
    Solve the optimization problem:
    .. math::
        min_x \quad  f(x)+\lambda g(x)
    where:
    - f is differentiable (df) and Lipshictz gradient
    - g is non-differentiable but has a proximal operator (prox_g)
    prox_g is a function providing the solution to problem
    .. math::
        min_x \quad  \frac{1}{2}\|x_0-x\|^2+\lambda*g(x)
    Several proximal operators are available at optim.prox
    Parameters
    ----------
    f : function
        Smooth function f: R^d -> R
    df : function
        Gradient of f, df:R^d -> R^d
    g : function
        Nonsmooth function g: R^d -> R
    prox_g : function
        Proximal of g, df:R^d -> R^d
    x_0 : (d,) numpy.array
        Initial point
    lambda : float
        Regularization parameter >0
    backtrack : boolean, optional
        Perform backtracking if true (default: True).
    bbrule : boolean, optional
        update step with bb rule.
    nbitermax : int, optional
        Max number of iterations in algorithm
    stopvarx : float, optional
        Stopping criterion for relative variation of the
        norm of x
    stopvarj : float, optional
        Stopping criterion for relative variation of the cost
    t0 : float, optional
        initial descent step
    verbose : boolean, optional
        prinrt optimization information
    m_back : int, optional
        Window size for backtrack (if <1 then non decreasing)
    sigma : float, optional
        descent parameter for backtrack
    eta : float, optional
        value multiplying t during backtrack
    nbitermax_back : int, optional
        Max number of backtrack iterations
    Returns
    -------
    x: (d,) ndarray
        Optimal solution x
    val: float
        optimal value of the objective (None if optimization error)
    log: dict
        Optional log output
    References
    ----------
    .. [1] Combettes, P. L., & Wajs, V. R. (2005). Signal recovery by proximal
        forward-backward splitting. Multiscale Modeling & Simulation, 4(4),
        1168-1200.
    .. [2] Barzilai, J., & Borwein, J. M. (1988). Two-point step size
        gradient methods. IMA journal of numerical analysis, 8(1), 141-148.
    See Also
    --------
    optim.prox : Module containing proximal operators
    optim.fmin_proj : Projected gradient (special case of proximal gradient)
    """
    x = x0.copy()
    grad = df(x, **kwargs)
    # grad[:]=0
    loss = list()
    loss.append(f(x, **kwargs) + g(x, lambd, **kwargs))
    if log:
        log = dict()
        log['loss'] = loss
    t = t0
    if verbose:
        print((prnt_str_name.format(it="It. ", loss='Loss ',
                                    dloss="Delta Loss ", step="Step ")))
        print((prnt_str_loop.format(it=0, loss=loss[-1], dloss=0, step=1 / t)))
    loop = True
    it = 1
    while loop:
        x_1 = x.copy()
        grad_1 = grad.copy()
        # gradient
        grad = df(x, **kwargs)
        # prox operator
        x = prox_g(x_1 - grad / t, lambd / t, **kwargs)
        # cost computation
        loss.append(f(x, **kwargs) + g(x, lambd, **kwargs))
        # line search backtrack
        it2 = 0
        thr_back = np.max([loss[-2 - k] - sigma / 2 * t * norm(x - x_1)
                           ** 2 for k in range(min(m_back, len(loss) - 1))])
        while loss[-1] > thr_back and it2 < nbitermax_back and backtrack:
            t = t * eta
            x = prox_g(x_1 - grad / t, lambd / t, **kwargs)
            loss[-1] = f(x, **kwargs) + g(x, lambd, **kwargs)
            thr_back = np.max([loss[-2 - k] - sigma / 2 * t * norm(x - x_1)
                               ** 2 for k in range(min(m_back, len(loss) - 1))
                               ])
            # print '\t',loss[-1],thr_back
            it2 += 1
        if it2 == nbitermax_back:
            print("Warning: backtrack failed")
        # print loss[-1],t
        # print information
        if verbose:
            if not (it) % 20:
                print((prnt_str_name.format(it="It. ", loss='Loss ',
                                            dloss="Delta Loss ",
                                            step="Step ")))
            print(prnt_str_loop.format(it=it,
                                       loss=loss[-1],
                                       dloss=(loss[-1] - loss[-2]
                                              ) / abs(loss[-2]),
                                       step=1 / t))
        # BB rule
        xbb = x - x_1
        ybb = grad - grad_1
        if it >= 1 and norm(xbb) > 1e-12 and norm(ybb) > 1e-12 and bbrule:
            t = abs(np.sum(xbb * ybb) / np.sum(xbb * xbb))
#        else:
#            t=t0
        # test convergence
        if norm(x - x_1) / norm(x) < stopvarx:
            loop = False
            if verbose:
                print("delta x convergence")
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
    if log:
        log['loss'] = loss
        return x, loss[-1], log
    else:
        return x, loss[-1]
####

#### utility functions
def kernel_dist(K):
    return np.asarray([[K[i,i]+K[j,j]-2*K[i,j] for j in range(0, K.shape[0])] for i in range(0, K.shape[0])])

def kernel_loss(w,x,K,kernel_sel,**kwargs):
    return np.sum(np.square(K-kernel_sel(w,x,x,**kwargs)))
    
def kernel_loss_log(w,x,K,kernel_sel,**kwargs):
    return np.log(np.sum(np.square(K-kernel_sel(w,x,x,**kwargs))) + 1)

def kpca_loss(w,x,X_kpca,kernel_sel,**kwargs):
    Kpca_dist = euclidean_dist(X_kpca,X_kpca)
    Kpca_dist[Kpca_dist < 0] = 0
    return np.sum( np.square( kernel_dist(kernel_sel(w,x,x,**kwargs)) - Kpca_dist ) )

def graph_loss(w,x,K,Lg,mu,kernel_sel,**kwargs):
    return np.sum(np.square(K-kernel_sel(w,x,x,**kwargs))) + mu * np.dot(np.dot(w.T,Lg),w)

def kpca_func(K, n_components=None):
    
    # first center the kernel
    n_samples = K.shape[0]
    K_fit_rows = np.sum(K, axis=0) / n_samples
    K = K - K_fit_rows - (np.sum(K, axis=1) / K_fit_rows.shape[0])[:, np.newaxis] + K_fit_rows.sum() / n_samples
    
    if n_components is None:
        n_components = K.shape[1]
    else:
        n_components = min(K.shape[1], int(n_components))
    lambdas, alphas = eigh(K, eigvals=(K.shape[1] - n_components, K.shape[1] - 1))

    # sort eigenvectors in descending order
    indices = lambdas.argsort()[::-1]
    lambdas = lambdas[indices]
    alphas = alphas[:, indices]
    
    # remove eigenvectors with a zero eigenvalue
    tolerance = 10**(-10)
    lambdas = np.round_(lambdas, int(-np.log10(tolerance)))
    
    alphas = alphas / np.sqrt(lambdas)
    X_transformed = alphas * lambdas
    return (lambdas, alphas, X_transformed)
####

#### main function
which = lambda lst:list(np.where(lst)[0])
def ukfspy(x, kernel_func, method, keepX, reg, n_components, Lg, mu, max_iter, nstep, kwargs):
    
    # force x to be a numpy array
    x = np.array(x)
    keepX = int(keepX)
    nstep = int(nstep)
    max_iter = int(max_iter)
    
    # optim parameters
    params=dict()
    params['nbitermax']=max_iter
    params['stopvarx']=1e-9
    params['stopvarj']=1e-9
    params['t0']=10.
    params['m_back']=1
    params['verbose']=False
    params['bbrule']=True
    params['log']=False
    
    if kernel_func == "linear":
        kernel = linear
        kernel_sel = linear_sel
        kwargs = {}
    elif kernel_func == "gaussian.radial.basis":
        kernel = gaussian_radial_basis
        kernel_sel = gaussian_radial_basis_sel
        # if no value provided for sigma parameter, find the sigma
        # that maximizes projections distances in the KPCA space
        if "sigma" not in kwargs:
            def opt_rbf_sigma (sigma):
                K = gaussian_radial_basis(x,x,sigma[0])
                (lambdas, alphas, X_kpca) = kpca_func(K, 2)
                all_dist = pdist(X_kpca[:,0:2], "euclidean")
                #all_dist = pdist(X_kpca[:,0:2], "euclidean")
                #return -np.sum(all_dist)
                return -np.sum(lambdas)
            kwargs = {"sigma" : minimize(opt_rbf_sigma, [0.1], method="L-BFGS-B", bounds=[(1e-12, 1)]).x[0]}
    elif kernel_func == "bray":
        if method == "kpca":
            # if kpca use the kernel version
            kernel = brayK
            kernel_sel = brayK_sel
        else:
            # if kernel and graph use the dissimilarity
            kernel = bray
            kernel_sel = bray_sel

        # only if differentiable
        params['bbrule']=False
        kwargs = {}
    
    K = kernel(x,x,**kwargs)
    
    if method == "kernel":
        grad_loss_K = autograd.grad(kernel_loss,0)
        # parameters for optim (loss/reg)
        f=lambda w:kernel_loss(w,x,K,kernel_sel,**kwargs) #  loss
        df=lambda w:grad_loss_K(w,x,K,kernel_sel,**kwargs) # grad  loss
    elif method == "kpca":
        (lambdas, alphas, X_kpca) = kpca_func(K, n_components)
        grad_loss_K = autograd.grad(kpca_loss,0)
        # parameters for optim (loss/reg)
        f=lambda w:kpca_loss(w,x,X_kpca,kernel_sel,**kwargs) #  loss
        df=lambda w:grad_loss_K(w,x,X_kpca,kernel_sel,**kwargs) # grad  loss
    elif method == "graph":
        grad_loss_K = autograd.grad(graph_loss,0)
        Lg = np.array(Lg)
        # parameters for optim (loss/reg)
        f=lambda w:graph_loss(w,x,K,Lg,mu,kernel_sel,**kwargs) #  loss
        df=lambda w:grad_loss_K(w,x,K,Lg,mu,kernel_sel,**kwargs) # grad  loss
    
    g=reg_l1
    prox_g=lambda x,lamb, **kw  : prox_l1(np.maximum(0,x),lamb,**kw)
    
    w0=np.ones(x.shape[1])

    if keepX == None:
        w, log = fmin_prox(f, df, g, prox_g, w0, lambd=reg, **params)
    else:
        # warm restart until number of variables ok
        val_regs =  np.logspace(-5, 2, num=nstep, endpoint=True, base=10.0)
        for i, reg in enumerate(val_regs):
            neww, log = fmin_prox(f, df, g, prox_g, w0, lambd=reg, **params)
            if (len(which(neww>0)) <= keepX):
                w = w0
                break
            else:
                w = neww
            w0 = neww
    index = sorted(range(len(w)), key=lambda k: w[k], reverse=True)
    if keepX is not None:
        index = index[:keepX]
    return (index)
####
