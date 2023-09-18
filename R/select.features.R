if(getRversion() >= "2.15.1")  utils::globalVariables(c("ukfspy"))
if(getRversion() >= "2.15.1")  utils::globalVariables(c("kokfspy"))
#' Select important features
#' 
#' Select features using supervised or unsupervised kernel method. A 
#' supervised feature selection method is performed if \code{Y} is provided.
#'
#' @param X a numeric matrix (or data frame) used to select variables. 
#' \code{NA}s not allowed.
#' @param Y a numeric matrix (or data frame) used to select variables. 
#' \code{NA}s not allowed.
#' @param kx.func the kernel function name to use on \code{X}. Widely used 
#' kernel functions are pre-implemented, and can be directly used by setting 
#' \code{kx.func} to one of the following values: \code{"linear"}, 
#' \code{"gaussian.radial.basis"} or \code{"bray"}. Default: \code{"linear"}. If 
#' \code{Y} is provided, the kernel \code{"bray"} is not allowed.
#' @param ky.func the kernel function name to use on \code{Y}. Available 
#' kernels are: \code{"linear"}, and \code{"gaussian.radial.basis"}. Default: 
#' \code{"linear"}. This value is ignored when \code{Y} is not provided.
#' @param keepX the number of variables to select.
#' @param method the method to use. Either an unsupervised variable selection
#' method (\code{"kernel"}), a kernel PCA oriented variable selection method 
#' (\code{"kpca"}) or a structure driven variable selection selection 
#' (\code{"graph"}). Default: \code{"kernel"}.
#' @param lambda the penalization parameter that controls the trade-off between 
#' the minimization of the distorsion and the sparsity of the solution 
#' parameter.
#' @param n_components how many principal components should be used with method
#' \code{"kpca"}. Required with method \code{"kpca"}. Default: \code{2}.
#' @param Lg the Laplacian matrix of the graph representing relations between 
#' the input dataset variables. Required with method \code{"graph"}.
#' @param mu the penalization parameter that controls the trade-off between the
#' the distorsion and the influence of the graph. Default: \code{1}.
#' @param max_iter the maximum number of iterations. Default: \code{100}.
#' @param nstep the number of values used for the regularization path. Default: 
#' \code{50}. 
#' @param ... the kernel function arguments. In particular, 
#' \code{sigma}(\code{"gaussian.radial.basis"}): double. The inverse kernel 
#' width used by \code{"gaussian.radial.basis"}.
#' 
#' @return \code{ukfs} returns a vector of sorted selected features indexes.
#' 
#' @author Celine Brouard <celine.brouard@@inrae.fr>
#' Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Brouard C., Mariette J., Flamary R. and Vialaneix N. (2022). 
#' Feature selection for kernel methods in systems biology. \emph{NAR Genomics
#' and Bioinformatics}, \bold{4}(1), lqac014. DOI: \doi{10.1093/nargab/lqac014}.
#' @seealso \code{\link{compute.kernel}}
#' @export select.features
#' @examples
#' ## These examples require the installation of python modules
#' ## See installation instruction at: http://mixkernel.clementine.wf
#'
#' data("Koren.16S")
#' \dontrun{
#'  sf.res <- select.features(Koren.16S$data.raw, kx.func = "bray", lambda = 1,
#'                            keepX = 40, nstep = 1)
#'  colnames(Koren.16S$data.raw)[sf.res]
#' }
#'
#' data("nutrimouse")
#' \dontrun{
#'  grb.func <- "gaussian.radial.basis"
#'  genes <- center.scale(nutrimouse$gene)
#'  lipids <- center.scale(nutrimouse$lipid)
#'  sf.res <- select.features(genes, lipids, kx.func = grb.func, 
#'                            ky.func = grb.func, keepX = 40)
#'  colnames(nutrimouse$gene)[sf.res]
#' }
#' 
select.features <- function(X, Y=NULL, kx.func=c("linear", "gaussian.radial.basis", "bray"),
                            ky.func=c("linear", "gaussian.radial.basis"), keepX=NULL,
                            method=c("kernel", "kpca", "graph"),
                            lambda=NULL, n_components=2, Lg=NULL,
                            mu=1, max_iter=100, nstep=50, ...) {

  source_python(file.path(system.file(package = "mixKernel"), "python", "ukfs.py"))
  source_python(file.path(system.file(package = "mixKernel"), "python", "kokfs.py"))
  kernel.args <- list(...)
  kx.func <- match.arg(kx.func)
  ky.func <- match.arg(ky.func)
  method <- match.arg(method)
  
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  if (method == "graph" && is.null(Lg)) {
  	stop("'Lg' is required with method 'graph'")
  }
  if (method == "graph") {
    if (is.null(dim(Lg))) {
      stop("'Lg' should be a matrix")
    }
    if (dim(Lg)[1] != dim(Lg)[2] || dim(Lg)[1] != dim(X)[2]) {
  	  stop(paste0("'Lg' should be a squared matrix of size "), dim(X)[2], ".")
  	}
    if (any(is.na(Lg))) {
  	  stop("'Lg' should not contains NA values")
  	}
  }
  if (!is.numeric(n_components)) {
    stop("'n_components' must be numeric")
  }
  if (!is.numeric(mu)) {
    stop("'mu' must be numeric")
  }
  if (is.null(keepX) && is.null(lambda)) {
    stop("One of parameters 'keepX' or 'lambda' is required")
  }
  if (!is.null(keepX) && !is.numeric(keepX)) {
    stop("'keepX' must be numeric")
  }
  if (!is.null(keepX) && keepX>dim(X)[2]) {
    stop(paste0("'keepX' must be numeric and inferior to the number of variables, ", dim(X)[2], "."))
  }
  if (!is.null(lambda) && !is.numeric(lambda)) {
    stop("'lambda' must be numeric")
  }
  if (!is.numeric(max_iter)) {
    stop("'max_iter' must be numeric")
  }
  if (!is.numeric(nstep)) {
    stop("'nstep' must be numeric")
  }

  # force X, Y and Lg to be a matrix for numpy
  X <- as.matrix(X)
  if (!is.null(Y)) {
    Y <- as.matrix(Y)
  }
  if (!is.null(Lg)) {
  	Lg <- as.matrix(Lg)
  }

  if (is.null(Y)) {
    w <- ukfspy(X, kx.func, method, keepX, lambda, n_components, Lg, mu, max_iter, nstep, kernel.args)
  } else {
    w <- kokfspy(t(X), t(Y), kx.func, ky.func, keepX, nstep)
  }
  return (w)

}