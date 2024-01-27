#' Kernel Principal Components Analysis
#' 
#' Performs a kernel PCA.
#'
#' @param K a kernel object obtained using either \code{compute.kernel} or
#' \code{combine.kernels}.
#' @param ncomp integer. Indicates the number of components to return..
#' 
#' @return \code{kernel.pca} returns an object of classes \code{"kernel.pca"} 
#' and \code{"pca"}, which is a list containing the following entries: 
#'   \item{ncomp}{: the number of principal components;}
#'   \item{X}{: the input kernel matrix;} 
#'   \item{kernel}{: the input kernel object provided by the user;}
#'   \item{sdev}{: the singular values (square root of the eigenvalues);} 
#'   \item{rotation}{: the matrix of variable loadings (\emph{i.e.}, a matrix 
#'   whose columns contain the eigenvectors);}
#'   \item{loadings}{: same as 'rotation' to keep the mixOmics spirit;}
#'   \item{x}{: same as 'rotation' to keep the mixOmics spirit;}
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Scholkopf B., Smola A. and Muller K.R. (1998) Nonlinear component 
#' analysis as a kernel eigenvalue problem. \emph{Neural Computation}, 
#' \bold{10}, 1299-1319.
#' @seealso \code{\link{compute.kernel}}, \code{\link{combine.kernels}}
#' @export
#' @examples
#' data(TARAoceans)
#' phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
#' kernel.pca.result <- kernel.pca(phychem.kernel, ncomp = 3)
#' 
kernel.pca <- function(K, ncomp = nrow(K$kernel)) {
  
  #-- checking general input parameters --#
  if (!inherits(K, "kernel")) {
    stop("'K' should be an instance of class 'kernel'.", call. = FALSE)
  }
  
  X <- K$kernel
  if (!is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp)) {
    stop("invalid value for 'ncomp'.", call. = FALSE)
  }
  if (ncomp > nrow(X)) {
    stop("use smaller 'ncomp'", call. = FALSE)
  }
  
  samples.names <- rownames(X)
  
  # test if the input matrix is symmetric and positive semidefinite
  if (!is.kernel(X)) {
    stop("Kernel matrix has to be symmetric and positive semidefinite")
  }
  
  #-- center the input kernel matrix -----#
  m <- dim(X)[1]
  X.centered <- t(t(X - colSums(X) / m) -  rowSums(X) / m) + sum(X) / m^2

  #-- compute eigenvectors ---------------#
  # eigen.result <- eigen(X.centered/m, symmetric=TRUE)
  eigen.result <- eigen(X.centered, symmetric = TRUE)
  # tolerance is used to solve numeric instabilities
  tolerance=10^(-10)
  eigen.result$values <- round(eigen.result$values, -log10(tolerance))
  pcv <- t(t(eigen.result$vectors[,1:ncomp]) / sqrt(eigen.result$values[1:ncomp]))
  colnames(pcv) <- c(1:ncomp)
  rownames(pcv) <- samples.names
  eig <- eigen.result$values[1:ncomp]
  
  cl <- match.call()
  cl[[1]] <- as.name('kernel.pca')

  result <- list(call         = cl, 
                 X            = K$kernel,
                 ncomp        = ncomp,	
                 rotation     = NULL,
                 kernel       = K,
                 sdev         = eig,
                 names        = list(X = rownames(pcv), sample = rownames(pcv)),
                 loadings     = list(X = NULL),
                 cc.distances = c(),
                 cc.variables = c(),
                 cc.blocks    = c(),
                 x            = pcv,
                 variates     = list(X = pcv))
  
  class(result) <- c("kernel.pca", "pca")
  
  # calcul explained variance
  # explX=explained_variance(X, result$variates$X, ncomp)
  # result$explained_variance=explX
  result$prop_expl_var$X <- eigen.result$values / sum(eigen.result$values)
  result$cum.var = cumsum(result$prop_expl_var$X)
  
  return(invisible(result))
}
