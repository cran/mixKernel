#' Combine multiple kernels into a meta-kernel
#' 
#' Compute multiple kernels into a single meta-kernel
#'
#' @param ... list of kernels (called 'blocks') computed on different datasets 
#' and measured on the same samples.
#' @param scale boleean. If \code{scale = TRUE}, each block is standardized to 
#' zero mean and unit variance and cosine normalization is performed on the 
#' kernel. Default: \code{TRUE}.
#' @param method character. Which method should be used to compute the 
#' meta-kernel. Default: \code{"full-UMKL"}.
#' @param knn integer. If \code{method = "sparse-UMKL"} or
#' \code{method = "full-UMKL"}, number of neighbors used to get a proxy of the 
#' local topology of the datasets from each kernel. Default: \code{5}.
#' @param rho integer. Parameters for the augmented Lagrangian method. Default: 
#' \code{20}.
#' 
#' @return \code{combine.kernels} returns an object of classes \code{"kernel"} 
#' and \code{"metaKernel"}, a list that contains the following components: 
#'   \item{kernel}{: the computed meta-kernel matrix;}
#'   \item{X}{: the dataset from which the kernel has been computed, as given by
#'   the function \code{\link{compute.kernel}}. Can be \code{NULL} if a kernel
#'   matrix was passed to this function;}
#'   \item{weights}{: a vector containing the weights used to combine the 
#'   kernels.} 
#'   
#' @details
#' The arguments \code{method} allows to specify the Unsupervised Multiple
#' Kernel Learning (UMKL) method to use: \itemize{
#'   \item \code{"STATIS-UMKL"}: combines input kernels into the best 
#'   consensus of all kernels;
#'   \item \code{"full-UMKL"}: computes a kernel that minimizes the distortion 
#'   between the meta-kernel and the k-NN graphs obtained from all input 
#'   kernels;
#'   \item \code{"sparse-UMKL"}: a sparse variant of the \code{"full-UMKL"} 
#'   approach.
#' }
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Mariette J. and Villa-Vialaneix N. (2018). Unsupervised multiple 
#' kernel learning for heterogeneous data integration . \emph{Bioinformatics}, 
#' \bold{34}(6), 1009-1015. DOI: \doi{10.1093/bioinformatics/btx682}.
#' @seealso \code{\link{compute.kernel}}, \code{\link{kernel.pca}}
#' @export
#' @examples
#' data(TARAoceans)
#' 
#' # compute one kernel per dataset
#' phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
#' pro.phylo.kernel <- compute.kernel(TARAoceans$pro.phylo, kernel.func = "abundance")
#' pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, kernel.func = "abundance")
#' 
#' # compute the meta kernel
#' meta.kernel <- combine.kernels(phychem = phychem.kernel,
#'                                pro.phylo = pro.phylo.kernel,
#'                                pro.NOGs = pro.NOGs.kernel, 
#'                                method = "full-UMKL")
#' 
combine.kernels <- function(..., scale = TRUE, 
                            method = c("full-UMKL", "STATIS-UMKL", "sparse-UMKL"),
                            knn = 5, rho = 20) {
  
  #-- checking general input parameters --#
  X <- list(...)
  method <- match.arg(method)
  
  if (!is.logical(scale)) {
    stop("scale must be either TRUE or FALSE")
  }
  
  if (!is.list(X)) {
    stop("X must be a list")
  }
    
  # check names on X are unique
  if (length(unique(names(X))) != length(X)) {
    stop("Each block of 'X' must have a unique name.")
  }
  
  # check all inputs are kernels
  are.kernels <- unlist(lapply(X, function (x) {
    "kernel" %in% class(x)
  }))
  if (!all(are.kernels)) {
    stop("Each block of 'X' must be a kernel.")
  }
  # check all kernels have been computed on the same number of observation
  if (length(unique(unlist(lapply(X, function(x) { ncol(x$kernel) })))) != 1) {
    stop("Unequal number of observations among the kernels of 'X'")
  }
  if (!is.numeric(knn) || knn < 1 || !is.finite(knn)) {
    stop("invalid value for 'knn'.", call. = FALSE)
  }
  if (!is.numeric(rho) || rho < 1 || !is.finite(rho)) {
    stop("invalid value for 'rho'.", call. = FALSE)
  }
  #-- cosinus scaling --------------------#
  if (scale) {
    X.scaled <- lapply (X, function(x) {
      x.cosinus <- sweep(sweep(x$kernel, 2, sqrt(diag(x$kernel)), "/"), 
                         1, sqrt(diag(x$kernel)), "/")
      t(t(x.cosinus - colSums(x.cosinus) / nrow(x.cosinus)) - rowSums(x.cosinus) / 
          nrow(x.cosinus)) + sum(x.cosinus) / nrow(x.cosinus)^2
    })
  } else {
    X.scaled <- X
  }
  
  
  #-- UMKL approaches --------------------#
  beta <- 1 / length(X.scaled)
  if (method == 'STATIS-UMKL') {
    
    similarities <- outer(1:length(X.scaled), 1:length(X.scaled), 
                          FUN = Vectorize(function(i, j) {
      tr(X.scaled[[i]] %*% X.scaled[[j]]) / (norm(X.scaled[[i]], type="F") * 
                                               norm(X.scaled[[j]], type="F"))
    }))
    weights <- eigen(similarities, symmetric = TRUE)$vectors[ ,1]
    weights <- weights / sum(weights)
    
  } else { # for both full-UMKL and sparse-UMKL methods
    
    # build the adjacency matrix considering the input kernels
    all.adjacency <- lapply(X.scaled, function(x) {
      adjacency.x <- apply(x, 1, function(row) {
        adjacency.row <- rank(row)
        adjacency.row[which(adjacency.row < length(adjacency.row)-knn)] <- 0
        adjacency.row[which(adjacency.row > 0)] <- 1
        adjacency.row
      })
      diag(adjacency.x) <- 0
      adjacency.x + t(adjacency.x)
    })
    graph.weights <- Reduce("+", all.adjacency)
    
    # C.matrix = lapply(X.scaled, function(k1){
    #   lapply(X.scaled, function(k2){
    #     sum (graph.weights * apply(k1, 1, function(rowi) {
    #       apply(k1, 1, function(rowj){
    #         sum(abs(rowi-rowj))
    #       })
    #     }) %*% apply(k2, 1, function(rowi) {
    #       apply(k2, 1, function(rowj){
    #         sum(abs(rowi-rowj))
    #       })
    #     }))
    #   })
    # })
    C.matrix <- outer(1:length(X.scaled), 1:length(X.scaled), 
                      FUN = Vectorize(function(i, j) {
      sum(graph.weights * apply(X.scaled[[i]], 1, function(rowi) {
        apply(X.scaled[[i]], 1, function(rowj) {
          sum(abs(rowi-rowj))
        })
      }) %*% apply(X.scaled[[j]], 1, function(rowi) {
        apply(X.scaled[[j]], 1, function(rowj){
          sum(abs(rowi-rowj))
        })
      }))
    }))
    C.matrix <- matrix(unlist(C.matrix), nrow=length(X.scaled))

    if (method == 'sparse-UMKL') { # if sparse-UMKL, use quadprog
      
      dvec <- rep(0, length(X.scaled))
      # cosinus normalization
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), "/"), 1, 
                          sqrt(diag(C.matrix)), "/")
      # frobenius normalization
      # C.matrix.n <- C.matrix / norm(C.matrix, type="F")
      Dmat <- 2 * C.matrix.n
      Amat <- matrix(data = 1, nrow = nrow(Dmat), ncol = 1)
      Amat <- cbind(Amat, diag(1, nrow = nrow(Dmat), ncol = nrow(Dmat)))
      bvec <- rep(0, length(X.scaled) + 1)
      bvec[1] <- 1
      weights <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)$solution
      
    } else { # else, solve the problem using ADMM
      
      k <- 1
      Z <- rep(1/sqrt(length(X.scaled)), length(X.scaled))
      Y <- rep(0,length(X.scaled))
      threshold <- 10^(-2)
      # cosinus normalization
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), "/"), 1, 
                          sqrt(diag(C.matrix)), "/")
      # frobenius normalization
      # C.matrix.n <- C.matrix / norm(C.matrix, type="F")
      Dmat <- C.matrix.n + diag(rho/2, ncol=length(X.scaled), nrow=length(X.scaled))
      Dmat <- 2 * Dmat
      Amat <- diag(1, nrow=dim(Dmat)[1], ncol=dim(Dmat)[1])
      bvec <- rep(0, length(X.scaled))
      repeat {
        # solve x using quadprog
        dvec <- rho*Z - Y
        X.ADMM <- solve.QP(Dmat, dvec, Amat, bvec, meq=0)$solution
        # project Z
        newZ <- X.ADMM / norm(matrix(X.ADMM), "2")
        if (k != 1 && norm(matrix(Z - newZ)) < threshold) break
        Z <- newZ
        # update Y
        Y <- Y + rho*(X.ADMM-Z)
        k <- k + 1
      }
      weights <- Z / norm(matrix(Z))
      
    }
  }
  
  # weights <- round(weights, -log10(1e-08))
  
  #-- compute the meta-kernel ------------#
  meta.kernel <- lapply(as.list(1:length(X.scaled)), function(x){
    X.scaled[[x]]*weights[x]
  })
  meta.kernel <- Reduce("+",meta.kernel)
  
  # output the meta-kernel, the coefficients used
  cl <- match.call()
  cl[[1]] <- as.name('combine.kernels')
  result <- list(call = cl, kernel = meta.kernel, X = X, weights = weights)
  class(result) <- c("kernel", "metaKernel")
  
  return(invisible(result))
}
