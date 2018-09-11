#############################################################################################################
# Author :
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
  X.classes <- unique(unlist(lapply(X, class)))
  if ( length(X.classes) != 1 || X.classes[1] != "kernel") {
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
