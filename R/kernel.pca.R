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


kernel.pca <- function(K, ncomp = nrow(K$kernel)) {
  
  #-- checking general input parameters --#
  if (!"kernel" %in% class(K)) {
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
  result$explained_variance <- eigen.result$values / sum(eigen.result$values)
  result$cum.var = cumsum(result$explained_variance)
  
  return(invisible(result))
}
