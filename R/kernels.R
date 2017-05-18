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


# Checks whether a matrix is a kernel
is.kernel <- function (K, test.pos.semidef = FALSE) {
  
  if (!is.matrix(K)) {
    K <- as.matrix(K)
  }
  
  if (test.pos.semidef) {
    eval <- eigen(K, only.values = TRUE, symmetric = TRUE)$values
    tolerance <- 10^(-10)
    
    if (is.complex( eval )) {
      stop("Kernel matrix has complex eigenvalues")
    }
    
    if (isSymmetric(K)) {
      return(all(eval[which(eval > tolerance)] > 0))
    } else {
      return(FALSE)
    }
  } else {
    return(isSymmetric(K))
  }
  
}


linear <- function(X, scale = TRUE) {

  #-- checking general input parameters --#
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  if (! is.logical(scale)) {
  	stop("'scale' should be a logical value")
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (! is.matrix(X) || is.character(X)) {
    stop("'X' must be a numeric matrix.", call. = FALSE)
  }
  if (any(apply(X, 1, is.infinite))) {
    stop("infinite values in 'X'.", call. = FALSE)
  }
  if (scale) {
	  X.scaled <- scale(X)
  }
  
  similarities <- tcrossprod(X.scaled)
  
  return(invisible(similarities))
  
}


abundance <- function(X, method = c("bray", "manhattan", "euclidean", "canberra",
                                           "kulczynski", "jaccard", "gower", "altGower", "morisita", 
                                           "horn",	"mountford", "raup" , "binomial", "chao", "cao")) {

  #-- checking general input parameters --#
  method <- match.arg(method)
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  if (is.data.frame(X)) {
    X = as.matrix(X)
  }
  if (!is.matrix(X) || is.character(X)) {
    stop("'X' must be a numeric matrix.", call. = FALSE)
  }
  if (any(apply(X, 1, is.infinite))) {
    stop("infinite values in 'X'.", call. = FALSE)
  }
  
  X.OTU <- otu_table(X, taxa_are_rows = FALSE)
  X.phyloseq <- phyloseq(X.OTU)
  
  dissimilarities <- as.matrix(distance(X.phyloseq, method=method))
    
  # convert dissimilarities into similarities
  tolerance <- 10^(-10)
  similarities <- -.5* (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities)) %*%
    dissimilarities %*% (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities))
  similarities <- round(similarities, -log10(tolerance))
  
  return(invisible(similarities))

}


phylogenetic <- function(X, method = c("wunifrac", "unifrac"), phylogenetic.tree = NULL) {

  #-- checking general input parameters --#
  method <- match.arg(method)
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || is.character(X)) {
    stop("'X' must be a numeric matrix.", call. = FALSE)
  }
  if (any(apply(X, 1, is.infinite))) {
    stop("infinite values in 'X'.", call. = FALSE)
  }
  if (is.null(phylogenetic.tree)) {
  	stop("'phylogenetic.tree' is required when computing a 'phylogenetic' kernel.", call. = FALSE)
  }
  if (!class(phylogenetic.tree) == "phylo") {
  	stop("'phylogenetic.tree' should be an instance of 'phylo-class' as defined in the ape-package.", call. = FALSE)
  }
  
  X.OTU <- otu_table(X, taxa_are_rows=FALSE)
  X.phyloseq <- phyloseq(X.OTU, phylogenetic.tree)  
  
  dissimilarities <- as.matrix(distance(X.phyloseq, method = method))

  # convert dissimilarities into similarities
  tolerance <- 10^(-10)
  similarities <- -.5* (diag(1, nrow(dissimilarities))- 1 / nrow(dissimilarities)) %*%
    dissimilarities %*% (diag(1, nrow(dissimilarities))- 1 / nrow(dissimilarities))
  similarities <- round(similarities, -log10(tolerance))
  
  return(invisible(similarities))

}
