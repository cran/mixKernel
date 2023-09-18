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

#############################################################################################################
#
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


#############################################################################################################
# linear kernel
#############################################################################################################

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
  } else {
    X.scaled <- X
  }
  
  similarities <- tcrossprod(X.scaled)
  
  return(invisible(similarities))
  
}


#############################################################################################################
# gaussian.radial.basis kernel
#############################################################################################################

gaussian.radial.basis <- function(X, scale = TRUE, sigma = 1) {
  
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
  } else {
    X.scaled <- X
  }

  n <- dim(X.scaled)[1]
  dota <- rowSums(X.scaled * X.scaled) / 2
  
  similarities <- tcrossprod(X.scaled)
  for (i in 1:n) {
    similarities[i, ] <- exp(2 * sigma * (similarities[i,] - dota - rep(dota[i], n)))
  }
  return(invisible(similarities))
  
}


#############################################################################################################
# abundance kernel
#############################################################################################################

abundance <- function(X, method = c("bray", "manhattan", "euclidean", "canberra", 
                                    "kulczynski", "jaccard", "gower", "altGower", 
                                    "morisita", "horn",	"mountford", "raup" , 
                                    "binomial", "chao", "cao")) {

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
  
  dissimilarities <- as.matrix(vegdist(X, method=method))
    
  # convert dissimilarities into similarities
  tolerance <- 10^(-10)
  similarities <- -.5* (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities)) %*%
    dissimilarities %*% (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities))
  similarities <- round(similarities, -log10(tolerance))
  
  return(invisible(similarities))

}


#############################################################################################################
# phylogenetic kernel
#############################################################################################################

phylogenetic <- function(X, method = c("wunifrac", "unifrac"), 
                         phylogenetic.tree = NULL) {

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
  	stop("'phylogenetic.tree' is required when computing a 'phylogenetic' kernel.", 
  	     call. = FALSE)
  }
  if (!inherits(phylogenetic.tree, "phylo")) {
  	stop("'phylogenetic.tree' should be an instance of 'phylo-class' as defined",
  	     " in the ape-package.", call. = FALSE)
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


#############################################################################################################
# poisson kernel
#############################################################################################################

null.model <- function(x, normalization = c("deseq", "mle", "quantile")) {
  
  normalization <- match.arg(normalization)
  rowsum <- rowSums(x)
  colsum <- colSums(x)
  mle <- outer(rowsum, colsum, "*") / sum(rowsum)
  
  if(normalization=="mle"){
    return(list(n = mle, sizes = rowSums(x) / sum(x)))
  } else if (normalization=="quantile"){
    # This is quantile normalization idea of Bullard et al 2010
    # quantile-normalize using 75th quantile of observed counts 
    # for each sample, excluding zero-counts
    sample.qts <- apply(x, 1, quantile, .75)
    # Don't wait to standardize by 0... min allowed is 1
    sample.qts <- pmax(sample.qts, 1) 
    sample.qts <- sample.qts / sum(sample.qts)
    fit <- outer(sample.qts, colsum, "*")
    return(list(n = fit, sizes = sample.qts))
  } else if (normalization=="deseq"){
    #Trying to implement idea from Anders and Huber Genome Biology 2010 paper.
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))
    sizes <- apply(counts, 2, function(cnts) median((cnts / geomeans)[geomeans > 0]))
    rawsizestr <- sizes
    sizes <- sizes / sum(sizes)
    fit <- outer(sizes, colsum, "*")
    return(list(n = fit, sizes = sizes, geomeans = geomeans, rawsizestr = rawsizestr))
  }
  
}
goodness.of.fit <- function(x, normalization) {
  ns <- null.model(x, normalization = normalization)$n
  return(sum(((x - ns)^2) / ns, na.rm = TRUE))
}
find.best.transform <- function(x) {
  alphas <- seq(.01, 1, len=50)
  gof <- rep(NA, length(alphas))
  for(alpha in alphas){
    gof[alphas==alpha] <- goodness.of.fit(x^alpha, normalization="mle")
  }
  return(alphas[which.min(abs(gof - (nrow(x) - 1) * (ncol(x) - 1)))])
}

poisson <- function(X, normalization = c("deseq", "mle", "quantile")) {
  
  beta <- 1
  normalization <- match.arg(normalization)
  alpha <- find.best.transform(X)
  X <- X^alpha
  dd <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  for(i in 2:nrow(dd)){
    Xi <- X[i,]
    for(j in 1:(i - 1)){
      Xj <- X[j,]
      n <- null.model(X[c(i, j), ], normalization = normalization)$n
      ni <- n[1,]
      nj <- n[2,]
      di <- (Xi + beta) / (ni + beta)
      dj <- (Xj + beta) / (nj + beta)
      dd[i, j] <- sum(ni + nj - ni * di - nj * dj + Xi * log(di) + Xj * log(dj))
    }
  }
  
  dissimilarities <- as.matrix(as.dist(dd+t(dd)))
  
  # convert dissimilarities into similarities
  tolerance <- 10^(-10)
  similarities <- -.5* (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities)) %*%
    dissimilarities %*% (diag(1, nrow(dissimilarities)) - 1 / nrow(dissimilarities))
  similarities <- forceSymmetric(round(similarities, -log10(tolerance)))
  
  return(invisible(as.matrix(similarities)))
  
}
