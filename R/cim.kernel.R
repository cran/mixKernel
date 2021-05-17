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


cim.kernel <- function(..., scale = TRUE, 
                       method = c("circle", "square", "number", 
                                  "shade", "color", "pie")) {
  
  #-- checking general input parameters --#
  
  K <- list(...)
  method <- match.arg(method)
  
  if (!is.logical(scale)) {
    stop("scale must be either TRUE or FALSE")
  }
  if (!is.list(K)) {
    stop("K must be a list")
  }
  # check names on K are unique
  if (length(unique(names(K))) != length(K)) {
    stop("Each block of 'K' must have a unique name.")
  }
  # check all inputs are kernels
  are.kernels <- unlist(lapply(K, function (k) {
    "kernel" %in% class(k)
  }))
  if (!all(are.kernels)) {
    stop("Each block of 'K' must be a kernel.")
  }
  # check all kernels have been computed on the same number of observation
  if (length(unique(unlist(lapply(K, function(x) { ncol(x$kernel) })))) != 1) {
    stop("Unequal number of observations among the kernels of 'K'")
  }
  
  #-- cosinus scaling --------------------#
  
  if (scale) {
    K.scaled <- lapply (K, function(x) {
      x.cosinus <- sweep(sweep(x$kernel, 2, sqrt(diag(x$kernel)), "/"),
                         1, sqrt(diag(x$kernel)), "/")
      t(t(x.cosinus - colSums(x.cosinus) / nrow(x.cosinus)) - 
          rowSums(x.cosinus) / nrow(x.cosinus)) + 
        sum(x.cosinus) / nrow(x.cosinus)^2
    })
  } else {
    K.scaled <- K
  }
  
  
  #-- STATIS similarities ----------------#
  
  similarities <- outer(1:length(K.scaled), 1:length(K.scaled), 
                        FUN = Vectorize(function(i, j) {
    out <- tr(K.scaled[[i]] %*% K.scaled[[j]]) 
    out <- out / (norm(K.scaled[[i]], type="F") * norm(K.scaled[[j]], type="F"))
    return(out)
  }))
  
  rownames(similarities) <- colnames(similarities) <- names(K)
  corrplot(similarities, type = "full", tl.col = "black", 
           tl.srt = 45, method = method)
  # corrplot.mixed(similarities, upper = method, cl.lim = c(0, 1))
  return(invisible(similarities))
  
}

