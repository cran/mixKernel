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

compute.kernel <- function(X, kernel.func = "linear", ..., test.pos.semidef = FALSE) {
  
  #-- checking general input parameters --#
  
  kernel.args <- list(...)
  
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  if (is(kernel.func, "function")) {
  	kernel.func <- deparse(substitute(kernel.func))
  }
  if (!is.logical(test.pos.semidef)) {
    stop("test.pos.semidef must be either TRUE or FALSE")
  }
  
  if (kernel.func == "kidentity") {
    similarities <- X
    colnames(similarities) <- rownames(similarities) <- rownames(X)
    X <- NULL
  } else {
    all.kernel.args <- kernel.args
    all.kernel.args$X <- X
    similarities <- do.call(kernel.func, all.kernel.args) 
    colnames(similarities) <- rownames(similarities) <- rownames(X)
  }
  
  # test if the resulting matrix is symmetric and positive semidefinite
  if (!is.kernel(similarities, test.pos.semidef)) {
    stop("Kernel matrix has to be symmetric and positive semidefinite")
  }
  
  cl <- match.call()
  cl[[1]] <- as.name('compute.kernel')
  result <- (list(call = cl, kernel = similarities, X = X,
                  kernel.func = kernel.func, kernel.args = kernel.args))
  class(result) <- c("kernel", "matrix", "array")
  
  return(invisible(result))
  
}
