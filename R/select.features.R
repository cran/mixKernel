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

if(getRversion() >= "2.15.1")  utils::globalVariables(c("ukfspy"))
if(getRversion() >= "2.15.1")  utils::globalVariables(c("kokfspy"))

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