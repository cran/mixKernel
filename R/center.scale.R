if(getRversion() >= "2.15.1")  utils::globalVariables(c("centerscalepy"))
#' Center and scale
#' 
#' Center and scale a dataset.
#'
#' @param X a numeric matrix (or data frame) to center and scaled.
#' \code{NA}s not allowed.
#' 
#' @return \code{center.scale} returns a centered and scaled matrix.
#' 
#' @author Celine Brouard <celine.brouard@@inrae.fr>
#' Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @seealso \code{\link{compute.kernel}}, \code{\link{combine.kernels}}
#' @export
#' @examples
#' data("nutrimouse")
#' \dontrun{
#'  nutrimouse.sc <- center.scale(nutrimouse$gene)
#' }
#' 
center.scale <- function(X) {

  source_python(file.path(system.file(package = "mixKernel"), "python", "center_scale.py"))
  
  if (any(is.na(X))) {
    stop("'X' should not contain NAs")
  }
  # force X, Y and Lg to be a matrix for numpy
  X <- as.matrix(X)
  return (centerscalepy(X))

}