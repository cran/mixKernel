#' Compute a kernel
#' 
#' Compute a kernel from a given data matrix.
#'
#' @param X a numeric matrix (or data frame) used to compute the kernel. 
#' \code{NA}s not allowed.
#' @param kernel.func the kernel function to use. This parameter can be set to 
#' any user defined kernel function. Widely used kernel functions are 
#' pre-implemented, that can be used by setting \code{kernel.func} to one of the
#'  following strings: \code{"kidentity"}, \code{"abundance"}, \code{"linear"}, 
#' \code{"gaussian.radial.basis"}, \code{"poisson"} or \code{"phylogenetic"}. 
#' Default: \code{"linear"}.
#' @param ... the kernel function arguments. Valid parameters for 
#' pre-implemented kernels are: \itemize{
#'     \item \code{phylogenetic.tree} (\code{"phylogenetic"}): an instance of 
#'     phylo-class that contains a phylogenetic tree (required).
#'     \item \code{scale} (\code{"linear"} or \code{"gaussian.radial.basis"}): 
#'       logical. Should the variables be scaled to unit variance prior the 
#'       kernel computation? Default: \code{TRUE}.
#'     \item \code{sigma} (\code{"gaussian.radial.basis"}): double. The inverse 
#'     kernel width used by \code{"gaussian.radial.basis"}.
#'     \item \code{method} (\code{"phylogenetic"} or \code{"abundance"}): 
#'     character. Can be \code{"unifrac"} or \code{"wunifrac"} for 
#'     \code{"phylogenetic"}. Which dissimilarity to use for \code{"abundance"}:
#'     one of \code{"bray"}, \code{"euclidean"}, \code{"canberra"}, 
#'     \code{"manhattan"}, \code{"kulczynski"}, \code{"jaccard"}, 
#'     \code{"gower"}, \code{"altGower"}, \code{"morisita"}, \code{"horn"}, 
#'     \code{"mountford"}, \code{"raup"}, \code{"binomial"}, \code{"chao"} and 
#'     \code{"cao"}.
#'     \item \code{normalization} (\code{"poisson"}): character. Can be 
#'     \code{"deseq"} (more robust), \code{"mle"} (less robust) or 
#'     \code{"quantile"}.
#' }
#' @param test.pos.semidef boleean. If \code{test.pos.semidef = TRUE}, the 
#' positive semidefiniteness of the resulting matrix is checked.
#' 
#' @return \code{compute.kernel} returns an object of classes \code{"kernel"}, a
#' list that contains the following components:
#' \item{kernel}{: the computed kernel matrix.}
#' \item{X}{: the original dataset. If \code{"kidentity"}, \code{X} is set to 
#' \code{NULL}.}
#' \item{kernel.func}{: the kernel function used.}
#' \item{kernel.args}{: the arguments used to compute the kernel.} 
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Lozupone C. and Knight R. (2005). UniFrac: a new phylogenetic 
#' method for comparing microbial communities. \emph{Applied and Environmental 
#' Microbiology}, \bold{71}(12), 8228-8235.
#'
#' Lozupone C., Hamady M., Kelley S.T. and Knight R. (2007). Quantitative and 
#' qualitative beta diversity measures lead to different insights into factors 
#' that structure microbial communities. \emph{Applied and Environmental 
#' Microbiology}, \bold{73}(5), 1576-1585.
#' 
#' Witten D. (2011). Classification and clustering of sequencing data using a 
#' Poisson model. \emph{Annals of Applied Statistics}, \bold{5}(4), 2493-2518.
#' @seealso \code{\link{combine.kernels}}, \code{\link{kernel.pca}}
#' @export
#' @examples
#' data(TARAoceans)
#' pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, 
#'                                   kernel.func = "abundance")
#' 
compute.kernel <- function(X, kernel.func = "linear", ..., 
                           test.pos.semidef = FALSE) {
  
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
