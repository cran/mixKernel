#' Compute and display similarities between multiple kernels
#' 
#' Compute cosine from Frobenius norm between kernels and display the 
#' corresponding correlation plot.
#'
#' @details
#' The displayed similarities are the kernel generalization of the 
#' RV-coefficient described in Lavit \emph{et al.}, 1994.
#' 
#' The plot is displayed using the \code{\link[corrplot]{corrplot}} package. 
#' Seven visualization methods are implemented: \code{"circle"} (default), 
#' \code{"square"}, \code{"number"}, \code{"pie"}, \code{"shade"} and 
#' \code{"color"}. Circle and square areas are proportional to the absolute 
#' value of corresponding similarities coefficients.
#'
#' @param ... list of kernels (called 'blocks') computed on different datasets 
#' and measured on the same samples.
#' @param scale boleean. If \code{scale = TRUE}, each block is standardized to 
#' zero mean and unit variance and cosine normalization is performed on the 
#' kernel. Default: \code{TRUE}.
#' @param method character. The visualization method to be used. Currently, 
#' seven methods are supported (see Details).
#' 
#' @return \code{cim.kernel} returns a matrix containing the cosine from 
#' Frobenius norm between kernels.
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Lavit C., Escoufier Y., Sabatier R. and Traissac P. (1994). The 
#' ACT (STATIS method). \emph{Computational Statistics and Data Analysis}, 
#' \bold{18}(1), 97-119.
#' 
#' Mariette J. and Villa-Vialaneix N. (2018). Unsupervised multiple kernel 
#' learning for heterogeneous data integration. \emph{Bioinformatics}, 
#' \bold{34}(6), 1009-1015.
#' @seealso \code{\link{compute.kernel}}
#' @export
#' @examples
#' data(TARAoceans)
#' 
#' # compute one kernel per dataset
#' phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
#' pro.phylo.kernel <- compute.kernel(TARAoceans$pro.phylo, 
#'                                    kernel.func = "abundance")
#' pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, 
#'                                   kernel.func = "abundance")
#' 
#' # display similarities between kernels
#' cim.kernel(phychem = phychem.kernel,
#'            pro.phylo = pro.phylo.kernel,
#'            pro.NOGs = pro.NOGs.kernel, 
#'            method = "square")
#' 
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

