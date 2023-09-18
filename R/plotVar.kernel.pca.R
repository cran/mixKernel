#' Plot importance of variables in kernel PCA
#' 
#' Provides a representation of variable importance in kernel PCA.
#'
#' @details
#' \code{plotVar.kernel.pca} produces a barplot for each block. The variables for which the 
#' importance has been computed with \code{\link{kernel.pca.permute}} are 
#' displayed. The representation is limited to the \code{ndisplay} most important 
#' variables.
#'
#' @param object : a kernel.pca object returned by \code{\link{kernel.pca}}.
#' @param blocks a numerical vector indicating the block variables to display.
#' @param ndisplay integer. The number of important variables per blocks shown in 
#' the representation. Default: \code{5}.
#' @param ncol integer. Each block of variables is displayed in a separate 
#' subfigure. \code{ncol} sets the number of columns for the global figure. 
#' Default: \code{2}.
#' @param ... external arguments.
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Crone L. and Crosby D. (1995). Statistical applications of a metric on subspaces
#' to satellite meteorology. \emph{Technometrics}, \bold{37}(3), 324-328.
#' @seealso \code{\link{kernel.pca}}, \code{\link{kernel.pca.permute}}
#' @export
#' @examples
#' data(TARAoceans)
#' 
#' # compute one kernel for the psychem dataset
#' phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
#' # perform a KPCA
#' kernel.pca.result <- kernel.pca(phychem.kernel)
#' # compute importance for all variables in this kernel
#' kernel.pca.result <- kernel.pca.permute(kernel.pca.result, phychem = colnames(TARAoceans$phychem))
#' 
#' \dontrun{plotVar.kernel.pca(kernel.pca.result, ndisplay = 10)}
#' 
plotVar.kernel.pca <- function(object, blocks = unique(object$cc.blocks), 
                               ndisplay = 5, ncol = 2, ...) {
  
  df.values <- df.variables <- df.blocks <- vector(length = ndisplay * length(blocks))
  i <- 1
  for (block in blocks) {
    values.block <- object$cc.distances[which(object$cc.blocks == block)]
    variables.block <- object$cc.variables[which(object$cc.blocks == block)]
    ordered.ids <- order(values.block, decreasing = TRUE)
    next.i <- i + ndisplay - 1
    df.values[i:next.i] <- values.block[ordered.ids][1:ndisplay]
    df.variables[i:next.i] <- variables.block[ordered.ids][1:ndisplay]
    df.blocks[i:next.i] <- rep(block, ndisplay)
    i <- i + ndisplay
  }

  df <- data.frame("variables" = df.variables, "values" = df.values,
                   "blocks" = df.blocks)
  df$variables <- reorder(df$variables, -df$values)

  ggplot(df, aes_string(x="variables", y = "values", fill="blocks")) +
    geom_bar(stat = "identity") + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + ylab("") + xlab("") +
    facet_wrap(~ blocks, ncol = ncol, scales = "free_x") + 
    theme(legend.position = "none")
  
}
