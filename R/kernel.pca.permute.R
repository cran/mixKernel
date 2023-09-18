#' Assess variable importance
#' 
#' Assess importance of variables on a given PC component by computing the 
#' Crone-Crosby distance between original sample positions and sample positions 
#' obtained by a random permutation of the variables.
#'
#' @details
#' \code{plotVar.kernel.pca} produces a barplot for each block. The variables 
#' for which the importance has been computed with 
#' \code{\link{kernel.pca.permute}} are displayed. The representation is limited 
#' to the \code{ndisplay} most important variables.
#'
#' @param kpca.result a kernel.pca object returned by the
#' \code{\link{kernel.pca}} function.
#' @param ncomp integer. Number of KPCA components used to compute the 
#' importance. Default: \code{1}.
#' @param ... list of character vectors. The parameter name must be the kernel 
#' name to be considered for permutation of variables. Provided vectors length 
#' has to be equal to the number of variables of the input dataset. A kernel is 
#' performed on each unique variables values. Crone-Crosby distances are 
#' computed on each KPCA performed on resulted kernels or meta-kernels and can 
#' be displayed using the \code{\link{plotVar.kernel.pca}}.
#' @param directory character. To limit computational burden, this argument 
#' allows to store / read temporary computed kernels.
#' 
#' @return \code{kernel.pca.permute} returns a copy of the input 
#' \code{kpca.result} results and add values in the three entries: 
#' \code{cc.distances}, \code{cc.variables} and \code{cc.blocks}.
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @references Mariette J. and Villa-Vialaneix N. (2018). Unsupervised multiple 
#' kernel learning for heterogeneous data integration. \emph{Bioinformatics}, 
#' \bold{34}(6), 1009-1015. DOI: \doi{10.1093/bioinformatics/btx682}
#' 
#' Crone L. and Crosby D. (1995). Statistical applications of a metric on 
#' subspaces to satellite meteorology. \emph{Technometrics}, \bold{37}(3), 
#' 324-328.
#' @seealso \code{\link{compute.kernel}}, \code{\link{kernel.pca}}
#' @export
#' @examples
#' data(TARAoceans)
#' 
#' # compute one kernel for the psychem dataset
#' phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
#' # perform a KPCA
#' kernel.pca.result <- kernel.pca(phychem.kernel)
#' 
#' # compute importance for all variables in this kernel
#' kernel.pca.result <- kernel.pca.permute(kernel.pca.result, 
#'                                         phychem = colnames(TARAoceans$phychem))
#' 
kernel.pca.permute = function(kpca.result, ncomp = 1, ..., directory = NULL) {
  
  #-- checking general input parameters --#
  
  permutations <- list(...)
  blocks <- names(permutations)
  
  if (!"kernel.pca" %in% class(kpca.result)) {
    stop(paste0("'kpca.result' should be an instance of 'kernel.pca' object", 
                " returned by the kernel.pca function."), call. = FALSE)
  }
  
  if (length(blocks) == 0) {
    stop("At least permutation variables of 1 block should be provided.", 
         call. = FALSE)
  }
  
  # it the kpca has been performed on a meta kernel
  is.metakernel <- "metaKernel" %in% class(kpca.result$kernel)

  if (is.metakernel) {
    # test if there is no kidentity kernels, i.e. kernel with X = NULL
    sapply(blocks, function(b.label) {
      if (is.null(kpca.result$kernel$X[[b.label]]$X)) {
        stop(paste0("No permutation can be done on block '", b.label, 
                    "' as there is no dataset available for this block in the",
                    " 'kpca.result' object."), call. = FALSE)
      }
    })
    sapply(blocks, function(b.label) {
      if (!b.label %in% names(kpca.result$kernel$X)) {
        stop(paste0("Block '", b.label, "' does not exists as a block in the ",
                    "'kpca.result' object."), call. = FALSE)
      }
    })
    sapply(blocks, function(b.label) {
      if (!ncol(kpca.result$kernel$X[[b.label]]$X) == length(permutations[[b.label]])) {
        stop(paste0("'permutations' vector length for block '", b.label, 
                    "' should be equal to the number of variables in '", 
                    b.label, "' dataset: ", 
                    ncol(kpca.result$kernel$X[[b.label]]$X), "."), call. = FALSE)
      }
    })
  } else {
    # test if there its a kidentity kernel, i.e. kernel with X = NULL
    if (is.null(kpca.result$kernel$X)) {
      stop(paste0("No permutation can be done as no dataset is available in ",
                  "'kpca.result' object."), call. = FALSE)
    }
    if (length(blocks) > 1) {
      stop(paste0("Only 1 'permutations' vector should be provided as the ",
                  "kernel.pca has been performed on a single kernel."), 
           call. = FALSE)
    }
  	if (!ncol(kpca.result$kernel$X) == length(permutations[[1]])) {
  	  stop(paste0("'permutations' vector length should be equal to the number ",
  	              "of variables in the input dataset: ", 
  	              ncol(kpca.result$kernel$X), "."), call. = FALSE)
  	}
  }
  
  if (!is.null(directory) && !file.exists(file.path(directory))) {
    stop(paste0("directory '", directory, "' does not exist."), call. = FALSE)
  }
  if (ncomp > kpca.result$ncomp) {
  	ncomp <- kpca.result$ncomp
  }
  
  #-- make permutations ------------------#

  all.pdist <- all.variable <- all.block <- c()

  for (block in blocks) {

    if (is.metakernel) {
      K <- kpca.result$kernel$X[[block]]
    } else {
      K <- kpca.result$kernel
    }
      
    for (variable in as.vector(unique(permutations[[block]]))) {
  
      if (!is.na(variable)) {
  
        if (!variable %in% kpca.result$cc.variables) {
  
          # permutes values of variable
          data.permuted <- K$X
          for (col in which(permutations[[block]] == variable)) {
            data.permuted[, col] <- sample(data.permuted[, col], replace = FALSE)
          }
  
          # if the directory parameter is set, try to load the permuted kernel
          f.path <- file.path(directory, paste0(variable, ".rds"))
          if (!is.null(directory) && file.exists(f.path)) {
            kernel.permuted <- readRDS(f.path)
          } else {
            # recompute the kernel
  
            compute.args <- K$kernel.args
            compute.args$X <- data.permuted
            compute.args$kernel.func <- K$kernel.func
            kernel.permuted <- do.call("compute.kernel", compute.args)
  
            if (!is.null(directory)) {
              saveRDS(kernel.permuted, f.path)
            }
          }
  
      	if (is.metakernel) {
  
          # change the considered kernel
          all.kernels <- kpca.result$kernel$X
          all.kernels[[block]] <- kernel.permuted
  
          all.kernels.scaled <- lapply (all.kernels, function(x) {
            x.cosinus <- sweep(sweep(x$kernel, 2, sqrt(diag(x$kernel)), "/"), 1, 
                               sqrt(diag(x$kernel)), "/")
            t(t(x.cosinus - colSums(x.cosinus) / nrow(x.cosinus)) - 
                rowSums(x.cosinus) / nrow(x.cosinus)) + sum(x.cosinus) / 
              nrow(x.cosinus)^2
          })
  
          # compute the new meta.kernel with permuted variables
          meta.kernel.permuted <- lapply(as.list(1:length(all.kernels.scaled)), 
                                         function(x) {
            all.kernels.scaled[[x]] * kpca.result$kernel$weights[x]
          })
          meta.kernel.permuted <- Reduce("+", meta.kernel.permuted)
  
          kernel.tmp <- list(kernel = meta.kernel.permuted,
                             weights = kpca.result$kernel$weights)
          class(kernel.tmp) <- c("kernel", "metaKernel")
  
      	} else {
            kernel.tmp <- kernel.permuted
      	}
  
      	# re performe the KPCA
          kpca.permuted <- kernel.pca(kernel.tmp, ncomp = kpca.result$ncomp)
  
          # compute the crone and crosby distance
          pdist <- Pdist(list(O2P(kpca.result$x[,1:ncomp]),
                              O2P(kpca.permuted$x[,1:ncomp])), 
                         weights = "constant")
          
          all.pdist <- c(all.pdist, pdist)
          all.variable <- c(all.variable, variable)
          all.block <- c(all.block, block)
  
        } else {
          stop(paste0("'", variable, "'", "variable has already been permuted",
                      " for this kernel.pca object!"))
        }
  
      }
    }
    
  }
    

  kpca.result$cc.distances <- c(kpca.result$cc.distances, all.pdist)
  kpca.result$cc.variables <- c(kpca.result$cc.variables, all.variable)
  kpca.result$cc.blocks <- c(kpca.result$cc.blocks, all.block)

  return(invisible(kpca.result))
  
}