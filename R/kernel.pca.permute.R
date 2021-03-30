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