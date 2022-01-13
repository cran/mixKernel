## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(dpi = 100, echo = TRUE, warning = FALSE, message = FALSE)

## ----load_lib-----------------------------------------------------------------
## required python modules: autograd, numpy, scipy, sklearn
## To properly install packages, run:
# install.packages("BiocManager")
# BiocManager::install("mixOmics")
# BiocManager::install("phyloseq")
# install.packages("mixKernel")
library(mixOmics)
library(mixKernel)

## ----input_data---------------------------------------------------------------
data(TARAoceans)
# more details with: ?TARAOceans
# we check the dimension of the data:
lapply(list("phychem" = TARAoceans$phychem, "pro.phylo" = TARAoceans$pro.phylo, 
            "pro.NOGs" = TARAoceans$pro.NOGs), dim)

## ----compute_kernel, echo=TRUE------------------------------------------------
phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
pro.phylo.kernel <- compute.kernel(TARAoceans$pro.phylo, kernel.func = "abundance")
pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, kernel.func = "abundance")

# check dimensions
dim(pro.NOGs.kernel$kernel)

## ----cim_kernel, fig.width=4--------------------------------------------------
cim.kernel(phychem = phychem.kernel,
           pro.phylo = pro.phylo.kernel,
           pro.NOGs = pro.NOGs.kernel, 
           method = "square")

## ----meta_kernel--------------------------------------------------------------
meta.kernel <- combine.kernels(phychem = phychem.kernel,
                               pro.phylo = pro.phylo.kernel,
                               pro.NOGs = pro.NOGs.kernel, 
                               method = "full-UMKL")

## ----KPCA---------------------------------------------------------------------
kernel.pca.result <- kernel.pca(meta.kernel, ncomp = 10)

## ----plotIndiv_PCA, fig.keep='all'--------------------------------------------
all.depths <- levels(factor(TARAoceans$sample$depth))
depth.pch <- c(20, 17, 4, 3)[match(TARAoceans$sample$depth, all.depths)]
plotIndiv(kernel.pca.result,
          comp = c(1, 2),
          ind.names = FALSE,
          legend = TRUE,
          group = as.vector(TARAoceans$sample$ocean),
          col.per.group = c("#f99943", "#44a7c4", "#05b052", "#2f6395", 
                            "#bb5352", "#87c242", "#07080a", "#92bbdb"),
          pch = depth.pch,
          pch.levels = TARAoceans$sample$depth,
          legend.title = "Ocean / Sea",
          title = "Projection of TARA Oceans stations",
          size.title = 10,
          legend.title.pch = "Depth")

## ----tune_pca-----------------------------------------------------------------
plot(kernel.pca.result)

## ----permute_kpca-------------------------------------------------------------
head(TARAoceans$taxonomy[ ,"Phylum"], 10)
head(TARAoceans$GO, 10)
# here we set a seed for reproducible results with this tutorial
set.seed(17051753)
kernel.pca.result <- kernel.pca.permute(kernel.pca.result, ncomp = 1,
                                        phychem = colnames(TARAoceans$phychem),
                                        pro.phylo = TARAoceans$taxonomy[, "Phylum"],
                                        pro.NOGs = TARAoceans$GO)

## ----display_var--------------------------------------------------------------
plotVar.kernel.pca(kernel.pca.result, ndisplay = 10, ncol = 3)

## ----proteobacteria_display, fig.keep='all'-----------------------------------
selected <- which(TARAoceans$taxonomy[, "Phylum"] == "Proteobacteria")
proteobacteria.per.sample <- apply(TARAoceans$pro.phylo[, selected], 1, sum) /
  apply(TARAoceans$pro.phylo, 1, sum)
colfunc <- colorRampPalette(c("royalblue", "red"))
col.proteo <- colfunc(length(proteobacteria.per.sample))
col.proteo <- col.proteo[rank(proteobacteria.per.sample, ties = "first")]
plotIndiv(kernel.pca.result,
          comp = c(1, 2),
          ind.names = FALSE,
          legend = FALSE,
          col = col.proteo,
          pch = depth.pch,
          pch.levels = TARAoceans$sample$depth,
          legend.title = "Ocean / Sea",
          title = "Representation of Proteobacteria abundance",
          legend.title.pch = "Depth")

## ----temperature_display, fig.keep='all'--------------------------------------
col.temp <- colfunc(length(TARAoceans$phychem[, 4]))
col.temp <- col.temp[rank(TARAoceans$phychem[, 4], ties = "first")]
plotIndiv(kernel.pca.result,
          comp = c(1, 2),
          ind.names = FALSE,
          legend = FALSE,
          col = col.temp,
          pch = depth.pch,
          pch.levels = TARAoceans$sample$depth,
          legend.title = "Ocean / Sea",
          title = "Representation of mean temperature",
          legend.title.pch = "Depth")

## ----dependencies-------------------------------------------------------------
have_depend <- reticulate::py_module_available("autograd") &
  reticulate::py_module_available("scipy") &
  reticulate::py_module_available("numpy") &
  reticulate::py_module_available("sklearn") 

## ----select_ukfs--------------------------------------------------------------
if (have_depend) {
  ukfs.res <- select.features(TARAoceans$pro.phylo, kx.func = "bray", lambda = 1, 
                              keepX = 5, nstep = 1)
  selected <- sort(ukfs.res, decreasing = TRUE, index.return = TRUE)$ix[1:5]
  TARAoceans$taxonomy[selected, ]
}

## ----comput_correlation_graph, eval=FALSE-------------------------------------
#  library("MASS")
#  library("igraph")
#  library("correlationtree")
#  
#  pro.phylo.alist <- data.frame("names" = colnames(TARAoceans$pro.phylo),
#                                t(TARAoceans$pro.phylo))
#  L <- mat2list(df2mat(pro.phylo.alist, 1))
#  corr.mat <- as.matrix(cross_cor(L, remove = TRUE))
#  pro.phylo.graph <- graph_from_adjacency_matrix(corr.mat,
#                                                 mode = "undirected",
#                                                 weighted = TRUE)
#  Lg <- laplacian_matrix(pro.phylo.graph, sparse=TRUE)

## ----select_ukfsg-------------------------------------------------------------
if (have_depend) {
  load(file = file.path(system.file(package = "mixKernel"), "loaddata", "Lg.rda"))
  ukfsg.res <- select.features(TARAoceans$pro.phylo, kx.func = "bray", 
                               lambda = 1, method = "graph", Lg = Lg, keepX = 5,
                               nstep = 1)
  
  selected <- sort(ukfsg.res, decreasing = TRUE, index.return = TRUE)$ix[1:5]
  TARAoceans$taxonomy[selected, ]
}

## ----session_information------------------------------------------------------
sessionInfo()

