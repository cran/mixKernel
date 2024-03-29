% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.kernels.R
\name{combine.kernels}
\alias{combine.kernels}
\title{Combine multiple kernels into a meta-kernel}
\usage{
combine.kernels(
  ...,
  scale = TRUE,
  method = c("full-UMKL", "STATIS-UMKL", "sparse-UMKL"),
  knn = 5,
  rho = 20
)
}
\arguments{
\item{...}{list of kernels (called 'blocks') computed on different datasets 
and measured on the same samples.}

\item{scale}{boleean. If \code{scale = TRUE}, each block is standardized to 
zero mean and unit variance and cosine normalization is performed on the 
kernel. Default: \code{TRUE}.}

\item{method}{character. Which method should be used to compute the 
meta-kernel. Default: \code{"full-UMKL"}.}

\item{knn}{integer. If \code{method = "sparse-UMKL"} or
\code{method = "full-UMKL"}, number of neighbors used to get a proxy of the 
local topology of the datasets from each kernel. Default: \code{5}.}

\item{rho}{integer. Parameters for the augmented Lagrangian method. Default: 
\code{20}.}
}
\value{
\code{combine.kernels} returns an object of classes \code{"kernel"} 
and \code{"metaKernel"}, a list that contains the following components: 
  \item{kernel}{: the computed meta-kernel matrix;}
  \item{X}{: the dataset from which the kernel has been computed, as given by
  the function \code{\link{compute.kernel}}. Can be \code{NULL} if a kernel
  matrix was passed to this function;}
  \item{weights}{: a vector containing the weights used to combine the 
  kernels.}
}
\description{
Compute multiple kernels into a single meta-kernel
}
\details{
The arguments \code{method} allows to specify the Unsupervised Multiple
Kernel Learning (UMKL) method to use: \itemize{
  \item \code{"STATIS-UMKL"}: combines input kernels into the best 
  consensus of all kernels;
  \item \code{"full-UMKL"}: computes a kernel that minimizes the distortion 
  between the meta-kernel and the k-NN graphs obtained from all input 
  kernels;
  \item \code{"sparse-UMKL"}: a sparse variant of the \code{"full-UMKL"} 
  approach.
}
}
\examples{
data(TARAoceans)

# compute one kernel per dataset
phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
pro.phylo.kernel <- compute.kernel(TARAoceans$pro.phylo, kernel.func = "abundance")
pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, kernel.func = "abundance")

# compute the meta kernel
meta.kernel <- combine.kernels(phychem = phychem.kernel,
                               pro.phylo = pro.phylo.kernel,
                               pro.NOGs = pro.NOGs.kernel, 
                               method = "full-UMKL")

}
\references{
Mariette J. and Villa-Vialaneix N. (2018). Unsupervised multiple 
kernel learning for heterogeneous data integration . \emph{Bioinformatics}, 
\bold{34}(6), 1009-1015. DOI: \doi{10.1093/bioinformatics/btx682}.
}
\seealso{
\code{\link{compute.kernel}}, \code{\link{kernel.pca}}
}
\author{
Jerome Mariette <jerome.mariette@inrae.fr>
Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>
}
