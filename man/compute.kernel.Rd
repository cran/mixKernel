\name{compute.kernel}
\alias{compute.kernel}
\title{Compute a kernel}

\description{Compute a kernel from a given data matrix.}

\usage{
compute.kernel(X, kernel.func = "linear", ..., test.pos.semidef = FALSE)
}

\arguments{
\item{X}{a numeric matrix (or data frame) used to compute the kernel. 
\code{NA}s not allowed.}
\item{kernel.func}{the kernel function to use. This parameter can be set to any
user defined kernel function. Widely used kernel functions are pre-implemented,
that can be used by setting \code{kernel.func} to one of the following strings:
\code{"kidentity"}, \code{"abundance"}, \code{"linear"} or 
\code{"phylogenetic"}. Default: \code{"linear"}.}
\item{...}{the kernel function arguments. Valid parameters for pre-implemented 
kernels are: 
\itemize{
  \item \code{phylogenetic.tree}(\code{"phylogenetic"}) an instance of 
  phylo-class that contains a phylogenetic tree (required).
  \item \code{scale} (\code{"linear"}): logical. Should the variables be 
  scaled to unit variance prior the kernel computation? Default: \code{TRUE}.
  \item \code{method}(\code{"phylogenetic"} or 
  computing the kernel? Can be \code{"unifrac"} or \code{"wunifrac"} for
  \code{"abundance"}): character. Which dissimilarity to use for 
  \code{"phylogenetic"} and one of \code{"bray"},
  \code{"euclidean"}, \code{"canberra"}, \code{"manhattan"}, \code{"kulczynski"}, 
  \code{"jaccard"}, \code{"gower"}, \code{"altGower"}, \code{"morisita"},
  \code{"horn"}, \code{"mountford"}, \code{"raup"}, \code{"binomial"}, 
  \code{"chao"} and \code{"cao"} for \code{"abundance"}.}
}
\item{test.pos.semidef}{boleean. If \code{test.pos.semidef = TRUE}, the resulting
matrix is tested to be positive-semidefinite.}
}

\value{
\code{compute.kernel} returns an object of classes \code{"kernel"}, a list that 
contains the following components:
  \item{kernel}{: the computed kernel matrix.}
  \item{X}{: the original dataset. If \code{"kidentity"}, \code{X} is set to \code{NULL}.}
  \item{kernel.func}{: the kernel function used.}
  \item{kernel.args}{: the arguments used to compute the kernel.} 
}

\references{
Lozupone C. and Knight R. (2005). UniFrac: a new phylogenetic method for 
comparing microbial communities. \emph{Applied and Environmental Microbiology},
\bold{71}(12), 8228-8235.

Lozupone C., Hamady M., Kelley S.T. and Knight R. (2007). Quantitative and 
qualitative beta diversity measures lead to different insights into factors that
structure microbial communities. \emph{Applied and Environmental Microbiology}, 
\bold{73}(5), 1576-1585.
}

\author{Jerome Mariette <jerome.mariette@inra.fr>

Nathalie Villa-Vialaneix <nathalie.villa-vialaneix@inra.fr>}

\seealso{
\code{\link{combine.kernels}}, \code{\link{kernel.pca}}
}

\examples{
data(TARAoceans)
pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, kernel.func = "abundance")
}
