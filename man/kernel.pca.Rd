\name{kernel.pca}
\alias{kernel.pca}
\title{Kernel Principal Components Analysis}

\description{Performs a kernel PCA.}

\usage{
kernel.pca(K, ncomp = nrow(K$kernel))
}

\arguments{
\item{K}{a kernel object obtained using either \code{compute.kernel} or
\code{combine.kernels}.}
\item{ncomp}{integer. Indicates the number of components to return.}
}

\value{\code{kernel.pca} returns an object of classes \code{"kernel.pca"} and 
\code{"pca"}, which is a list containing the following entries: \itemize{
  \item{ncomp}{: the number of principal components;}
  \item{X}{: the input kernel matrix;} 
  \item{kernel}{: the input kernel object provided by the user;}
  \item{sdev}{: the singular values (square root of the eigenvalues);} 
  \item{rotation}{: the matrix of variable loadings (\emph{i.e.}, a matrix whose
  columns contain the eigenvectors);}
  \item{loadings}{: same as 'rotation' to keep the mixOmics spirit;}
  \item{x}{: same as 'rotation' to keep the mixOmics spirit;}
}}

\section{References}{
Scholkopf B., Smola A. and Muller K.R. (1998) Nonlinear component analysis as a
kernel eigenvalue problem. \emph{Neural Computation}, \bold{10}, 1299-1319.
}

\author{Jerome Mariette <jerome.mariette@inra.fr>

Nathalie Villa-Vialaneix <nathalie.villa-vialaneix@inra.fr>}

\seealso{
\code{\link{compute.kernel}}, \code{\link{combine.kernels}}
}

\examples{
data(TARAoceans)
phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
kernel.pca.result <- kernel.pca(phychem.kernel, ncomp = 3)
}
