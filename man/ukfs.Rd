\name{ukfs}
\alias{ukfs}
\title{UKFS}

\description{Select variables using unsupervised kernel method.}

\usage{
ukfs(X, kernel.func=c("linear", "gaussian.radial.basis", "bray"),
     method=c("kernel", "kpca", "graph"),
     keepX=NULL, lambda=NULL, n_components=2, Lg=NULL,
     mu=1, max_iter=100, ...)
}

\arguments{
\item{X}{a numeric matrix (or data frame) used to select variables.
\code{NA}s not allowed.}
\item{kernel.func}{the kernel function name to use. Widely used kernel functions
are pre-implemented, that can be used by setting \code{kernel.func} to one of the
following strings: \code{"linear"}, \code{"gaussian.radial.basis"} or \code{"bray"}. 
Default: \code{"linear"}.}
\item{method}{the method to use. Either an unsupervised variable selection
method (\code{"kernel"}), a kernel PCA oriented variable selection method (\code{"kpca"})
or a structure driven variable selection selection (\code{"graph"}).
Default: \code{"kernel"}.}
\item{keepX}{the number of variables to select.}
\item{lambda}{the penalization parameter that controls the trade-off between the
minimization of the distorsion and the sparsity of the solution parameter.}
\item{n_components}{how many principal components should be used with method
\code{"kpca"}.
Required with method \code{"kpca"}.
Default: \code{2}.}
\item{Lg}{the Laplacian matrix of the graph representing relations between the input
dataset variables.
Required with method \code{"graph"}.}
\item{mu}{the penalization parameter that controls the trade-off between the
the distorsion and the influence of the graph.
Default: \code{1}.}
\item{max_iter}{the maximum number of iterations.
Default: \code{100}.}
\item{...}{the kernel function arguments. In particular \code{sigma}
(\code{"gaussian.radial.basis"}): double. The inverse kernel width used by
\code{"gaussian.radial.basis"}.}
}

\value{
\code{ukfs} returns a vector of size.
}

\references{
Mariette J., Brouard C. Flamary R. and Vialaneix N. (2020). Unsupervised variable selection
for kernel methods in systems biology. \emph{Preprint} 
}

\author{Jerome Mariette <jerome.mariette@inra.fr>

Nathalie Villa-Vialaneix <nathalie.villa-vialaneix@inra.fr>}

\seealso{\code{\link{compute.kernel}}}

\examples{
data("Koren.16S")
\dontrun{
ukfs.res <- ukfs(Koren.16S$data.raw, kernel.func = "bray", lambda=1)
}
}
