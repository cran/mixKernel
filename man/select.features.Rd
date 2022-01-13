\name{select.features}
\alias{select.features}
\title{Select important features}

\description{Select features using supervised or unsupervised kernel method. A 
supervised feature selection method is performed if \code{Y} is provided.}

\usage{
select.features(X, Y = NULL, 
                kx.func = c("linear", "gaussian.radial.basis", "bray"),
                ky.func = c("linear", "gaussian.radial.basis"), keepX = NULL,
                method = c("kernel", "kpca", "graph"), lambda = NULL,
                n_components = 2, Lg = NULL, mu = 1, max_iter = 100, nstep = 50,
                ...)
}

\arguments{
\item{X}{a numeric matrix (or data frame) used to select variables.
\code{NA}s not allowed.}
\item{Y}{a numeric matrix (or data frame) used to select variables.
\code{NA}s not allowed.}
\item{kx.func}{the kernel function name to use on \code{X}. Widely used kernel
functions are pre-implemented, and can be directly used by setting 
\code{kx.func} to one of the following values: \code{"linear"}, 
\code{"gaussian.radial.basis"} or \code{"bray"}. Default: \code{"linear"}. If 
\code{Y} is provided, the kernel \code{"bray"} is not allowed.}
\item{ky.func}{the kernel function name to use on \code{Y}. Available kernels 
are: \code{"linear"}, and \code{"gaussian.radial.basis"}. Default: 
\code{"linear"}. This value is ignored when \code{Y} is not provided.}
\item{keepX}{the number of variables to select.}
\item{method}{the method to use. Either an unsupervised variable selection
method (\code{"kernel"}), a kernel PCA oriented variable selection method 
(\code{"kpca"}) or a structure driven variable selection selection 
(\code{"graph"}). Default: \code{"kernel"}.}
\item{lambda}{the penalization parameter that controls the trade-off between the
minimization of the distorsion and the sparsity of the solution parameter.}
\item{n_components}{how many principal components should be used with method
\code{"kpca"}. Required with method \code{"kpca"}. Default: \code{2}.}
\item{Lg}{the Laplacian matrix of the graph representing relations between the 
input dataset variables. Required with method \code{"graph"}.}
\item{mu}{the penalization parameter that controls the trade-off between the
the distorsion and the influence of the graph. Default: \code{1}.}
\item{max_iter}{the maximum number of iterations. Default: \code{100}.}
\item{nstep}{the number of values used for the regularization path. Default: 
\code{50}.}
\item{...}{the kernel function arguments. In particular \code{sigma}
(\code{"gaussian.radial.basis"}): double. The inverse kernel width used by
\code{"gaussian.radial.basis"}.}
}

\value{
\code{ukfs} returns a vector of sorted selected features indexes.
}

\references{
Brouard C., Mariette J., Flamary R. and Vialaneix N. (2022). Feature selection
for kernel methods in systems biology. \emph{NAR Genomics and Bioinformatics},
\emph{Forthcoming}.
}

\author{Celine Brouard <celine.brouard@inrae.fr>

Jerome Mariette <jerome.mariette@inrae.fr>

Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>}

\seealso{\code{\link{compute.kernel}}}

\examples{
## These examples require the installation of python modules
## See installation instruction at: http://mixkernel.clementine.wf

data("Koren.16S")
\dontrun{
sf.res <- select.features(Koren.16S$data.raw, kx.func = "bray", lambda = 1,
                          keepX = 40, nstep = 1)
colnames(Koren.16S$data.raw)[sf.res]
}

data("nutrimouse")
\dontrun{
grb.func <- "gaussian.radial.basis"
genes <- center.scale(nutrimouse$gene)
lipids <- center.scale(nutrimouse$lipid)
sf.res <- select.features(genes, lipids, kx.func = grb.func, ky.func = grb.func,
                          keepX = 40)
colnames(nutrimouse$gene)[sf.res]
}
}
