\name{kernel.pca.permute}
\alias{kernel.pca.permute}
\title{Assess variable importance}

\description{Assess importance of variables on a given PC component by computing 
the Crone-Crosby distance between original sample positions and sample positions 
obtain by a random permutation of the variables.}

\usage{
kernel.pca.permute(kpca.result, ncomp = 1, ..., directory = NULL)
}

\arguments{
\item{kpca.result}{a kernel.pca object returned by the \code{\link{kernel.pca}}
function.}
\item{ncomp}{integer. Number of KPCA components used to compute the importance. 
Default: \code{1}.}
\item{...}{list of character vectors. The parameter name must be the kernel name 
to be considered for permutation of variables. Provided vectors length has to be 
equal to the number of variables of the input dataset. A kernel is performed on 
each unique variables values. Crone-Crosby distances are computed on each KPCA 
performed on resulted kernels or meta-kernels and can be displayed using the 
\code{\link{plotVar.kernel.pca}}.}
\item{directory}{character. To limit computational burden, this argument allows
to store / read temporary computed kernels.}
}

\value{\code{kernel.pca.permute} returns a copy of the input \code{kpca.result}
results and add values in the three entries: \code{cc.distances}, 
\code{cc.variables} and \code{cc.blocks}.}

\section{References}{
Mariette J. and Villa-Vialaneix N. (2017). Integrating \emph{TARA} Oceans 
datasets using unsupervised multiple kernel learning. \emph{Preprint}

Crone L. and Crosby D. (1995). Statistical applications of a metric on subspaces
to satellite meteorology. \emph{Technometrics}, \bold{37}(3), 324-328.
}

\author{Jerome Mariette <jerome.mariette@inra.fr>

Nathalie Villa-Vialaneix <nathalie.villa-vialaneix@inra.fr>}

\seealso{\code{\link{compute.kernel}}, \code{\link{kernel.pca}}, 
\code{\link{plotVar.kernel.pca}}}

\examples{
data(TARAoceans)

# compute one kernel for the psychem dataset
phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
# perform a KPCA
kernel.pca.result <- kernel.pca(phychem.kernel)

# compute importance for all variables in this kernel
kernel.pca.result <- kernel.pca.permute(kernel.pca.result, phychem = colnames(TARAoceans$phychem))
}
