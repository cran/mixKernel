\name{cim.kernel}
\alias{cim.kernel}
\title{Compute and display similarities between multiple kernels}

\description{Compute cosine from Frobenius norm between kernels and display the 
corresponding correlation plot.}

\usage{
cim.kernel(..., scale = TRUE, 
           method = c("circle", "square", "number", "shade", "color", "pie"))
}

\arguments{
\item{...}{list of kernels (called 'blocks') computed on different datasets and
measured on the same samples.}
\item{scale}{boleean. If \code{scale = TRUE}, each block is standardized to zero
mean and unit variance and cosine normalization is performed on the kernel.
Default: \code{TRUE}.}
\item{method}{character. The visualization method to be used. Currently, seven
methods are supported (see Details).}
}

\details{The displayed similarities are the kernel generalization of the 
RV-coefficient described in Lavit \emph{et al.}, 1994.

The plot is displayed using the \code{\link[corrplot]{corrplot}} package. Seven 
visualization methods are implemented: \code{"circle"} (default), 
\code{"square"}, \code{"number"}, \code{"pie"}, \code{"shade"} and 
\code{"color"}. Circle and square areas are proportional to the absolute value
of corresponding similarities coefficients.}

\value{
\code{cim.kernel} returns a matrix containing the cosine from Frobenius norm 
between kernels.
}

\references{
Lavit C., Escoufier Y., Sabatier R. and Traissac P. (1994). The ACT (STATIS
method). \emph{Computational Statistics and Data Analysis}, \bold{18}(1), 
97-119.

Mariette J. and Villa-Vialaneix N. (2017). Integrating \emph{TARA} Oceans 
datasets using unsupervised multiple kernel learning. \emph{Preprint}
}

\author{Jerome Mariette <jerome.mariette@inrae.fr>

Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>}

\seealso{\code{\link{compute.kernel}}}

\examples{
data(TARAoceans)

# compute one kernel per dataset
phychem.kernel <- compute.kernel(TARAoceans$phychem, kernel.func = "linear")
pro.phylo.kernel <- compute.kernel(TARAoceans$pro.phylo, kernel.func = "abundance")
pro.NOGs.kernel <- compute.kernel(TARAoceans$pro.NOGs, kernel.func = "abundance")

# display similarities between kernels
cim.kernel(phychem = phychem.kernel,
           pro.phylo = pro.phylo.kernel,
           pro.NOGs = pro.NOGs.kernel, 
           method = "square")
}
