\name{center.scale}
\alias{center.scale}
\title{Center and scale}

\description{Center and scale a dataset.}

\usage{
center.scale(X)
}

\arguments{
\item{X}{a numeric matrix (or data frame) to center and scaled.
\code{NA}s not allowed.}
}

\value{\code{center.scale} returns a centered and scaled matrix.
}

\author{Celine Brouard <celine.brouard@inrae.fr>
     
Jerome Mariette <jerome.mariette@inrae.fr>

Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>}

\seealso{
\code{\link{compute.kernel}}, \code{\link{combine.kernels}}
}

\examples{
data("nutrimouse")
\dontrun{
nutrimouse.sc <- center.scale(nutrimouse$gene)
}
}
