\name{tiecorr}
\alias{tiecorr}
\title{correction for tied ranks}
\description{
 calculates a correction for tied ranks
}
\usage{
 tiecorr(rankarray)
}
\arguments{
  \item{rankarray}{matrix of ranks}
}
\details{
 \code{tiecorr} calculates a correction for tied ranks that is appropriate for
 Kendall's W and the Spearman rank-order correlation coefficient.
}
\value{
  \item{tie3margsum}{correction factor}
}
\references{Siegel, S. & Castellan, N.J.Jr. (1988) Nonparametric
 statistics for the behavioral sciences. Boston, MA: McGraw-Hill.}
\author{Jim Lemon}
\seealso{\code{\link{kendall.w}}}
\keyword{misc}
