\name{mcnemar.mh}
\alias{mcnemar.mh}
\title{2 way reliability coefficient for nominal data}
\description{
 calculates a coefficient of reliability for nominal data
}
\usage{
mcnemar.mh(x)
}
\arguments{
  \item{x}{2x2 classification matrix or matrix of dichotomous
   classification scores.}
}
\details{
 \code{mcnemar.mh} calculates a reliability coefficient for two raters
  classifying n objects into two categories. It will accept either a 2x2
  classification matrix of counts of objects falling into two categories
  or a 2xn or nx2 matrix of classification scores.
}
\value{
  \item{kappa}{value of kappa}
  \item{Z}{the Z-score approximation}
  \item{p}{the probability of Z}
}
\references{Siegel, S. & Castellan, N.J.Jr. (1988) Nonparametric
 statistics for the behavioral sciences. Boston, MA: McGraw-Hill.}
\author{Jim Lemon}
\seealso{\code{\link{stuart.maxwell.mh}}}
\examples{
 xt<-table(sample(0:1,50,TRUE),sample(0:1,50,TRUE))
 mcnemar.mh(xt)
}
\keyword{misc}
