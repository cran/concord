\name{kendall.w}
\alias{kendall.w}
\title{Kendall's W coefficient of concordance}
\description{
 calculates Kendall's W coefficient of concordance
}
\usage{
 kendall.w(x,lambda,descending=TRUE,ranks=FALSE)
}
\arguments{
  \item{x}{matrix of scores or ranks}
  \item{lambda}{optional contrast coefficient vector or matrix}
  \item{descending}{whether high (default) or low scores represent top ranks}
  \item{ranks}{whether the values in \code{x} are scores or ranks}
}
\details{
 \code{kendall.w} will accept either a matrix or data frame of scores or ranks
 where the rows represent methods (usually raters) and the columns represent
 data objects. It will handle ties, but not missing values. By default it 
 assumes that numerically higher scores represent numerically lower ranks. A
 vector or matrix of contrast coefficients (which each must sum to zero) may 
 be supplied. A z-score approximation of the significance of each contrast will
 be displayed if \code{lambda} is present. The \code{ranks} argument allows the
 user to pass ranks directly to the function. If ranks are passed without
 setting \code{ranks} to TRUE and with \code{descending} TRUE, the order of
 the ranks will be reversed. For small values of k (methods), \code{kendall.w} 
 will try to lookup the tabled values for significance. For \code{k} greater
 than 7, a chi-squared approximation is returned. Only one of these values will
 be returned.
}
\value{
 \item{W}{value of W}
 \item{p.table}{whether the obtained W exceeded the table value for small N}
 \item{p.chisq}{the probability of the obtained chi-squared value for larger N}
}
\references{Siegel, S. & Castellan, N.J.Jr. (1988) Nonparametric
 statistics for the behavioral sciences. Boston, MA: McGraw-Hill.}
\author{Jim Lemon}
\note{Kendall's W may not be appropriate for nominal class data.}
\seealso{\code{\link{kripp.alpha}},\code{\link{cohen.kappa}}}
\keyword{misc}
