\name{coincidence.matrix}
\alias{coincidence.matrix}
\title{calculate the concordance/discordance matrix for
 Krippendorff's alpha}
\description{
 calculates the coincidence matrix for Krippendorff's alpha.
}
\usage{
coincidence.matrix(x)
}
\arguments{
  \item{x}{classifier by object matrix of classifications}
}
\details{
 probably of no use by itself.
}
\value{
  \item{statistic}{always NA - value is calculated by \code{kripp.alpha}}
  \item{coincidence.matrix}{the concordance/discordance matrix
   used in the calculation of alpha}
  \item{data.values}{a character vector of the original data values}
  \item{nmatchval}{the count of matches, used in calculation}
}
\references{Krippendorff, K. (1980) Content analysis: An introduction to
 its methodology. Beverly Hills, CA: Sage.}
\author{Jim Lemon}
\note{This is the initial step in calculating Krippendorff's alpha coefficient.}
\seealso{\code{\link{kripp.alpha}}}
\examples{
 # the "C" data from Krippendorff
 nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
 coincidence.matrix(nmm)
}
\keyword{misc}
