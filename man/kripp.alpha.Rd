\name{kripp.alpha}
\alias{kripp.alpha}
\title{calculate Krippendorff's alpha reliability coefficient}
\description{
 calculates the alpha coefficient of reliability proposed by Krippendorff.
}
\usage{
kripp.alpha(x, method="nominal"))
}
\arguments{
  \item{x}{classifier by object matrix of classifications}
  \item{method}{data level of x}
}
\details{
}
\value{
  \item{statistic}{value of alpha}
  \item{coincidence.matrix}{the concordance/discordance matrix
   used in the calculation of alpha}
  \item{data.values}{a character vector of the unique data values}
  \item{nmatchval}{the count of matches, used in calculation}
}
\references{Krippendorff, K. (1980) Content analysis: An introduction to
 its methodology. Beverly Hills, CA: Sage.}
\author{Jim Lemon}
\note{Krippendorff's alpha coefficient is particularly useful where
 the level of measurement of classification data is higher than nominal
 or ordinal.}
\seealso{\code{\link{coincidence.matrix}}}
\examples{
 # the "C" data from Krippendorff
 nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
 # first assume the default nominal classification
 kripp.alpha(nmm)
 # now use the same data with the other three methods
 kripp.alpha(nmm,"ordinal")
 kripp.alpha(nmm,"interval")
 kripp.alpha(nmm,"ratio") 
}
\keyword{misc}
