\name{scores.to.counts}
\alias{scores.to.counts}
\title{transform a score matrix to a count matrix}
\description{
 transforms an object by classifier matrix of classification scores to
 an object by category matrix of classification counts as used by
 \code{cohen.kappa}.
}
\usage{
scores.to.counts(scores)
}
\arguments{
  \item{scores}{object by classifier matrix of classification scores}
}
\details{
 usually called by \code{cohen.kappa}
}
\value{
  \item{counts}{an object by category matrix of counts}
  \item{scores.to.counts}{the concordance/discordance matrix
   used in the calculation of alpha}
  \item{data.values}{a character vector of the original data values}
  \item{nmatchval}{the count of matches, used in calculation}
}
\author{Jim Lemon}
\note{There seems to be no established convention for whether objects(cases)
 should be rows (as they are here) or columns. Make sure to transpose as
 below if necessary.}
\seealso{\code{\link{cohen.kappa}}}
\examples{
 # the "C" data from Krippendorff
 nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
 # remember to transpose the data to the expected format
 scores.to.counts(t(nmm))
}
\keyword{misc}
