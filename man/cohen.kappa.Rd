\name{cohen.kappa}
\alias{cohen.kappa}
\title{kappa reliability coefficient for nominal data}
\description{
 calculates the kappa coefficient of reliability for nominal data
}
\usage{
cohen.kappa(classif, type=c("count","score"))
}
\arguments{
  \item{classif}{matrix of classification counts or scores}
  \item{type}{whether classif is an object by category matrix of counts
   or an object by classifier matrix of scores}
}
\details{
 \code{cohen.kappa} will accept either an object by category matrix of counts
 in which the numbers represent how many classifiers have placed the object
 in each category, or an object by classifier matrix of categories in which
 the numbers represent each classifier's categorization of that object. The
 default is to assume counts and the operator \bold{must} specify if scores
 are used.
}
\value{
  \item{kappa}{value of kappa}
  \item{Z}{the Z-score approximation}
  \item{p}{the probability of Z}
}
\references{Siegel, S. & Castellan, N.J.Jr. (1988) Nonparametric
 statistics for the behavioral sciences. Boston, MA: McGraw-Hill.}
\author{Jim Lemon}
\note{This is sometimes called Cohen's kappa. The name also avoids confusion
 with the kappa estimate of the conditioning number of a matrix. For a 
 contingency table version of this statistic, see \code{classAgreement}
 in package e1071}
\seealso{\code{\link{scores.to.counts}},\code{\link{wtpc}}}
\examples{
 # the "C" data from Krippendorff
 nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
 # first show the score to count transformation, remembering that
 # Krippendorff's data is classifier by object and must be transposed
 scores.to.counts(t(nmm))
 # now calculate kappa
 cohen.kappa(t(nmm),"score")
}
\keyword{misc}
