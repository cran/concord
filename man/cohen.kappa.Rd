\name{cohen.kappa}
\alias{cohen.kappa}
\title{kappa reliability coefficient for nominal data}
\description{
 calculates the kappa coefficient of reliability for nominal data
}
\usage{
cohen.kappa(classif, type=c("score","count"))
}
\arguments{
  \item{classif}{matrix of classification counts or scores}
  \item{type}{whether classif is an object by method matrix of scores
   or an object by category matrix of counts}
}
\details{
 \code{cohen.kappa} will accept either an object by category matrix of counts
 in which the numbers represent how many methods have placed the object
 in each category, or an object by method matrix of categories in which
 the numbers represent each method's categorization of that object. The
 default is to assume scores and the operator \bold{must} specify if counts
 are used.
 \code{cohen.kappa} reports two or three kappa values. If the classification
 matrix is composed of scores, the first is the original calculation from 
 Cohen(1960) which does not assume equal classification proportions for the 
 different methods. The next value is calculated as in Siegel & Castellan 
 (1988) and uses pooled classification proportions. This method provides an 
 adjustment for bias, where the different methods systematically differ in 
 their categorization. The third value is adjusted for prevalence using the 
 method proposed by Byrt, Bishop and Carlin (1993). An approximate 
 distribution of this statistic does not seem to be available, so there is 
 no approximation or probability reported.
}
\value{
  \item{kappa.c}{value of kappa (Cohen)}
  \item{kappa.sc}{value of kappa (Siegel & Castellan)}
  \item{kappa.bbc}{value of kappa (Byrt, Bishop & Carlin)}
  \item{Zc}{the Z-score approximation for kappa.c}
  \item{Zsc}{the Z-score approximation for kappa.sc}
  \item{pc}{the probability for Zc}
  \item{psc}{the probability for Zsc}
}
\references{
 {Byrt, T., Bishop, J. & Carlin, J.B. (1993) Bias, Prevalence and Kappa.
  Journal of Clinical Epidemiology, 46(5): 423-429.}
 {Cohen, J. (1960) A coefficient of agreement for nominal scales.
  Educational and Psychological Measurement, 20: 37-46.}
 {Siegel, S. & Castellan, N.J.Jr. (1988) Nonparametric
 statistics for the behavioral sciences. Boston, MA: McGraw-Hill.}
}
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
 # now calculate kappa - note that Cohen's method does not work with NAs
 cohen.kappa(t(nmm),"score")
}
\keyword{misc}
