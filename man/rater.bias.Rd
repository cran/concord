\name{rater.bias}
\alias{rater.bias}
\title{coefficient of rater bias}
\description{
 calculates a coefficient of systematic bias between two raters
}
\usage{
rater.bias(x)
}
\arguments{
 \item{x}{CxC classification matrix or matrix of classification scores
  into C categories.}
}
\details{
 \code{rater.bias} calculates a reliability coefficient for two raters
  classifying n objects into any number of categories. It will accept either 
  a cxc classification matrix of counts of objects falling into c categories
  or a 2xn or nx2 matrix of classification scores.
}
\value{
  \item{statistic}{value of coefficient}
  \item{p}{the probability of the coefficient as a df 1 Chi-square variable}
}
\references{
 Bishop Y.M.M., Fienberg S.E. & Holland P.W. (1975) Discrete multivariate 
 analysis: theory and practice. Cambridge, Massachusetts: MIT Press.
}
\author{Jim Lemon}
\seealso{\code{\link{mcnemar.mh}}}
\examples{
 # fake a 2xn matrix of three way classification scores
 ratings<-matrix(sample(1:3,60,TRUE),nrow=2)
 rater.bias(ratings)
}
\keyword{misc}
