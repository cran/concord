\name{stuart.maxwell.mh}
\alias{stuart.maxwell.mh}
\title{between rater concordance coefficient for classification}
\description{
 calculates the coefficient of concordance for two raters
}
\usage{
stuart.maxwell.mh(x)
}
\arguments{
 \item{x}{cxc classification matrix or matrix of classification scores
  into c categories.}
}
\details{
 \code{stuart.maxwell.mh} calculates a reliability coefficient for two raters
  classifying n objects into any number of categories. It will accept either 
  a cxc classification matrix of counts of objects falling into c categories
  or a cxn or nxc matrix of classification scores.
}
\value{
  \item{statistic}{value of coefficient}
  \item{p}{the probability of the coefficient as a Chi-square variable}
}
\references{
 Stuart A.A. (1955) A test for homogeneity of the marginal 
 distributions in a two-way classification. Biometrika, 42, 412-416.
 
 Maxwell A.E. (1970) Comparing the classification of subjects by two 
 independent judges. British Journal of Psychiatry, 116, 651-655.
}
\author{Jim Lemon}
\seealso{\code{\link{rater.bias}}}
\examples{
 # fake a 2xn matrix of three way classification scores
 ratings<-matrix(sample(1:3,60,TRUE),nrow=2)
 stuart.maxwell.mh(ratings)
}
\keyword{misc}
