\name{quest.reliability}
\alias{quest.reliability}
\title{Questionnaire reliability measures}
\description{Calculates a few common measures of questionnaire reliability}
\usage{
 quest.reliability(x,nsplits=10)
}
\arguments{
 \item{x}{a matrix of numeric scores with subjects as rows and items
  as columns}
 \item{nsplits}{the number of splits to make in calculating split-half
  reliability}
}
\value{
 A list containing five components:
 \item{cronbach.alpha}{Cronbach's alpha for the entire scale.}
 \item{split.half}{An approximate split-half reliability is
  calculated by dividing the items into halves on a random basis and 
  calculating the correlation between the summed scores for the halves.
  This is repeated \code{nsplits} times.}
 \item{item.whole}{The correlation of each item with the sum of the
  remaining items is calculated.}
 \item{item.alpha}{Cronbach's alpha is calculated for each item on all other
  items, i.e. the scale alpha if that item were dropped.}
 \item{item.var}{The variance of scores for each item is calculated.
  Items with relatively low variance are likely to be poor discriminators.}
}
\details{
 \code{quest.reliability} will only handle numeric scores. Character
 scores can be transformed into numeric using \code{as.numeric(as.factor(...}
 but take care that the transformation preserves the order of the scores.
 Any rows with NA values in the score matrix will be removed before the
 calculations are made.
}
\author{Jim Lemon}
\seealso{\code{\link{cor}},\code{\link{var}}}
\examples{
 # fake a data set with one bad item
 fakedata<-matrix(c(rep(c(2,2:9,9),9)+sample(c(-1,0,1),90,TRUE),
  rep(5,10)+sample(c(-1,0,1),10,TRUE)),nrow=10)
 quest.reliability(fakedata)
}
\keyword{misc}
