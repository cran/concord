\name{wtpc}
\alias{wtpc}
\title{calculate weighted percentages for nominal classification data}
\description{
 calculates weighted percentages for nominal classification data
}
\usage{
wtpc(x, n.methods, n.objects, type=c("count","score"))
}
\arguments{
  \item{x}{vector, data frame or matrix of ratings}
  \item{n.methods}{the number of methods (e.g. raters) used}
  \item{n.objects}{the number of objects rated or classified}
  \item{type}{whether the data represent counts of objects in categories
  or scores for objects by methods}
}
\details{
 \code{wtpc} will accept a vector (if only one method is used) matrix or 
 data frame in which values represent counts of objects (rows) classified 
 as categories (columns), or the same objects in which values represent 
 categorizations of those objects. The default is to assume counts 
 and the operator \bold{must} specify if scores are used. The formula
 used is: \code{wtpc = 100/n.methods * sum(category counts)/n.objects}.
}
\value{
  \item{wtpc}{vector of weighted percentages of occurrence of categories}
}
\author{Jim Lemon}
\note{This is a convenience function for calculating the frequencies of
category occurrences in ratings where all methods may not recognize the
same objects as belonging to the available categories.}
\examples{
 # the "C" data from Krippendorff
 nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
 # Krippendorff's data is classifier by object and must be transposed
 wtpc(t(nmm),4,12,"score")
}
\keyword{misc}
