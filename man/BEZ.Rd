\name{BEZ}
\alias{BEZ}
\title{Z score approximation of contrast significance}
\description{
 calculates a Z score approximation of contrast significance
}
\usage{
 BEZ(rankarray,lambda)
}
\arguments{
  \item{rankarray}{matrix of ranks}
  \item{lambda}{contrast coefficient vector}
}
\details{
 \code{BEZ} calculates a Z score approximation of the probability that the
 deviation of the mean rank scores from the contrast coefficient vector would
 occur given randomly allocated ranks. It is only useful in conjuction with
 the /code{kendall.w} function.
}
\value{
  \item{Z}{Z score}
}
\author{Jim Lemon}
\seealso{\code{\link{kendall.w}}}
\keyword{misc}
