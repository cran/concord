# mcnemar.mh computes the simple 2x2 McNemar test for marginal
# homogeneity.

mcnemar.mh<-function(x) {
 if(is.matrix(x) || is.data.frame(x)) {
  dimx<-dim(x)
  if(length(dimx) == 2) {
   if(any(dimx>2)) {
    if(dimx[1] == 2) mnx<-as.matrix(table(x[1,],x[2,]))
    else mnx<-as.matrix(table(x[,1],x[,2]))
   }
   else mnx<-as.matrix(x)
   mns<-(mnx[1,2]-mnx[2,1])^2/(mnx[1,2]+mnx[2,1])
   if((mnx[1,2]+mnx[2,1]) < 10)
    warning("low cell counts - consider binomial test")
   return(list(statistic=mns,p=1-pchisq(mns,1)))
  }
  else cat("Dimension higher than 2, cannot compute\n")
 }
 cat("Usage: mcnemar.mh(x)\n")
 cat("\twhere x is an nx2 matrix or data frame of scores\n")
 cat("\tor a 2x2 matrix or data frame of rater agreement on 2 categories\n")
}
