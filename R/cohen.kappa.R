# accepts a "score matrix" like:
#
# 		rater1	rater2	...
# object1	1	2
# object2	1	1
# ...
#
# and returns a "count matrix" like:
#
#		1	2	...
# object1	1	1
# object2	2	0
# ...
#
# as used by kappa.nom()

scores.to.counts<-function(scores) {
 if(is.data.frame(scores))
  score.names<-levels(as.factor(as.vector(unlist(scores))))
 if(is.matrix(scores))
  score.names<-levels(as.factor(as.vector(scores)))
 if(missing(score.names)) stop("scores must be a data frame or matrix")
 score.levels<-as.numeric(score.names)
 nlevels<-length(score.levels)
 nobj<-length(scores[,1])
 counts<-matrix(0,nrow=nobj,ncol=nlevels)
 colnames(counts)<-score.names
 for(i in 1:nobj) {
  for(j in 1:nlevels) counts[i,j]<-sum(scores[i,]==score.levels[j],na.rm=TRUE)
 }
 return(counts)
}

cohen.kappa<-function(classif,type=c("count","score")) {
 if(!missing(classif)) {
  if(type == "score") classif.mat<-scores.to.counts(classif)
  else classif.mat<-as.matrix(classif)
  k<-apply(classif.mat,1,sum)
  # check that all the row sums are equal
  if(any(k != k[1])) {
   # stick on an extra column of no-classification counts
   classif.mat<-cbind(classif.mat,max(k)-k)
   # recalculate the row sums
   k<-apply(classif.mat,1,sum)
   # let the user know
   cat("Different row sums, a no-classification category was added.\n\n")
  }
  matdim<-dim(classif.mat)
  N<-matdim[1]
  Cj<-apply(classif.mat,2,sum)
  pj<-Cj/(N*k[1])
  PE<-sum(pj^2)
  Si<-1/(k[1]*(k[1]-1))*sum(classif.mat*(classif.mat-1))
  PA<-(1/N)*sum(Si)
  K<-(PA-PE)/(1-PE)
  varK<-(2/(N*k[1]*(k[1]-1)))*
   (PE-(2*k[1]-3)*PE^2+2*(k[1]-2)*sum(pj^3))/(1-PE)^2
  Z<-K/sqrt(varK)
  p<-1-pnorm(Z)
  cat("Kappa test for nominally classified data\n")
  cat(matdim[2],"categories -",k[1],"methods\n\n")
  cat("kappa =",K,", z =",Z,", p =",p,"\n")
  invisible(list(kappa=K,Z=Z,p=p))
 }
 else {
  cat("Usage: cohen.kappa(classif,type=c(\"count\",\"score\"))\n")
  cat("\twhere classif is a data frame or matrix of counts\n")
  cat("\tof data objects (as rows) into categories (as columns)\n")
  cat("\tor of data objects (as rows) into category scores of\n")
  cat("\traters (as columns). If the latter, specify type=\"score\"\n")  
  cat("\tNote: if all row sums (number of classifiers) are not equal,\n")
  cat("\ta no-classification category will be added to make them so.\n")
 }
}

# wtpc calculates the weighted percentages using the following formula:
# <weighted pc><-(100/<n methods>)*<n ratings>/<n data objects>
# The format of the data is the same as that used for calculating the
# kappa for nominal data kappa.nom()

wtpc<-function(x,n.methods,n.objects,type=c("count","score")) {
 if(!missing(x) && !missing(n.methods) && !missing(n.objects)) {
  if(type == "score") x<-scores.to.counts(x)
  if(is.data.frame(x)) sumx<-sapply(x,sum)
  if(is.matrix(x)) sumx<-apply(x,2,sum)
  else sumx<-sum(x)
  return((100/n.methods)*sumx/n.objects)
 }
 else {
  cat("Usage: wtpc(x,n.methods,n.objects,type=c(\"count\",\"score\"))\n")
  cat("\twhere x is a vector, data frame or matrix of ratings,\n")
  cat("\tif x is scores rather than counts, specify type=score\n")
  cat("\tn.methods is the number of rating methods (e.g. raters)\n")
  cat("\tand n.objects is the number of data objects (e.g. subjects)\n")
 }
}
