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
# as used by the *.kappa() functions

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

cohen.kappa<-function(classif,type=c("score","count")) {
 if(missing(classif))
  stop("Usage: cohen.kappa(classif, type=\"score\")\n")
 if(is.character(classif))
  classif<-apply(apply(classif,2,as.factor),2,as.numeric)
 if(type[1] == "score") classif.mat<-scores.to.counts(classif)
 else classif.mat<-as.matrix(classif)
 minclassif<-min(classif.mat) 
 # bump the minimum value up to 1 for tabulate
 #if(minclassif < 1) classif.mat<-classif.mat+1-minclassif
 classdim<-dim(classif)
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
 N<-matdim[1] # number of data objects
 ncat<-matdim[2] # number of categories
 if(type[1] == "score") {
  if(any(is.na(classif))) {
   cat("Can't use Cohen's method with NAs\n")
   PEc<-NA
  }
  else {
   pj<-apply(apply(classif,2,tabulate)/N,1,prod)
   PEc<-sum(pj)
  }
 }
 else PEc<-NA
 Cj<-apply(classif.mat,2,sum)
 pj<-Cj/(N*k[1])
 PEsc<-sum(pj^2)
 Si<-1/(k[1]*(k[1]-1))*sum(classif.mat*(classif.mat-1))
 PA<-(1/N)*sum(Si)
 Ksc<-(PA-PEsc)/(1-PEsc)
 if(type[1] == "score") {
  Kc<-(PA-PEc)/(1-PEc)
  varKc<-(2/(N*k[1]*(k[1]-1)))*
   (PEc-(2*k[1]-3)*PEc^2+2*(k[1]-2)*sum(pj^3))/(1-PEc)^2
  Zc<-Kc/sqrt(varKc)
  pc<-1-pnorm(Zc)
 }
 else {
  Kc<-NA
  Zc<-NA
 }
 varKsc<-(2/(N*k[1]*(k[1]-1)))*
  (PEsc-(2*k[1]-3)*PEsc^2+2*(k[1]-2)*sum(pj^3))/(1-PEsc)^2
 Zsc<-Ksc/sqrt(varKsc)
 psc<-1-pnorm(Zsc)
 Kbbc<-2*PA-1
 c.k<-list(kappa.c=Kc,kappa.sc=Ksc,kappa.bbc=Kbbc,
  Zc=Zc,Zsc=Zsc,pc=pc,psc=psc,categories=matdim[2],methods=k[1])
 class(c.k)<-"cohen.kappa"
 return(c.k)
}

print.cohen.kappa<-function(x,...) {
 cat("Kappa test for nominally classified data\n")
 cat(paste(x$categories,"categories -",x$methods,"methods\n"))
 if(!is.na(x$kappa.c)) {
  cat(paste("kappa (Cohen) =",signif(x$kappa.c),", Z =",signif(x$Zc),
   ", p =",signif(x$pc),"\n"))
 }
 cat(paste("kappa (Siegel) =",signif(x$kappa.sc),", Z =",signif(x$Zsc),
  ", p =",signif(x$psc),"\n"))
 cat(paste("kappa (2*PA-1) =",signif(x$kappa.bbc),"\n\n"))
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
