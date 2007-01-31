quest.reliability<-function(x,nsplits=10) {
 if(any(is.na(x))) x<-x[complete.cases(x),]
 nitems<-dim(x)[2]
 oa<-(nitems/(nitems-1))*(1-sum(apply(x,2,var))/var(apply(x,1,sum)))
 ivar<-iir<-iwcors<-rep(0,nitems)
 shcors<-rep(0,nsplits)
 for(i in 1:nsplits) {
  allitems<-1:nitems
  hi1<-sample(allitems,nitems/2)
  hi2<-allitems[!(allitems%in%hi1)]
  shcors[i]<-cor(apply(x[,hi1],1,sum),apply(x[,hi2],1,sum))
 }
 for(i in 1:nitems) {
  iwcors[i]<-cor(x[,i],apply(x[,-i],1,sum))
  iir[i]<-
   ((nitems-1)/(nitems-2))*(1-sum(apply(x[,-i],2,var))/var(apply(x[,-i],1,sum)))
  ivar[i]<-var(x[,i])
 }
 sr.list<-list(cronbach.alpha=oa,split.half=shcors,item.whole=iwcors,
  item.alpha=iir,item.var=ivar)
 class(sr.list)<-"scale.reliability"
 return(sr.list)
}

print.quest.reliability<-function(x,...) {
 cat("Questionnaire reliability analysis\n\n")
 if(!is.null(x$cronbach.alpha))
  cat("Cronbach's alpha =",format(x$cronbach.alpha,digits=3),"\n")
 if(!is.null(x$split.half))
  cat("Split-half reliability (",length(x$split.half),
   " samples) = ",format(mean(x$split.half),digits=2),"\n",sep="")
 if(!is.null(x$item.whole)) {
  nitems<-length(x$item.whole)
  cat("Item Item/Total Alpha-Item Item var.\n")
  for(i in 1:nitems) {
   cat(format(as.character(i),width=8,justify="left"),
    format(as.character(format(x$item.whole[i],digits=2)),width=10,justify="left"),
    format(as.character(format(x$item.alpha[i],digits=2)),width=10,justify="left"),
    format(as.character(format(x$item.var[i],digits=2)),justify="left"),"\n")
  }
  cat("\n")
 }
}
