LL_oprob_bounded.probs <- function(probs, Yoprob, Ylevel) {
  
  pdown <- probs[,"-1"]
  pup <- probs[,"1"]
  pstay <- probs[,"0"]
  
  boundlow <- Ylevel==min(Ylevel)
  boundhigh <- Ylevel==max(Ylevel)
  
  pstay[boundlow] <- pstay[boundlow]+pdown[boundlow]
  pdown[boundlow==1] <- 0
  pstay[boundhigh] <- pstay[boundhigh]+pup[boundhigh]
  pup[boundhigh] <- 0
  LL <- sum(log(pdown[Yoprob==-1])) + 
    sum(log(pstay[Yoprob==0])) + 
    sum(log(pup[Yoprob==1]))
  
  if (is.na(LL)){
    print("NaN likelihood")
    print(summary(pdown[Yoprob==-1]))
    print(summary(pup[Yoprob==1]))
    print(summary(pstay[Yoprob==0]))
  }
  return(LL)
}