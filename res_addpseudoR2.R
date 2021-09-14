res_addpseudoR2 <- function(res,Yoprob,Yprob,Ylevel,posc.low, posc.high,Nobs,ziop=FALSE){
  
  fv <- res$fitted.values
  posannounce <- Yprob==1
  N <- length(Yprob)
  
  res$pred.Yprob <- (fv[,"1.prob"]>0.5)*1
  res$pred.Yoprob <- (fv[,"1.oprob"]>0.5)*1 - (res$fitted.values[,"-1.oprob"]>0.5)*1
  res$corrpred.prob <- mean(res$pred.Yprob == Yprob)
  res$corrpred.oprob <- mean(res$pred.Yoprob[posannounce] == Yoprob[posannounce])
  res$corrpred.full <- (sum(fv[Yoprob==-1,"-1"]>0.5) + sum(fv[Yprob==0,"00"]>0.5) + 
                          sum(fv[Yoprob==0 & Yprob==1,"01"]>0.5) + sum(fv[Yoprob==1,"1"]>0.5)) / N
  
  
  # polr_null <- polr(formula = "Yoprob~Inter",data=data.frame(Yoprob=as.factor(Yoprob[posannounce]),Inter=1),method="probit")
  # thres_polr <- polr_null$zeta
  # coeffs_null <- c(qnorm(mean(Yprob)),-thres_polr[1],log(thres_polr[2]-thres_polr[1]))
  coeffs_null <- c(res$prob_coeff[(length(res$prob_coeff)-length(Nobs)+1):length(res$prob_coeff)],
                   res$oprob_coeff[(length(res$oprob_coeff)-2*length(Nobs)+1):length(res$oprob_coeff)])
  print(coeffs_null)
  if (ziop){
    resnull <- optimx(coeffs_null,LL_ziopc_bounded_null,Yoprob=Yoprob,Yprob=NULL,Ylevel=Ylevel,posc.low=posc.low, posc.high=posc.high,Nobs=Nobs,method="BFGS")
  }else{
    resnull <- optimx(coeffs_null,LL_ziopc_bounded_null,Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,posc.low=posc.low, posc.high=posc.high,Nobs=Nobs,method="BFGS")
  }
  print(resnull)
  
  res$LLnull <- -resnull[,"value"]
  res$pseudoR2 <- 1-res$LL/res$LLnull
  if (!ziop){
    res$LLnull.oprob <- -polr(formula = "Yoprob~Inter",data=data.frame(Yoprob=as.factor(Yoprob),Inter=1),method="probit")$deviance
    res$pseudoR2.oprob <- 1-res$LL.oprob/res$LLnull.oprob
    res$pseudoR2.oprob.inflated <- 1-res$LL.oprob.inflated/res$LLnull.oprob
  }
  
  return(res)
}

