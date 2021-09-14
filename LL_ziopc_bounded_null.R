LL_ziopc_bounded_null <- function(coeff, Yoprob, Yprob = NULL, Ylevel = NULL, Ybounds = NULL, posc.low, posc.high, Nobs=1, onlyLL = TRUE) {

  # LL_ZIOPC_bounded computes the negative log likelihood of a zero inflated ordered probit
  # following Brooks/Harris/Spencer (2007). This function is not meant to be called manually
  # but is used as objective function for the optimization in ziop.r
  # It adjusts for boundaries as in Hantzsche 2015
  #
  # coeff     vector of coefficientss (probit eq, threshold probit, oprobit eq, thresholds oprob, <correlation>)
  # Yoprob    ordinal variable (has to be -1, 0, 1)
  # Yprob     binary variable. If missing, Yprob=|Yoprob|
  # Ylevel    Boundary adjustment, if necessary
  # X         matrix of exogenous variables for probit equation 
  # Z         matrix of exogenous variables for ordered probit equation 
  # posc.low  position of lower bound threshold in (-Inf,thresh1,thresh2,Inf) for every observation
  # posc.high position of upper bound threshold in (-Inf,thresh1,thresh2,Inf) for every observation
  # Nobs      Number of observations per agency
  # onlyLL (logical, default = TRUE) for use in optimization routines. LL_ziopc only returns log likelihood
  #
  # last changed: 12/18/2015
  #
  # author: mei, gsz 
  
  check = require(pbivnorm)
  if (!check) {
    install.packages("pbivnorm")
    require(pbivnorm)
  }
  
  if (is.null(Yprob)){Yprob<-abs(Yoprob)}
  
  lc = length(coeff)
  lN <- length(Nobs)
  if (lc==3*lN+1) {z = coeff[lc]} else {z = 0}
  rho <- (exp(2*z)-1)/(exp(2*z)+1)

  r_star = rep(coeff[1:lN],Nobs)
  con_ord <- c(-10e9,coeff[lN+1:(2*lN)],10e9)
  oprobmat <- cbind(con_ord[posc.low],con_ord[posc.high])
  
  # Pvec is the probability of the observed outcome
  Pvec.ann <- pbivnorm(cbind(r_star,oprobmat[,2]),rho=-rho)-pbivnorm(cbind(r_star,oprobmat[,1]),rho=-rho)
  Pvec.noann <- 1-pnorm(r_star)
  # Pvec.noann[Yprob==1] <- 0
  Pvec <- Pvec.noann
  Pvec[Yprob==1] <- Pvec.ann[Yprob==1]
  
  LL <- -sum(log(Pvec))
  if (onlyLL) {return(LL)}
  
  # the following assumes that constants are ordered into lower-upper threshold, separately by agency
  posc.low.down <- 1
  posc.high.down <- rep(seq(2,2*lN+1,by=2),Nobs) # sequence of first thresholds
  posc.low.up <- rep(seq(3,2*lN+1,by=2),Nobs)    # sequence of second thresholds
  posc.high.up <- 2*lN+2
  
  if (!is.null(Ylevel)){
    # boundary adjustment for boundary cases, if provided
    if (!is.null(Ybounds)){
      posc.high.down[(Ylevel==Ybounds[1])] <- 1
      posc.low.up[(Ylevel==Ybounds[2])] <- 2*lN+2
    }else{
      posc.high.down[(Ylevel==min(Ylevel))] <- 1
      posc.low.up[(Ylevel==max(Ylevel))] <- 2*lN+2
    }
  }
  
  oprobmat.down <- cbind(con_ord[posc.low.down]-y_star,con_ord[posc.high.down]-y_star)
  oprobmat.stay <- cbind(con_ord[posc.high.down]-y_star,con_ord[posc.low.up]-y_star)
  oprobmat.up <- cbind(con_ord[posc.low.up]-y_star,con_ord[posc.high.up]-y_star)
  
  p <- cbind(pbivnorm(cbind(r_star,oprobmat.down[,2]),rho=-rho)-pbivnorm(cbind(r_star,oprobmat.down[,1]),rho=-rho),
             1-pnorm(r_star),
             pbivnorm(cbind(r_star,oprobmat.stay[,2]),rho=-rho)-pbivnorm(cbind(r_star,oprobmat.stay[,1]),rho=-rho),
             pbivnorm(cbind(r_star,oprobmat.up[,2]),rho=-rho)-pbivnorm(cbind(r_star,oprobmat.up[,1]),rho=-rho),
             pnorm(oprobmat.down[,2])-pnorm(oprobmat.down[,1]),
             pnorm(oprobmat.stay[,2])-pnorm(oprobmat.stay[,1]),
             pnorm(oprobmat.up[,2])-pnorm(oprobmat.up[,1]),
             pnorm(r_star))
  colnames(p) = c("-1","00","01","1","-1.oprob","0.oprob","1.oprob","1.prob")
  
  res = list()
  res$LL = -LL
  res$p_case <- Pvec
  res$fitted.values = p
  res$latent.prob <- r_star
  res$latent.oprob <- oprobmat
  
  
  # res$fitted.values.oprob <- cbind(pnorm(-y_star),1-pnorm(-y_star)-pnorm(y_star-mu),pnorm(y_star-mu))
  # res$fitted.values.prob <- cbind(1-pchange,pchange)
  return(res)
}