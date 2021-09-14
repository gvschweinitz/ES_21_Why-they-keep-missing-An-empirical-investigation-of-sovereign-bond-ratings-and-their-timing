siop <- function(Yoprob,Yprob=NULL,Ylevel=NULL,X,Z, posc.low.oprob, posc.high.oprob, corr = FALSE,conf = "none",method="BFGS",reps = NULL,b.constant.pool=TRUE,Nobs=NULL,agencies=NULL) {

  # SIOP estimates a selection inflated ordered probit as described in El-Shagi and von Schweinitz (2021). 
  # The difference between SIOP and the closely related ZIOP (zero inflated) model by Brooks and Spencer (2007) is that we observe the time when a reevaluation decision is taken
  #
  # Yoprob          observed outcome of ordered probit equation
  # Yprob           observed outcome of probit equation. If not provided, the model assumes a possible reevaluation at all times and estimates a ZIOP model
  # X               matrix of exogenous variables for probit equation (no constant column)
  # Z               matrix of exogenous variables for ordered probit equation (no constant column)
  # posc.low.oprob  lower threshold of every observation in the ordered probit equation, potentially accounting for boundary conditions
  # posc.high.oprob upper threshold of every observation in the ordered probit equation, potentially accounting for boundary conditions
  # corr            (boolean) determines whether equation errors may be correlated or not. WARNING: corr=TRUE takes considerably longer
  # conf            Confidence bounds using a hessian matrix ("hess"), bootstrapping ("boot") or omits confidence bounds (default = "none")
  # method          method used in optimx. Default = "BFGS"
  # reps            repetitions used in bootstrap (default = 100)
  # b.constant.pool (boolean) defining of constants are pooled or not
  # Nobs            vector of observations by agency
  # agencies        Agency names
  #
  # last changed: 09/14/2021
  #
  # author: mei, gsz
  

  # set calculation method for confidence bounds
  if (conf == "none") {
    hess = FALSE
    boot = FALSE
  }
  if (conf == "boot") {
    hess = FALSE
    boot = TRUE
    stop("Bootstrap not yet implemented")
  }
  if (conf == "hess") {
    hess = TRUE
    boot = FALSE
  }
  
  
  if ((!boot) && (!is.null(reps))) {warning("reps is only used when bootstrapping.")}
  if (is.null(reps)) {reps = 100}
  
  if (is.null(Nobs)){Nobs <- length(Yoprob)}
  if (b.constant.pool){Nobs <- sum(Nobs)}
  
  if (!is.null(Yprob)){
    pos <- Yprob==1
  }else{
    Yprob <- abs(Yoprob)
    pos <- rep(TRUE,length(Yprob))
    print("No separate information on Yprob. Setting Yprob to Yoprob")
  }
  
  if (method == "all"){
    cl.optim <- list(kkt=FALSE,follow.on=TRUE,trace=0,maxit=50000,all.methods = TRUE,starttests=FALSE)
    method=c("BFGS")
  }else{
    cl.optim <- list(kkt=FALSE,follow.on=TRUE,trace=0,maxit=50000,starttests=FALSE)
    method=method
  }
  ##################################
  # Set up matrix of probit constants and oprob thresholds
  pos_startend <- cumsum(c(0,Nobs))
  constmat <- matrix(0,sum(Nobs),length(Nobs))
  for (k in 1:length(Nobs)){constmat[(pos_startend[k]+1):pos_startend[k+1],k]<-1}
  
  if (length(Nobs)==1){
    # one agency or common thresholds
    start_oprob <- setNames(c(rep(0,dim(Z)[2]),-1,1),c(colnames(X),"X0.1","X.1.0"))
    N.const.prob <- 1
    N.const.oprob <- 2
    colnames.p <- c(-1,0,1)
    colnames(constmat) <- "Intercept"
  }else{
    # multiple agencies with non-pooled thresholds
    start_oprob <- setNames(c(rep(0,dim(Z)[2]),rep(c(-1,1),length(Nobs))),
                            c(colnames(X),unlist(lapply(agencies,FUN=function(x) paste0(c("X0.1.","X.1.0."),x)))))
    N.const.prob <- length(Nobs)
    N.const.oprob <- 2*length(Nobs)
    colnames.p <- rep(list(c(-1,0,1)),3)
    colnames(constmat) <- paste("Intercept",agencies,sep=".")
  }
  
  ################################
  # Estimate model equations separately
  print("Estimating probit model")
  data_glm = as.data.frame(cbind(Yprob,X,constmat)) 
  rhs_glm = colnames(data_glm)[-1]
  formula_glm = as.formula(paste("Yprob  ~",paste(rhs_glm,sep="",collapse = "+"), "- 1"))
  res_glm <- glm(data=data_glm,formula=formula_glm,family =  binomial(link = "probit"))
  prob_coeff <- res_glm$coeff
  
  pos.names <- names(start_oprob) %in% names(prob_coeff)
  names(start_oprob)[pos.names] <- paste(names(start_oprob)[pos.names],".1",sep="")
  print("SIOP: Estimating ordered probit model on announcement observations")
  out_oprob <- optimx(start_oprob, fn=LL_oprob_bounded, gr=grad_oprob_bounded, hess=hess_oprob_bounded, 
                       X=Z[pos,], posc.low=posc.low.oprob[pos], posc.high=posc.high.oprob[pos], N.const=N.const.oprob,
                       method=method,control=cl.optim)
  row <- which(out_oprob$value == min(out_oprob$value))
  oprob_coeff <- unlist(out_oprob[row[1],names(start_oprob)])
  
  # Calculate ordered probit model on all observations
  print("ZIOP: Estimating ordered probit model on all observations")
  out_polr <- optimx(oprob_coeff, fn=LL_oprob_bounded, gr=grad_oprob_bounded, hess=hess_oprob_bounded, 
                     X=Z, posc.low=posc.low.oprob, posc.high=posc.high.oprob, N.const=N.const.oprob,
                     method=method,control=cl.optim)
  row <- which(out_polr$value == min(out_polr$value))
  polr_coeff <- unlist(out_polr[row[1],names(oprob_coeff)])
  
  if (corr){
    # estimate probit and oprob jointly
    print("Estimating ZIOP model with correlated errors")
    start_coeff = setNames(c(prob_coeff,oprob_coeff,0),c(names(prob_coeff),names(oprob_coeff),"rho"))
    out_ziopc = optimx(start_coeff, LL_ziopc_bounded, Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,X=X,Z=Z,posc.low=posc.low.oprob, posc.high=posc.high.oprob,Nobs=Nobs,
                       control=cl.optim)
    row <- which(out_ziopc$value == min(out_ziopc$value))
    opt_coeff <- as.matrix(out_ziopc[row[1],1:(length(start_coeff)-1)])
    rho <- as.matrix(out_ziopc[row[1],length(start_coeff)])
  }else{
    rho <- 0
    opt_coeff <- as.matrix(t(c(prob_coeff,oprob_coeff)))
  }
  
  # get optimum from all methods and calculate all output
  res_opt <- LL_ziopc_bounded(opt_coeff,Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,X=X,Z=Z,posc.low=posc.low.oprob, posc.high=posc.high.oprob,Nobs=Nobs,onlyLL=FALSE)
  res_polr <- LL_oprob_bounded(polr_coeff, X=Z, posc.low=posc.low.oprob, posc.high=posc.high.oprob,N.const=N.const.oprob, Nobs = Nobs, colnames.p = colnames.p, onlyLL = FALSE)
  
  fv.oprob <- res_polr$fitted.values
  colnames(fv.oprob) <- sort(unique(Yoprob))
  fv.oprob.inflated <- cbind(res_opt$fitted.values[,1],res_opt$fitted.values[,2]+res_opt$fitted.values[,3],res_opt$fitted.values[,4])
  colnames(fv.oprob.inflated) <- sort(unique(Yoprob))
  LL.oprob <- LL_oprob_bounded.probs(fv.oprob,Yoprob,Ylevel)
  LL.oprob.inflated<- LL_oprob_bounded.probs(fv.oprob.inflated,Yoprob,Ylevel)
  
  print(res_opt$LL)  
  out <- list()
  out$LL <- res_opt$LL                # max likelihood
  out$p_case <- res_opt$p_case        # outcome probabilities
  out$opt_coeff <- opt_coeff          # optimal coefficients
  out$prob_coeff <- prob_coeff        # coefficients of probit model
  out$oprob_coeff <- oprob_coeff      # coefficients of ordered probit model on announcement observations
  out$polr_coeff <- polr_coeff        # coefficients of ordered probit model on all observations
  out$rho <- rho                      # correlation
  out$fitted.values <- res_opt$fitted.values  # latent variables
  out$latent.prob <- res_opt$latent.prob
  out$latent.oprob <- res_opt$latent.oprob
  out$LL.oprob <- LL.oprob
  out$LL.oprob.inflated <- LL.oprob.inflated
  out$p.chi2 <- 1-pchisq(2*(LL.oprob.inflated-LL.oprob),length(opt_coeff)-length(res_polr$coefficients)-2)
  
  if (hess){
    print("Calculating Hessian") 
    if (corr){
      out$hessian <- hessian(LL_ziopc_bounded,opt_coeff,Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,X=X,Z=Z,posc.low=posc.low.oprob, posc.high=posc.high.oprob,Nobs=Nobs)
    }else{
      hX <- ginv(vcov(glm(data=data_glm,formula=formula_glm,family =  binomial(link = "probit"))))
      hZ <- hess_oprob_bounded(oprob_coeff,X=Z[pos,], posc.low=posc.low.oprob[pos], posc.high=posc.high.oprob[pos], N.const=N.const.oprob)
      out$hessian <- as.matrix(bdiag(hX,hZ))
      out$hessian_polr <- hess_oprob_bounded(polr_coeff,X=Z[pos,], posc.low=posc.low.oprob[pos], posc.high=posc.high.oprob[pos], N.const=N.const.oprob)
    }

    sd <- sqrt(diag(ginv(out$hessian)))
    out$p <- pt(opt_coeff/sd,df=length(Yoprob)-length(opt_coeff))
    out$p <- 2*pmin(out$p,1-out$p)
    out$sd=sd
    
  }
  
  return(out)
}