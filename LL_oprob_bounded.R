LL_oprob_bounded <- function(coeff, X, posc.low, posc.high,N.const=2, Nobs = NULL, colnames.p = NULL, onlyLL = TRUE) {
  # Likelihood of an ordered probit model
  # keep polr-coefficients as they are consistent with http://people.stern.nyu.edu/wgreene/OrderedChoiceModeling.pdf, p.142
  #################################
  # INPUT
  #   coeff     coefficient vector of size dim(X)[2] + N.const
  #   X         matrix of covariates
  #   posc.low  lower bound constant of latent variable in c(-Inf,c_1,...,c_N.const,Inf)
  #   posc.high upper bound constant of latent variable in c(-Inf,c_1,...,c_N.const,Inf)
  #   N.const   number of constants = number of observed states in the ordered probit
  #   onlyLL    Boolean if only log-likelihood should be returned

  lc = length(coeff)
  N = dim(X)[1]
  k = dim(X)[2]
  if (length(N.const)>1){stop("Differentiate number of constants per agency via coefficient names")}
  if (lc!=k+N.const){stop("dimensions (X and coeffs) do not match")}
  constants <- c(-Inf,coeff[k+1:N.const],Inf)
  G = coeff[1:k]
  y_star = X %*% (G)
  oprobmat <- cbind(constants[posc.low]-y_star,constants[posc.high]-y_star)   # mu_j - beta'x_j , mu_j-1 - beta'x_j
  Fmat <- pnorm(oprobmat)
  Pvec <- Fmat[,2]-Fmat[,1]
  LL <- -sum(log(Pvec))
  
  if (onlyLL) {return(LL)}
  
  #############################
  # Case else: onlyLL = FALSE
  if (is.null(Nobs)){Nobs <- dim(X)[1]}
  if (is.null(colnames.p)){colnames.p <- c(1:ceiling(N.const/length(Nobs)))}
  if (!is.list(colnames.p)){colnames.p <- rep(list(colnames.p),length(Nobs))}
  
  colnames.p.unique <- sort(unique(unlist(colnames.p)))
  N.const.max <- length(colnames.p.unique)
  cols.count <- 1:N.const.max
  lN <- length(Nobs)
  if (lN>1){
    thresh.exists <- matrix(unlist(lapply(colnames.p, FUN = function(x) setdiff(colnames.p.unique,max(colnames.p.unique)) %in% x)),N.const.max-1,lN)
    const.pos <- 1+apply(thresh.exists,2,cumsum)
    for (col in 2:lN){
      const.pos[const.pos[,col]>1,col] <- const.pos[const.pos[,col]>1,col] + const.pos[N.const.max-1,col-1]-1
    }
  }else{
    # pos.shift <- 0
    thresh.exists <- matrix(TRUE,N.const.max-1,lN)
    const.pos <- 1+apply(thresh.exists,2,cumsum)
  }
  
  oprobmat <- rep(-Inf,sum(Nobs))
  for (thresh in 1:N.const.max-1){
    pos.base <- const.pos[thresh,]  # Some agencies may have fewer rating classes than others. But we assume here that differences between agencies arise at low rating classes
    pos.by.ag <- rep(pos.base,Nobs)
    oprobmat <- cbind(oprobmat,constants[pos.by.ag] - y_star)
  }
  oprobmat <- cbind(oprobmat,rep(Inf,sum(Nobs)))
  p <- pnorm(oprobmat[,2:(N.const.max+1)]) - pnorm(oprobmat[,1:(N.const.max)])
  colnames(p) = colnames.p.unique
  
  res = list()
  res$LL = -LL
  res$p_case <- Pvec
  res$fitted.values = p
  res$latent.oprob <- oprobmat
  res$oprob_coeff <- c(coeff)
  return(res)
}