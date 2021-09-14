grad_prob <- function(coeff, Y, X,Nobs,do.contribution = FALSE) {
  
  if (length(unique(Y))==3){Y<-abs(Y)}
  
  lc = length(coeff)
  k_x = dim(X)[2]
  T <- dim(X)[1]
  N.const <- length(Nobs)
  if (lc!=k_x+N.const){stop("dimensions (X and coeffs) do not match")}
  
  X <- cbind(X,as.matrix(bdiag(apply(matrix(Nobs),1,FUN=function(x){return(rep(1,x))}))))
  r_star <- X%*%coeff
  # B = coeff[1:k_x]
  # con_pro = rep(coeff[(k_x+1):lc],Nobs)
  # r_star = X %*% (B)+con_pro
  
  dmat = Y*dnorm(r_star)/pnorm(r_star)- (1-Y)*dnorm(r_star)/pnorm(-r_star)
  grad.contribution <- X * matrix(dmat,T,lc,byrow=FALSE)
  grad <- colSums(grad.contribution)
  if (do.contribution){
    return(list(grad = -grad, grad.contribution = -grad.contribution))
  }else{
    return(-grad)
  }
  
  # grad.const <- rep(0,N.const)
  # Npos <- cumsum(c(0,Nobs))
  # for (j in 1:N.const){
  #   grad.const[j] <- sum(dmat[(Npos[j]+1):Npos[j+1]])
  # }
  # print(-c(grad.beta,grad.const))
  # return(c(grad.beta,grad.const))
}