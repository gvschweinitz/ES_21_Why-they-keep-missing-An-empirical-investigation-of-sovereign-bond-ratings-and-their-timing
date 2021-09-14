hess_prob <- function(coeff, Y, X,Nobs) {
  if (length(unique(Y))==3){Y<-abs(Y)}
  
  lc = length(coeff)
  k_x = dim(X)[2]
  T <- dim(X)[1]
  N.const <- length(Nobs)
  if (lc!=k_x+N.const){stop("dimensions (X and coeffs) do not match")}
  
  X <- cbind(X,as.matrix(bdiag(apply(matrix(Nobs),1,FUN=function(x){return(rep(1,x))}))))
  r_star <- X%*%coeff
  
  dmat <- dnorm(r_star) * (
            Y*(dnorm(r_star)/(pnorm(r_star)*pnorm(r_star))+r_star/pnorm(r_star))- 
            (1-Y)*(dnorm(r_star)/(pnorm(-r_star)*pnorm(-r_star))-r_star/pnorm(-r_star)))

  hess.mat <- matrix(0,lc,lc)
  for (t in 1:T){
    hess.mat <- hess.mat+dmat[t]*(X[t,]%*%t(X[t,]))
  }
  # Npos <- cumsum(c(0,Nobs))
  # for (j in 1:N.const){
  #   pos <- (Npos[j]+1):Npos[j+1]
  #   hess.mat[1:k_x,k_x+j] <- colSums(X[pos,] * matrix(dmat[pos],length(pos),k_x,byrow=FALSE))
  #   hess.mat[k_x+j,1:k_x] <- t(hess.mat[1:k_x,k_x+j])
  # }
  
  return(-hess.mat)
}