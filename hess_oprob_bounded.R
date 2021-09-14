hess_oprob_bounded <- function(coeff, X, posc.low, posc.high,N.const=2) {
  #keep polr-coefficients as they are consistent with http://people.stern.nyu.edu/wgreene/OrderedChoiceModeling.pdf, p.142
  # posc.low and posc.high give the position in the vector of 
  # constants c(-Inf,c_1,...,Inf) that give the upper and lower bound
  # This approach is quite flexible and allows for instances that cover multiple classes
  # (as in the boundary adjustment of this paper)
  
  lc = length(coeff)
  k = dim(X)[2]
  T <- dim(X)[1]
  if (lc!=k+N.const){stop("dimensions (X and coeffs) do not match")}
  constants <- c(-Inf,coeff[k+1:N.const],Inf)
  
  G = coeff[1:k]
  y_star = X %*% (G)
  pointmat <- cbind(constants[posc.low]-y_star,constants[posc.high]-y_star)
  Fmat <- pnorm(pointmat)
  Pvec <- Fmat[,2]-Fmat[,1]
  dmat <- dnorm(pointmat)
  ddmat <- ddnorm(pointmat)
  ddmat[is.nan(ddmat)] <-  0

  hess.mat <- matrix(0,lc,lc)
  for (t in 1:T){
    # Hessian of beta v beta
    hess.mat[1:k,1:k] <- hess.mat[1:k,1:k]+
      ((ddmat[t,2]-ddmat[t,1])/Pvec[t]-((dmat[t,2]-dmat[t,1])/Pvec[t])^2)*(X[t,]%*%t(X[t,]))
  }
  for (j in 1:N.const){
    phigh <- posc.high==j+1
    plow <- posc.low==j+1
    # Hessian of beta v constants
    hess.mat[1:k,k+j] <- colSums(-X[phigh,]*(ddmat[phigh,2]/Pvec[phigh] - (dmat[phigh,2]-dmat[phigh,1])*dmat[phigh,2]/(Pvec[phigh]^2))) + 
      colSums(-X[plow,]*(-ddmat[plow,1]/Pvec[plow]-(dmat[plow,2]-dmat[plow,1])*(-dmat[plow,1])/(Pvec[plow]^2)))
    hess.mat[k+j,1:k] <- t(hess.mat[1:k,k+j])
    
    # Hessian of constants v constants
    hess.mat[k+j,k+j] <- sum(ddmat[phigh,2]/Pvec[phigh]-(dmat[phigh,2]/Pvec[phigh])^2)+
      sum(-ddmat[plow,1]/Pvec[plow]-(-dmat[plow,1]/Pvec[plow])^2)
    
    for (j2 in setdiff(1:j,j)){
      poverlap <- (posc.low==j2+1) & phigh
      hess.mat[k+j,k+j2] <- sum(dmat[poverlap,2]*dmat[poverlap,1]/(Pvec[poverlap]^2))
      hess.mat[k+j2,k+j] <- hess.mat[k+j,k+j2]
    }
  }
  return(-hess.mat)
}