grad_oprob_bounded <- function(coeff, X, posc.low, posc.high,N.const=2,do.contribution = FALSE) {
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
  pointmat <- cbind(constants[posc.low]-y_star,constants[posc.high]-y_star)     # mu_j - beta'x_j , mu_j-1 - beta'x_j
  Fmat <- pnorm(pointmat)
  Pvec <- Fmat[,2]-Fmat[,1]
  dmat <- dnorm(pointmat)
  
  grad.contribution.beta <- -X * matrix((dmat[,2]-dmat[,1])/Pvec,T,k,byrow=FALSE)
  
  
  grad.contribution.const <- matrix(0,T,sum(N.const))
  for (j in 1:N.const){
    phigh <- posc.high==j+1
    plow <- posc.low==j+1
    grad.contribution.const[phigh,j] <- dmat[phigh,2]/Pvec[phigh]
    grad.contribution.const[plow,j] <- -dmat[plow,1]/Pvec[plow]
  }
  
  grad.beta <- colSums(grad.contribution.beta)
  grad.const <- colSums(grad.contribution.const)
  
  if (do.contribution){
    return(list(grad = -c(grad.beta,grad.const), grad.contribution = -cbind(grad.contribution.beta,grad.contribution.const)))
  }else{
    return(-c(grad.beta,grad.const))
  }
  
}