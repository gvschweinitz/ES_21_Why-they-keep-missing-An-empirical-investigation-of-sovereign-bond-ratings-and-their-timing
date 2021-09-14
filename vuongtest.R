vuongtest  <- function(LLcont_F,LLcont_G,grad_F,grad_G,hessian_F,hessian_G,ndraws,alpha){
  
  # Conducts a non-degenerate Vuong-test (Shi, Quantitative Economics, 2015; Vuong, Ectra 1989) which compares potentially overlapping models.
  # Matlab code is provided from the website of Xiaoxia Shi (https://www.ssc.wisc.edu/~xshi/research.html)
  # Requires from each model optimal coefficients, log-likelihoods, gradient contributions of individual observations and hessian matrices
  # Gregor von Schweinitz, 2021-05-03
  
  # Step 1: input arguments
  
  require(expm)


  n <- dim(grad_F)[1]
  d_F <- dim(grad_F)[2]
  d_G <- dim(grad_G)[2]
  k <- d_F + d_G
  
  # Step 2: Ahat and Bhat  
  Ahat <- as.matrix(bdiag(hessian_F,-hessian_G))
  
  d_log <- cbind(grad_F,-grad_G)
  meandL <- colMeans(d_log)
  dmeandL <- d_log - meandL
  Bhat <- t(dmeandL) %*% dmeandL/n + 1e-12*diag(1,k,k)
  sqrtmBhat <- sqrtm(Bhat)
  
  What <- (sqrtmBhat %*% ginv(Ahat)) %*% sqrtmBhat
  What <- 1/2*(t(What) + What)
  V <- eigen(What,only.values = TRUE)$values
  
  # Step 3
  # missing in the code: choose c and calculate \hat{T}_n^{mod}(c) according to eq. (4.3)
  
  # Step 4: simplified as rho_star takes on extreme values at maximum eigenvalue
  Z0 <- matrix(rnorm(ndraws*(k+1)),ndraws,k+1)
  rho_star <- 1*(abs(V) == max(abs(V)))
  VZ <- rbind(matrix(c(1,rho_star),1,k+1),
              cbind(matrix(rho_star,k,1),diag(1,k,k)))   # diag(1,k+1,k+1) in code, diag(V) in paper -> normalization?
  Z <- Z0 %*% sqrtm(VZ + 1e-12*diag(1,k+1,k+1))
  Z_L <- Z[,1]
  Z_p <- Z[,2:(k+1)]
  
  trVsq <- sum(V^2)
  Vnmlzd <- V/sqrt(trVsq)
  
  # J_Lmod <- function(sig,c){
  #   J_Lmod <- sig*Z_L - (Z_p^2) %*% Vnmlzd/2 + sum(Vnmlzd)/2
  #   return(J_Lmod)
  # }
  # 
  # J_omod <- function(sig,c){
  #   J_omod <- sig^2 - 2*sig*Z_p%*%diag(Vnmlzd)%*%rho_star + (Z_p^2) %*% Vnmlzd^2 + c
  #   return(J_omod)
  # }
  
  # Step 5-6: Computation of J_s(sig,rho_star,V,const), and the 1-alpha quantile thereof
  J_Lconst <- - (Z_p^2) %*% Vnmlzd/2 + sum(Vnmlzd)/2
  J_oconst <- (Z_p^2) %*% Vnmlzd^2
  J_fac <- -Z_p%*%diag(Vnmlzd)%*%rho_star
  
  quant <- function(sig,const,Z_L,J_fac,J_Lconst,J_oconst){
    J_L <- sig*Z_L + J_Lconst
    J_o <- sig^2 + 2*sig*J_fac + J_oconst + const
    # return(quantile(abs(J_Lmod(sig,c) / sqrt(J_omod(sig,c))),1-alpha))
    return(quantile(abs(J_L / sqrt(J_o)),1-alpha))
  }
  
  # Step 7: calculate supremum of quantiles over sigmas for given const
  sigstar <- function(const,Z_L,J_fac,J_Lconst,J_oconst){
    res <- optimx(1,fn=quant,lower=0,upper=5,const=const,Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst,
                  method="L-BFGS-B",control=list(maximize=TRUE,trace=0))
    return(list(sig=res[,1],critval=res[,"value"]))
  }
  
  sigstar_vec <- function(constvec,Z_L,J_fac,J_Lconst,J_oconst){
    sig <- rep(0,length(constvec))
    critval <- rep(0,length(constvec))
    for (k in 1:length(constvec)){
      res <- optimx(1,fn=quant,lower=0,upper=5,const=constvec[k],Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst,
                    method="L-BFGS-B",control=list(maximize=TRUE))
      sig[k] <- res[,1]
      critval[k] <- res[,"value"]
    }
    
    return(list(sig=sig,critval=critval))
  }
  
  quantdiff <- function(const,qshift,Z_L,J_fac,J_Lconst,J_oconst){
    outconst <- sigstar(const=const,Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst)
    qdiff <- (outconst$critval-qshift)^2
    print(c(const,outconst$critval,qdiff))
    return(qdiff)
  }
  
  #   (a) comparison to normal distribution -> if this is a close approximation, stick with const = 0
  out0 <- sigstar(0,Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst)
  znormal <- qnorm(1-alpha/2)
  znormal_sim <- max(znormal,quantile(abs(Z_L),1-alpha))
  
  if (out0$critval - znormal_sim <= 0.1){
    critval <- max(out0$critval,znormal)
    conststar <- 0
  }else{
    print("optimizing constant")
    qshift <- znormal_sim+0.1
    lb <- 0
    ub <- 10
    nsteps <- 11
    for (step in 1:10){
      grid <- seq(lb,ub,length.out = nsteps)
      out_step <- sigstar_vec(grid,Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst)
      qdiff <- (out_step$critval-qshift)^2
      lb <- grid[max(1,which.min(qdiff)-1)]
      ub <- grid[min(nsteps,which.min(qdiff)+1)]
    }
    
    conststar <- grid[which.min(grid)]
    outc <- sigstar(conststar,Z_L=Z_L,J_fac=J_fac,J_Lconst=J_Lconst,J_oconst=J_oconst)
    critval <- max(outc$critval,znormal)
  }
  
  # Computation of LR statistic
  nLR_hat <- sum(LLcont_F-LLcont_G)
  
  nomega2_hat = sum((LLcont_F-LLcont_G)^2) - nLR_hat^2/n
  
  # Non-degenerate Vuong test statistic
  Tnd <- (nLR_hat + sum(V)/2) / (sqrt(nomega2_hat) + conststar*sum(V^2))
  
  return(list(conststar = conststar, critval = critval, Tnd = Tnd))
}