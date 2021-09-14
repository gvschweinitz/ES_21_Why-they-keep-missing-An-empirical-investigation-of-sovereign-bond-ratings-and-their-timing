# SCRIPT for conducting different specification tests for the SIOP for rating data
# Output used for: Tables 1, 4 and 5

# last changed: 9/14/2021
#
# author: mei, gsz

wd <- getwd()
source(paste(wd,"/siop.R",sep=""))
source(paste(wd,"/fun_collection_derivatives.R",sep=""))
source(paste(wd,"/grad_oprob_bounded.R",sep=""))
source(paste(wd,"/hess_oprob_bounded.R",sep=""))
source(paste(wd,"/grad_prob.R",sep=""))
source(paste(wd,"/hess_prob.R",sep=""))
source(paste(wd,"/LL_oprob_bounded.R",sep=""))
source(paste(wd,"/LL_oprob_bounded.probs.R",sep=""))
source(paste(wd,"/LL_ziopc_bounded.R",sep=""))
source(paste(wd,"/LL_ziopc_bounded_null.R",sep=""))
source(paste(wd,"/res_addpseudoR2.R",sep=""))
source(paste(wd,"/vuongtest.R",sep=""))

packages <- c("optimx","MASS","Matrix","matrixcalc","miscTools","pbivnorm","numDeriv","R.utils","cgwtools")
for (p in packages){
  check = require(p,character.only=TRUE)
  if (!check) {
    install.packages(p)
    require(p,character.only=TRUE)
  }  
}

savefolder <- wd
keepvars.alltests <- c(ls(),"keepvars.alltests","output.table")

evalcrit <- c("m1","m2","#obs","#coeffs (m1)","#coeffs (m2)","LL (m1, 3 case/4 case)","LL (m2, 3 case/4 case)","LRstat","p(LR-test)","Hstat","p(Hausman)","Vuongstat","Vuongcritval (10%)")
output.table <- as.data.frame(matrix("",0,length(evalcrit),dimnames=list(NULL,evalcrit)))

##########################################
##########################################
##########################################
# Table 4: Model comparison


##########################################
# Model comparison, step 1: LR test of MIOP against oprob; Hausman test of SIOP against MIOP
folders <- c("","../R trend/")
names.add <- c("","trend")
mod <- "baseline"
nmodels.old <- dim(output.table)[1]
nmodels <- 2*length(folders)
output.table <- rbind(output.table,matrix(NA,nmodels,dim(output.table)[2],
                                          dimnames=list(nmodels.old + (1:nmodels),dimnames(output.table)[[2]])))

keepvars = c(ls(),"pos","keepvars")

for (pos in 1:length(folders)){
  modfile<- paste(folders[pos],mod,".pool.RData",sep="")
  out <- loadToEnv(modfile)    #[c("res","resfund","resziop","resziop.fund")]
  res <- out[["res"]]
  resziop <- out[["resziop"]]
  
  if (is.null(resziop)){
    print("calculating full ziop (may take some time)")
    Yprob <- out[["Yprob"]]
    Yoprob <- out[["Yoprob"]]
    Ylevel <- out[["Ylevel"]]
    X <- out[["X"]]
    Z <- out[["Z"]]
    
    kX <- dim(X)[2]
    kZ <- dim(Z)[2]
    posc.low.oprob <- out[["posc.low.oprob"]]
    posc.high.oprob <- out[["posc.high.oprob"]]
    Nobs <- out[["Nobs"]]
    corr <- out[["corr"]]
    
    allcoeff_ziop = optimx(c(res$opt_coeff), LL_ziopc_bounded,method=c("BFGS"),Yoprob=Yoprob,X=X,Z=Z,posc.low=posc.low.oprob,posc.high=posc.high.oprob,Nobs=Nobs,
                           control=list(kkt=FALSE,follow.on=TRUE,maxit=1000,trace=0,reltol=1e-6))
    ziopcoeffs <- as.numeric(allcoeff_ziop[which.min(allcoeff_ziop[,"value"]),1:(dim(X)[2]+dim(Z)[2]+3*length(Nobs)+corr)])
    resziop <- LL_ziopc_bounded(ziopcoeffs,Yoprob=Yoprob,X=X,Z=Z,posc.low=posc.low.oprob,posc.high=posc.high.oprob,Nobs=Nobs,onlyLL=FALSE)
    resziop$hessian <- hessian(LL_ziopc_bounded,ziopcoeffs,Yoprob=Yoprob,X=X,Z=Z,posc.low=posc.low.oprob, posc.high=posc.high.oprob,Nobs=Nobs)
    resziop$LL.SIOP <- -LL_ziopc_bounded(ziopcoeffs,Yoprob=Yoprob,Yprob=Yprob,X=X,Z=Z,posc.low=posc.low.oprob,posc.high=posc.high.oprob,Nobs=Nobs,onlyLL=TRUE)  # adding 4-case likelihood
    resave(resziop,file=modfile)
  }else{
    ziopcoeffs <- c(resziop$prob_coeff,resziop$oprob_coeff)
  }
  
  # LR test of MIOP against oprob
  LRstat <- 2*(resziop$LL-res$LL.oprob)
  p.LR <- 1-pchisq(LRstat,length(ziopcoeffs)-length(res$polr_coeff))
  output.table[nmodels.old+2*(pos-1)+1,] <- c(paste("MIOP",names.add[pos]),paste("oprob",names.add[pos]),length(out$Yoprob),
                                              length(ziopcoeffs),length(res$polr_coeff),
                                              round(resziop$LL,digits=2),round(res$LL.oprob,digits=2),
                                              round(LRstat,digits=2),round(p.LR,digits=2),rep(NA,4))
  
  # Hausman test of SIOP against MIOP
  H_res <- res$hessian
  H_resziop <- -resziop$hessian
  Hdiff <- H_res - H_res %*% ginv(H_resziop + H_res) %*% H_res
  Hstat <- (c(res$opt_coeff)-ziopcoeffs) %*% Hdiff %*% (c(res$opt_coeff)-ziopcoeffs)
  p.Hausman <- 1-pchisq(Hstat,length(res$opt_coeff)-1)
  output.table[nmodels.old+2*(pos-1)+2,] <- c(paste("SIOP",names.add[pos]),paste("MIOP",names.add[pos]),length(out$Yoprob),length(res$opt_coeff),length(ziopcoeffs),
                                              paste(round(res$LL.oprob.inflated,digits=2),round(res$LL,digits=2),sep="/"),
                                              paste(round(resziop$LL,digits=2),round(resziop$LL.SIOP,digits=2),sep="/"),
                                              rep(NA,2),round(Hstat,digits=2),round(p.Hausman,digits=2),rep(NA,2))
  rm(list=setdiff(ls(),keepvars))
}

print("Writing output")
# format(output.table,digits=3)
write.csv(output.table,file=paste(savefolder,"/LLR_tests_new.csv",sep=""),na="",row.names = FALSE)

rm(list=setdiff(ls(),keepvars.alltests))



##########################################
# Model comparison, step 2: Compare model to alternative setups ("submodels") with the same variables
# Table 4 uses only the first comparison of these (baseline vs. baseline fundamentals)
estimations <- c("baseline","NoRatingsq","NoDynamics","NoSpilloverCompetitor","NoCompetitors","NoSpillovers",
                 "RobustDefHist","political_full","Homebias_short","Homebias_long","oecd","eu","outlook","baseline.outlook")

conf.LL <- "hess"
b.constant.pool <- FALSE

nmodels.old <- dim(output.table)[1]
submodels <- c("fundamentals")
models <- expand.grid(estimations,submodels)
models <- models[order(models[,1]),]
nmodels <- dim(models)[1]
output.table <- rbind(output.table,matrix(NA,nmodels,dim(output.table)[2],
                                          dimnames=list(nmodels.old + (1:nmodels),dimnames(output.table)[[2]])))

output.table[nmodels.old + (1:nmodels),c(1:2)] <- cbind(as.character(models[,1]),as.character(models[,2]))

keepvars = c(ls(),"keepvars","pos1")

for (pos1 in 1:length(estimations)) {
  m1 <- estimations[pos1]
  print(m1)
  comparisons <- models[models[,1]==m1,2]
  out <- loadToEnv(paste(wd,"/",m1,".pool.RData",sep=""))#[c("res","resfund","resziop","resziop.fund")]
  res <- out[["res"]]
  
  Yprob <- out[["Yprob"]]
  Yoprob <- out[["Yoprob"]]
  Ylevel <- out[["Ylevel"]]
  X <- out[["X"]]
  Z <- out[["Z"]]
  
  kX <- dim(X)[2]
  kZ <- dim(Z)[2]
  posc.low.oprob <- out[["posc.low.oprob"]]
  posc.high.oprob <- out[["posc.high.oprob"]]
  Nobs <- out[["Nobs"]]
  corr <- out[["corr"]]
  kbase <- kX+kZ+3*length(Nobs)+corr
  H_res <- res$hessian
  var_res <- ginv(H_res)
  if (corr){stop("not adjusted for correlated models")}
  
  if ("fundamentals" %in% comparisons){
    print("fundamentals")
    usefund <- c(out[["colfund"]],"default")
    kXfund <- length(intersect(usefund,colnames(X)))
    kZfund <- length(intersect(usefund,colnames(Z)))
    kfund <- kXfund+kZfund+3*length(Nobs)+corr
    
    resfund <- out[["resfund"]]
    if (is.null(resfund)){
      Xfund <- X[,intersect(usefund,colnames(X))]
      Zfund <- Z[,intersect(usefund,colnames(Z))]
      resfund <- siop(Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,X=Xfund,Z=Zfund,posc.low.oprob=posc.low.oprob,posc.high.oprob=posc.high.oprob,corr=corr,
                          conf=conf.LL,method=c("BFGS"),b.constant.pool=b.constant.pool,Nobs=Nobs,
                          agencies=out[["agencies"]])
      resave(resfund,file=paste(wd,"/",m1,".pool.RData",sep=""))
    }
    LRstat <- 2*(res$LL-resfund$LL)
    p.LR <- 1-pchisq(LRstat,kbase-kfund)
    
    coeff_comp <- 0*res$opt_coeff
    m1.coeffpos.comp <- match(colnames(resfund$opt_coeff),colnames(res$opt_coeff))
    coeff_comp[m1.coeffpos.comp] <- resfund$opt_coeff
    # Calculate inverse of difference between variances (https://en.m.wikipedia.org/wiki/Woodbury_matrix_identity)
    m2.H<- -resfund$hessian
    U <- matrix(0,kbase,kfund)
    U[m1.coeffpos.comp,] <- diag(kfund)
    V <- t(U)
    Hdiff <- H_res - (H_res %*% U) %*% ginv(m2.H+ V %*% H_res %*% U) %*% (V %*% H_res)
    Hstat <- (coeff_comp-res$opt_coeff) %*% Hdiff %*% t(coeff_comp-res$opt_coeff)
    p.Hausman <- 1-pchisq(Hstat,kbase-1)
    
    output.table[output.table[,1]==m1 & output.table[,2]=="fundamentals",3:dim(output.table)[2]] <- 
      c(length(Yoprob),kbase,kfund,round(c(res$LL,resfund$LL,LRstat,p.LR,Hstat,p.Hausman),digits=2),NA,NA)
  }
  
  rm(list=setdiff(ls(),keepvars))
  
}

print("Writing output")
write.csv(output.table,file=paste(savefolder,"/LLR_tests_new.csv",sep=""),na="",row.names = FALSE)

rm(list=setdiff(ls(),keepvars.alltests))

##########################################
# Model comparison, step 3: Compare baseline model to alternative smaller specifications
print("Comparing baseline to smaller models")
wd <- getwd()
savefolder <- wd
est.baseline <- "baseline"
estimations <- c("NoRatingsq","NoDynamics","NoCompetitors","NoSpillovers","NoSpilloverCompetitor")
nmodels <- length(estimations)
nmodels.old <- dim(output.table)[1]
output.table <- rbind(output.table,matrix(NA,nmodels,dim(output.table)[2],
                                          dimnames=list(nmodels.old + (1:nmodels),dimnames(output.table)[[2]])))


m1 <- est.baseline
m1.res <- loadToEnv(paste(wd,"/",m1,".pool.RData",sep=""))[["res"]]
Yoprob <- loadToEnv(paste(wd,"/",m1,".pool.RData",sep=""))[["Yoprob"]]
k1 <- length(m1.res$opt_coeff)
m1.H <- m1.res$hessian
m1.var <- ginv(m1.res$hessian)

keepvars = c(ls(),"keepvars","pos2")
for (pos2 in 1:nmodels){
  m2 <- estimations[pos2]
  m2.res <- loadToEnv(paste(wd,"/",m2,".pool.RData",sep=""))[["res"]]
  k2 <- length(m2.res$opt_coeff)
  LRstat <- 2*(m1.res$LL-m2.res$LL)
  p.LR <- 1-pchisq(LRstat,k1-k2)
  
  coeff_comp <- 0*m1.res$opt_coeff
  m1.coeffpos.comp <- match(colnames(m2.res$opt_coeff),colnames(m1.res$opt_coeff))
  coeff_comp[m1.coeffpos.comp] <- m2.res$opt_coeff
  # Calculate inverse of difference between variances (https://en.m.wikipedia.org/wiki/Woodbury_matrix_identity)
  m2.H <- -m2.res$hessian
  U <- matrix(0,k1,k2)
  U[m1.coeffpos.comp,] <- diag(k2)
  V <- t(U)
  Hdiff <- m1.H - (m1.H %*% U) %*% ginv(m2.H + V %*% m1.H %*% U) %*% (V %*% m1.H)
  Hstat <- (coeff_comp-m1.res$opt_coeff) %*% Hdiff %*% t(coeff_comp-m1.res$opt_coeff)
  p.Hausman <- 1-pchisq(Hstat,k1-1)
  
  output.table[nmodels.old+pos2,] <- c(m1,m2,length(Yoprob),k1,k2,round(c(m1.res$LL,m2.res$LL,LRstat,p.LR,Hstat,p.Hausman),digits=2),NA,NA)
  
  
  rm(list=setdiff(ls(),keepvars))
}
# }
write.csv(output.table,file=paste(savefolder,"/LLR_tests_new.csv",sep=""),na="",row.names = FALSE)
rm(list=setdiff(ls(),keepvars.alltests))

######################
# Model comparison, step 4: Compare trend model to baseline model
print("Comparing baseline to trend model")

mod <- "baseline"
out <- loadToEnv(paste(wd,"/",mod,".pool.RData",sep=""))#[c("res","res_polr","resfund","res.ziop","res.ziop.fund")]
outtrend <- loadToEnv(paste(file.path("..","R trend"),"/",mod,".pool.RData",sep=""))#[c("res","res_polr","resfund","res.ziop","res.ziop.fund")]
vuong_SIOP_SIOPtrend <- out[["vuong_SIOP_SIOPtrend"]]

res.o <- out[["res"]]
res.t <- outtrend[["res"]]
pos_oint <- out$pos %in% outtrend$pos
pos_tino <- outtrend$pos %in% out$pos
Nobs_joint <- c(sum(pos_oint[1:out$Nobs[1]]),sum(pos_oint[out$Nobs[1]+1:out$Nobs[2]]),sum(pos_oint[out$Nobs[1]+out$Nobs[2]+1:out$Nobs[3]]))

Yoprob.joint <- out$Yoprob[pos_oint]
Yprob.joint <- out$Yprob[pos_oint]
X.o <- out$X[pos_oint,]
X.t <- outtrend$X[pos_tino,]
Z.o <- out$Z[pos_oint,]
Z.t <- outtrend$Z[pos_tino,]
posc.low.oprob <- out$posc.low.oprob[pos_oint]
posc.high.oprob <- out$posc.high.oprob[pos_oint]

posYoprob <- Yprob.joint == 1
fv.o <- res.o$fitted.values[pos_oint,]
fv.t <- res.t$fitted.values[pos_tino,]

LLcont_F <- rep(0,length(Yoprob.joint))
LLcont_G <- rep(0,length(Yoprob.joint))
LLcont_F[!posYoprob] <- log(fv.o[!posYoprob,"00"])
LLcont_G[!posYoprob] <- log(fv.t[!posYoprob,"00"])
LLcont_F[posYoprob & Yoprob.joint==-1] <- log(fv.o[posYoprob & Yoprob.joint==-1,"-1"])
LLcont_G[posYoprob & Yoprob.joint==-1] <- log(fv.t[posYoprob & Yoprob.joint==-1,"-1"])
LLcont_F[posYoprob & Yoprob.joint==0] <- log(fv.o[posYoprob & Yoprob.joint==0,"01"])
LLcont_G[posYoprob & Yoprob.joint==0] <- log(fv.t[posYoprob & Yoprob.joint==0,"01"])
LLcont_F[posYoprob & Yoprob.joint==1] <- log(fv.o[posYoprob & Yoprob.joint==1,"1"])
LLcont_G[posYoprob & Yoprob.joint==1] <- log(fv.t[posYoprob & Yoprob.joint==1,"1"])
if (is.null(vuong_SIOP_SIOPtrend)){
  grad_F <- matrix(0,length(Yoprob.joint),length(res.o$opt_coeff))
  grad_F[,1:length(res.o$prob_coeff)] <- -grad_prob(res.o$prob_coeff,Y = Yprob.joint, X=X.o, Nobs = Nobs_joint,do.contribution = TRUE)$grad.contribution
  grad_F[posYoprob,(length(res.o$prob_coeff)+1):length(res.o$opt_coeff)] <- 
    -grad_oprob_bounded(res.o$oprob_coeff, X=Z.o[posYoprob,], posc.low=posc.low.oprob[posYoprob], posc.high=posc.high.oprob[posYoprob],N.const=2*length(Nobs_joint),do.contribution = TRUE)$grad.contribution
  
  grad_G <- matrix(0,length(Yoprob.joint),length(res.t$opt_coeff))
  grad_G[,1:length(res.t$prob_coeff)] <- -grad_prob(res.t$prob_coeff,Y = Yprob.joint, X=X.t, Nobs = Nobs_joint,do.contribution = TRUE)$grad.contribution
  grad_G[posYoprob,(length(res.t$prob_coeff)+1):length(res.t$opt_coeff)] <- 
    -grad_oprob_bounded(res.t$oprob_coeff, X=Z.t[posYoprob,], posc.low=posc.low.oprob[posYoprob], posc.high=posc.high.oprob[posYoprob],N.const=2*length(Nobs_joint),do.contribution = TRUE)$grad.contribution
  
  hessian_F <- -1/2*(res.o$hessian + t(res.o$hessian))
  hessian_G <- -1/2*(res.t$hessian + t(res.t$hessian))
  
  vuong_SIOP_SIOPtrend <- vuongtest(LLcont_F,LLcont_G,grad_F,grad_G,hessian_F,hessian_G,ndraws=5000,alpha=0.1)
  resave(vuong_SIOP_SIOPtrend,file=paste(wd,"/",mod,".pool.RData",sep=""))
}

nmodels.old <- dim(output.table)[1]
output.table <- rbind(output.table,matrix(NA,1,dim(output.table)[2],
                                          dimnames=list(nmodels.old + 1,dimnames(output.table)[[2]])))
output.table[nmodels.old + 1,] <- c("SIOP","SIOPtrend",length(Yoprob.joint),length(res.o$opt_coeff),length(res.t$opt_coeff),
                                    round(c(sum(LLcont_F),sum(LLcont_G),rep(NA,4),vuong_SIOP_SIOPtrend$Tnd,vuong_SIOP_SIOPtrend$critval),digits=2))
write.csv(output.table,file=paste(savefolder,"/LLR_tests_new.csv",sep=""),na="",row.names = FALSE)
rm(list=setdiff(ls(),keepvars.alltests))

###########################################
###########################################
###########################################



###########################################
# Tables 1 and 5: Confusion matrices for different models (level, siop, miop)

print("Construct confusion matrix")
source(paste0(wd,"/translate_rating.R"))
keepvars.review <- c(keepvars.alltests)

mod <- "baseline"
out <- loadToEnv(paste(wd,"/",mod,".pool.RData",sep=""))#[c("res","res_polr","resfund","res.ziop","res.ziop.fund")]

res <- out[["res"]]
resziop <- out[["resziop"]]
res_oproblevel_agg <- out[["res_oproblevel_agg"]]
do_levelmods <- c(is.null(res_oproblevel_agg))
Yoprob <- out$Yoprob
Yprob <- out$Yprob

# Estimate level model if this has not been done already
if (any(do_levelmods)){
  Ylevel.first <- out$Ylevel
  X <- out$X
  Z <- out$Z
  kX <- dim(X)[2]
  kZ <- dim(Z)[2]
  posc.low.oprob <- out$posc.low.oprob
  posc.high.oprob <- out$posc.high.oprob
  N.const.oprob <- 2*length(out$Nobs)
  Nobs <- out[["Nobs"]]
  corr <- out[["corr"]]
  
  Nobs_cum <- cumsum(c(0,out$Nobs))
  
  # Construct level series from last rating of the month (Ylevel has first rating of the month)
  # Adapt to have different constants for different agencies
  pos.Moody <- (Nobs_cum[1]+1):Nobs_cum[2]
  pos.SP <- (Nobs_cum[2]+1):Nobs_cum[3]
  pos.Fitch <- (Nobs_cum[3]+1):Nobs_cum[4]
  Ylevel.Moody <- out$data.orig[out$pos[pos.Moody],"Rating.last.Moody"]
  Ylevel.SP <- out$data.orig[out$pos[pos.SP],"Rating.last.SP"]
  Ylevel.Fitch <- out$data.orig[out$pos[pos.Fitch],"Rating.last.Fitch"]
  
  print(paste("Estimating level model"))

  # Aggregating rating classes
  matchtable <- data.frame(orig = 0:24,new = c(rep(4,6),rep(seq(7,22,3),each=3),24))
  Ylevel.Moody <- matchtable[match(Ylevel.Moody,matchtable[,1]),2]
  Ylevel.SP <- matchtable[match(Ylevel.SP,matchtable[,1]),2]
  Ylevel.Fitch <- matchtable[match(Ylevel.Fitch,matchtable[,1]),2]
  Ylevel.first <- matchtable[match(Ylevel.first,matchtable[,1]),2]
  
  # Setting up variables, accounting for different thresholds by agency
  Ylevel.polr.orig <- c(Ylevel.Moody,Ylevel.SP,Ylevel.Fitch)
  Ylevels.unique <- sort(unique(Ylevel.polr.orig))
  Ylevel.polr.est <- c(Ylevel.Moody,24+Ylevel.SP,48+Ylevel.Fitch)
  N.const.Moody <- length(unique(Ylevel.Moody))
  Ylevel.count.Moody <- aggregate(data.frame(count = Ylevel.Moody), by = list(value = Ylevel.Moody), length)
  posc.low.Moody <- match(Ylevel.Moody,sort(unique(Ylevel.Moody)))
  posc.high.Moody <- posc.low.Moody+1
  N.const.SP <- length(unique(Ylevel.SP))
  Ylevel.count.SP <- aggregate(data.frame(count = Ylevel.SP), by = list(value = Ylevel.SP), length)
  posc.low.SP <- match(Ylevel.SP,sort(unique(Ylevel.SP))) + N.const.Moody - 1
  posc.high.SP <- posc.low.SP+1
  N.const.Fitch <- length(unique(Ylevel.Fitch))
  Ylevel.count.Fitch <- aggregate(data.frame(count = Ylevel.Fitch), by = list(value = Ylevel.Fitch), length)
  posc.low.Fitch <- match(Ylevel.Fitch,sort(unique(Ylevel.Fitch))) + N.const.Moody + N.const.SP - 2
  posc.high.Fitch <- posc.low.Fitch+1
  
  N.const.polr <- sum(c(N.const.Moody-1,N.const.SP-1,N.const.Fitch-1))
  posc.high.Moody[posc.high.Moody==max(posc.high.Moody)] <- N.const.polr+2
  posc.high.SP[posc.high.SP==max(posc.high.SP)] <- N.const.polr+2
  posc.low.SP[posc.low.SP==min(posc.low.SP)] <- 1
  posc.low.Fitch[posc.low.Fitch==min(posc.low.Fitch)] <- 1
  posc.low.polr <- c(posc.low.Moody,posc.low.SP,posc.low.Fitch)
  posc.high.polr <- c(posc.high.Moody,posc.high.SP,posc.high.Fitch)
  
  coeff_start <- setNames(c(rep(0,dim(Z)[2]),
                            qnorm(cumsum(Ylevel.count.Moody$count[1:(N.const.Moody-1)])/sum(Ylevel.count.Moody$count)),
                            qnorm(cumsum(Ylevel.count.SP$count[1:(N.const.SP-1)])/sum(Ylevel.count.SP$count)),
                            qnorm(cumsum(Ylevel.count.Fitch$count[1:(N.const.Fitch-1)])/sum(Ylevel.count.Fitch$count))),
                          c(colnames(Z),
                            paste0(setdiff(sort(unique(Ylevel.Moody)),24),".Moody"),
                            paste0(setdiff(sort(unique(Ylevel.SP)),24),".SP"),
                            paste0(setdiff(sort(unique(Ylevel.Fitch)),24),".Fitch")))
  colnames.p = list(sort(unique(Ylevel.Moody)),sort(unique(Ylevel.SP)),sort(unique(Ylevel.Fitch)))
  
  # Estimate ordered probit for rating levels
  coeff_polrlevel = optimx(coeff_start, fn=LL_oprob_bounded, gr=grad_oprob_bounded, hess=hess_oprob_bounded,
                           X=Z, posc.low=posc.low.polr, posc.high=posc.high.polr, N.const=N.const.polr,
                           method="BFGS",control=list(kkt=FALSE,follow.on=TRUE,trace=0,maxit=50000,starttests=FALSE))
  
  
  # Transform predictions of rating changes into predictions of rating levels
  res_oproblevel_agg <- LL_oprob_bounded(as.matrix(coeff_polrlevel[1:length(coeff_start)]),X=Z, posc.low=posc.low.polr, posc.high=posc.high.polr, N.const = N.const.polr, onlyLL = FALSE, Nobs = Nobs, colnames.p = colnames.p)
  fitted.values.changes <- matrix(0,length(Yoprob),3,dimnames=list(NULL,c("-1","0","1")))
  for (k in 1:length(Ylevels.unique)){
    posk <- Ylevel.first==Ylevels.unique[k]    # upgrades/downgrades relative to ratings at first day of the month
    if(k>1){fitted.values.changes[posk,1] <- rowSums(as.matrix(res_oproblevel_agg$fitted.values[posk,1:(k-1)]))}
    fitted.values.changes[posk,2] <- res_oproblevel_agg$fitted.values[posk,k]
    if(k<length(Ylevels.unique)){fitted.values.changes[posk,3] <- rowSums(as.matrix(res_oproblevel_agg$fitted.values[posk,(k+1):(length(Ylevels.unique))]))}
  }
  p_case.changes <- rep(0,length(Yoprob))
  for (k in unique(Yoprob)){
    p_case.changes[Yoprob==k] <- fitted.values.changes[Yoprob==k,as.character(k)]
  }
  
  res_oproblevel_agg$LL.changes <- LL_oprob_bounded.probs(fitted.values.changes,Yoprob,Ylevel.first)
  res_oproblevel_agg$fitted.values.changes <- fitted.values.changes
  res_oproblevel_agg$p_case.changes <- p_case.changes
  res_oproblevel_agg$Ylevel.first <- Ylevel.first
  res_oproblevel_agg$Ylevel.predict <- Ylevel.polr.orig
  
  resave(res_oproblevel_agg,file=paste(wd,"/",mod,".pool.RData",sep=""))
}

resnames <- c("res_oproblevel_agg","res","resziop")
matnames <- paste0("confusionmatrix_",c("level_agg","siop","ziop"),".csv")

# Table 1: Confusion matrix for level model
for (r in 1){
  restemp <- get(resnames[r])
  # if (r %in% c(1)){
    pred <- translate_rating(as.numeric(colnames(restemp$fitted.values))[max.col(restemp$fitted.values)])
    classnames <- translate_rating(sort(as.numeric(colnames(restemp$fitted.values)),decreasing=TRUE))
    Ypred <- translate_rating(restemp$Ylevel.predict)
    pos.change <- !(restemp$Ylevel.predict - restemp$Ylevel.first)==0
  # }
  nclasses <- length(classnames)
  cm <- matrix("",nclasses,nclasses,dimnames=list(classnames,classnames))
  obs.change.total <- 0
  for (cm.row in 1:nclasses){
    pos.act <- Ypred == classnames[cm.row]
    for (cm.col in 1:nclasses){
      pos.pred <- pred == classnames[cm.col]
      
      obs.all <- sum(pos.act & pos.pred)
      if ((r %in% c(1:2)) & (cm.col != cm.row)){
        obs.change <- sum(pos.act & pos.pred & pos.change)   # how many observations should have predicted a change but didn't?
      }else{
        obs.change <- 0
      }
      
      if (obs.change==0){
        cm[cm.row,cm.col] <- obs.all
      }else{
        cm[cm.row,cm.col] <- paste0(obs.all," (",obs.change,")")
      }
      obs.change.total <- obs.change.total + obs.change
    }
  }
  print(obs.change.total)
  write.csv(cm,matnames[r])
}

# Table 5: Elements of confusion matrix for siop/miop model
for (r in 3:4){
  restemp <- get(resnames[r])
  fv <- restemp$fitted.values[,c(2,1,3,4)]
  colnames(fv) <- c("00","-1","01","1")
  classnames <- c("00","-1","01","1")
  Ypred <- as.character(Yoprob)
  Ypred[Yprob==0 & Yoprob==0] <- "00"
  Ypred[Yprob==1 & Yoprob==0] <- "01"
  
  nclasses <- length(classnames)
  cm <- matrix(0,nclasses,nclasses+1,dimnames=list(classnames,c("count",classnames)))
  obs.change.total <- 0
  for (cm.row in 1:nclasses){
    pos.act <- Ypred == classnames[cm.row]
    cm[cm.row,1] = sum(pos.act)
    for (cm.col in 1:nclasses){
      cm[cm.row,cm.col+1] <- mean(fv[pos.act,classnames[cm.col]])
    }
  }
  write.csv(cm,matnames[r])
}

rm(list=setdiff(ls(),keepvars.alltests))


