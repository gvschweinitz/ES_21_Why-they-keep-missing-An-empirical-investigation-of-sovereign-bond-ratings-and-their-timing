# FUN_COLLECTION is a script that initializes auxilliary function for the ZIOP routines
#   All functions use a data.frame (data) and a column identifier (col)
#   Most functions allow for normalization or standardization
#   Many functions allow lagging
# 
# Functions:
# keep          "keeps" a variable from the dataset  (default = standardization)
# interact      creates interaction terms (default = no normalization standardization)
# square        generates square (default = no normalization standardization)
# lag           creates lag (default = standardization)
# calc_gr       creates growth rates ((default = standardization)
# calc_windgr   creates windsorized growth rates (default = standardization)
# calc_cumsum   creates rolling cumlative sums over previous lags (default = 12), for either one country (default) or all other countries (default = normalization)
# calc_windsor  windsorizes the data (default = standardization)
# norm_country  expresses every observation as a fraction of the reference country at the same time
#
# Internal functions (not meant for direct calling)
# normalize     normalization
# standardize   standardization
#
# last changed: 12/18/2015
#
# author: gsz,mei  

calc_fundchange <- function(data,colfund,countryvec,Yprob,normalize=TRUE){
  Xfund <- data[,colfund]
  Xfund[abs(Xfund)<10^-6] <- 10^-6
  T <- dim(Xfund)[1]
  K <- dim(Xfund)[2]
  countries <- unique(countryvec)
  Xbm <- matrix(NA,T,K)
  if (normalize) {normchange <- Xbm}
  for (coun in countries){
    pos <- which(countryvec==coun)
    if (length(pos)>1){
      if (normalize){
        normchange[pos[2:length(pos)],] <- (Xfund[pos[2:length(pos)],] - Xfund[pos[1:(length(pos)-1)],])/Xfund[pos[1:(length(pos)-1)],]
      }
    }
    posannounce <- unique(c(1,which(Yprob[pos]==1)+1,length(pos)+1))
    Xbm.c <- matrix(Xbm[pos,],length(pos),K)
    for (n in 1:(length(posannounce)-1)){
      rows <- posannounce[n]:(posannounce[n+1]-1)
      Xbm.c[rows,] <- matrix(Xfund[pos[posannounce[n]],],length(rows),K,byrow=TRUE)
    }
    Xbm[pos,] <- Xbm.c
  }
  if (normalize) {
    normvec <- colMeans(normchange,na.rm=TRUE)
  }else{
    normvec <- rep(1,K)
  }
  Xchange <- ((Xfund-Xbm)/Xbm) / matrix(normvec,T,K,byrow=TRUE)
  for (k in 1:K){
    a1 <- quantile(Xchange[,k],0.05)
    Xchange[Xchange[,k]<a1,k] <- a1
    a2 <- quantile(Xchange[,k],0.95)
    Xchange[Xchange[,k]>a2,k] <- a2
  }
  changefund <- rowMeans(Xchange^2)
  changefund <- (changefund-mean(changefund))/sd(changefund)
  return(changefund)
  # data$fundchange <- rowSums(Xchange^2)
  # res <- NULL
  # res$Xbm <- Xbm
  # res$normchange <- normchange
  # res$normvec <- normvec
  # res$Xfund <- Xfund
  # res$Xchange <- Xchange
  # res$changefund <- changefund
  # return(res)
}


keep <- function(data,col,standardize=TRUE,normalize=FALSE,...){
  if (standardize+normalize==2){
    warning("standardization AND normalization set to true. Data will ONLY be normalized")
    standardize = FALSE
  }
  out <- data[,col]
  if (standardize) {out <- standardize(out)}
  if (normalize) {out <- normalize(out,...)}
  return(out)
}

interact <- function(data,col,col2,standardize=FALSE,normalize=FALSE,factor=NULL,...){
  if (standardize+normalize==2){
    warning("standardization AND normalization set to true. Data will ONLY be standardized")
    normalize = FALSE
  }
  out <- data[,col]*data[,col2]
  if (standardize) {out <- standardize(out)}
  if (normalize) {out <- normalize(out,factor)}
  return(out)
}

square <- function(data, col,standardize=FALSE,normalize=FALSE,factor=NULL,...) {
  if (standardize+normalize==2){
    warning("standardization AND normalization set to true. Data will ONLY be standardized")
    normalize = FALSE
  }
  out <- data[,col]^2
  if (standardize) {out <- standardize(out)}
  if (normalize) {out <- normalize(out,factor)}
  return(out)
}

norm_country <- function(data,col,refcountry = "United.States",...){
  out <- data[,col]
  refpos <- which(data[,crossid]==refcountry)
  refdates <- data[refpos,dateid]
  countries <- unique(data[,crossid])  
  for (coun in countries){
    pos <- which(data[,crossid]==coun)
    posmatch <- pmatch(data[pos,dateid],data[refpos,dateid])  #positions in refpos that match the dates in pos
    if (length(posmatch)>0){
      out[pos[is.na(posmatch)]] <- NA
      pos <- pos[!is.na(posmatch)]
      posmatch <- posmatch[!is.na(posmatch)]
      out[pos] <- data[pos,col]/data[refpos[posmatch],col]
    }
  }
return(out)
}


lag <- function(data,col,lag,standardize=TRUE,normalize=FALSE,factor=NULL,...){
  if (standardize+normalize==2){
    warning("standardization AND normalization set to true. Data will ONLY be standardized")
    normalize = FALSE
  }
  out <- rep(NA,dim(data)[1])
  
  for (coun in unique(data[,crossid])){
    pos <- which(data[,crossid]==coun)
    N <- length(pos)
    x <- data[pos,col]
    if (lag>=0){out[pos[(lag+1):N]] <- x[1:(N-lag)]}
    if (lag<0){out[pos[1:(N-lag)]] <- x[(lag+1):N]}
#     if (length(x)>lag){
#       out[pos] <- c(rep(NA,lag),(exp(diff(log(x),lag=lag))-1)*100)
#     }
#     else{
#       out[pos] <- rep(NA,length(x))
#     }   
  }
  if (!is.null(real)){out <- out-data[,real]}
  if (standardize) {out <- standardize(out)}
  if (normalize) {out <- normalize(out,factor)}
  return(out)
}


calc_gr <- function(data,col,lag=lag_setting,real=NULL,standardize=TRUE,...){
  out <- rep(NA,dim(data)[1])
  
  for (coun in unique(data[,crossid])){
    pos <- which(data[,crossid]==coun)
    x <- data[pos,col]
    if (length(x)>lag){
      out[pos] <- c(rep(NA,lag),(exp(diff(log(x),lag=lag))-1)*100)
    }
    else{
      out[pos] <- rep(NA,length(x))
    }   
  }
  if (!is.null(real)){out <- out-data[,real]}
  if (standardize) {out <- standardize(out)}
  return(out)
}

calc_windgr <- function(data,col,lag=lag_setting,alpha=alpha_setting,real=NULL,...){
  out <- rep(NA,dim(data)[1])
  
  for (coun in unique(data[,crossid])){
    pos <- which(data[,crossid]==coun)
    x <- data[pos,col]
    if (length(x)>lag){
      out[pos] <- c(rep(NA,lag),(exp(diff(log(x),lag=lag))-1)*100)
    }
    else{
      out[pos] <- rep(NA,length(x))
    }   
  }
  if (!is.null(real)){out <- out-data[,real]}
  out <- calc_windsor(out,alpha=alpha,...)
  return(out)
}

calc_cumsum <- function(data,col,lag=lag_setting,current = FALSE, byCountry=TRUE,standardize=FALSE,normalize=TRUE,...){
  if (standardize+normalize==2){
    warning("standardization AND normalization set to true. Data will ONLY be standardized")
    normalize = FALSE
  }
  
  if (current) {add = 0} else {add = 1}
  
  if (lag==0){
    out <- data[,col]
  }else{
    out <- rep(NA,dim(data)[1])
    if (!byCountry){
      dates <- sort(unique(data[,dateid]))
      total <- mat.or.vec(length(dates),3)
      for (k in 1:length(dates)){
        pos <- data[,dateid]==dates[k]
        total[k,] <- c(k,sum(data[pos,col],na.rm=TRUE),sum(!is.na(data[pos,col])))
      }
    }
    
    for (coun in unique(data[,crossid])){
      pos <- which(data[,crossid]==coun)  
      N <- length(pos)
      if (byCountry){
        x <- data[pos,col]
        y <- mat.or.vec(N,1)+1/abs(lag-add+1)
      }else{
        pos_total <- match(data[pos,dateid],dates)
        x <- (total[pos_total,2]-data[pos,col])
        y <- total[pos_total,3]-1
      }
      
      if (N>abs(lag)){
        if (lag>0){
          out[pos[1:lag]] <- NA
          for (k in (lag+1):N){
            if (sum(!is.na(x[(k-lag):(k-add)]))==0){
              out[pos[k]] <- NA
            }else{
              out[pos[k]] <- sum(x[(k-lag):(k-add)],na.rm=TRUE) / sum(y[(k-lag):(k-add)],na.rm=TRUE)
            }
          }
        }else if (lag<0){
          lead <- -lag
          out[pos[(N-lead+1):N]] <- NA
          for (k in (1:(N-lead))){
            if (sum(!is.na(x[(k+add):(k+lead)]))==0){
              out[pos[k]] <- NA
            }else{
              out[pos[k]] <- sum(x[(k+add):(k+lead)],na.rm=TRUE) / sum(y[(k+add):(k+lead)],na.rm=TRUE)
            }
          }
        }
      }else{
        out[pos] <- rep(NA,N)
      } 
    }
  }
  if (standardize) {out <- standardize(out)}
  if (normalize) {out <- normalize(out,...)}
  return(out)
}

calc_windsor <- function(data,col = NULL,alpha=alpha_setting,real=NULL,standardize=TRUE,...){
  cl <- class(data)
  if (cl!="data.frame") {data <- as.data.frame(data)}
  if (is.null(col)) {col = colnames(data)} 
  for (c in col) {
    pos <- which(!is.na(data[,c]))
    temp <- data[pos,c]
    if (length(alpha)==1) {
      temp[temp>quantile(temp,alpha)] <- quantile(temp,alpha)
    } else {
      temp[temp<quantile(temp,alpha[1])] <- quantile(temp,alpha[1])
      temp[temp>quantile(temp,alpha[2])] <- quantile(temp,alpha[2])
    }
    if (standardize) {temp <- standardize(temp)}
    data[pos,c] <- temp
  }
  if (cl=="numeric"){data <- as.matrix(data)}
  data <- as(data,cl)
  if (cl!="numeric"){data <- data[,col]}
  if (!is.null(real)){
    for (c in colnames(data)){
      data[,c] <- data[,c]-data[,real]  
    }
  }
  return(data)
}

normalize <- function(vec,factor=NULL,demean=FALSE,fullinfo=FALSE,...){
  if (is.null(factor)) {std <- max(abs(vec),na.rm=TRUE)
  }else {std = factor}
  vec <- vec/std
  if (demean){
    add <- mean(vec,na.rm=TRUE)
  }else{
    add <- 0
  }
  vec <- vec-add
  
  if (!fullinfo){return(vec)}
  out <- list()
  out$add <- add
  out$vec <- vec
  out$std <- std
}

standardize <- function(vec,fullinfo=FALSE,...) {
  add = mean(vec,na.rm= TRUE)
  std <- sd(vec, na.rm = TRUE)
  vec = (vec - add)/std
  
  if (!fullinfo){return(vec)}
  out <- list()
  out$vec <- vec
  out$add <- add
  out$std <- std
  return(out)
}


colnames.add = function(cnames,add){
  # Function puts "add" between variable name and addendum ".sq" or interaction signifier ".X."
  colsq <- grep(".sq",cnames)
  colX <- grep(".X.",cnames)
  if (length(c(colsq,colX))==0){return(paste(cnames,add,sep="."))}
  
  cnames.new <- cnames
  cnames.new[-unique(c(colsq,colX))] <- paste(cnames.new[-unique(c(colsq,colX))],add,sep=".")
  cnames.new[setdiff(colsq,colX)] <- paste(unlist(strsplit(cnames.new[setdiff(colsq,colX)],".sq")),add,"sq",sep=".")
  for (k in setdiff(colX,colsq)){
    parts <- strsplit(cnames.new[k],".X.")[[1]]
    parts[parts %in% cnames] <- paste(parts[parts %in% cnames],add,sep=".")
    cnames.new[k] <- paste(parts,collapse=".X.")
  }
  for (k in intersect(colX,colsq)){
    parts <- strsplit(cnames.new[k],".sq.X.")[[1]]
    parts[parts %in% cnames] <- paste(parts[parts %in% cnames],add,sep=".")
    cnames.new[k] <- paste(parts,collapse=".sq.X.")
  }
  return(cnames.new)
}
