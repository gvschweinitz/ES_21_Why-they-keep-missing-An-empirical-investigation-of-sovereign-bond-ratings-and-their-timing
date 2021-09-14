# SCRIPT for estimating SIOP for rating data
# Runs the pooling tests necessary to decide which variables should enter X, Z or both.
# Output: Appendix table A.3
#
# last changed: 09/14/2021
#
# author: mei, gsz


rm(list=ls())
wd <- getwd()

load(paste(wd,"/rating_data_calc.RData",sep=""))

varlist <- c(ls(),"varlist")
data.orig <- data.m
rm(list=varlist)
wd <- getwd()
source(paste(wd,"/siop.R",sep=""))
source(paste(wd,"/grad_oprob_bounded.R",sep=""))
source(paste(wd,"/hess_oprob_bounded.R",sep=""))
source(paste(wd,"/LL_oprob_bounded.R",sep=""))
source(paste(wd,"/LL_oprob_bounded.probs.R",sep=""))
source(paste(wd,"/LL_ziopc_bounded.R",sep=""))
source(paste(wd,"/loadpartial.R",sep=""))
source(paste(wd,"/fun_collection.R",sep=""))
# source(paste(wd,"/fun_collection_derivatives.R",sep=""))
source(paste(wd,"/select_data.R",sep=""))
source(paste(wd,"/select_mean_sd.R",sep=""))
source(paste(wd,"/make_result_table.R",sep=""))

packages <- c("optimx","MASS","Matrix","pbivnorm","numDeriv")
for (p in packages){
  check = require(p,character.only=TRUE)
  if (!check) {
    install.packages(p)
    require(p,character.only=TRUE)
  }  
}

###############################
# Data identifiers; rating agencies and agency-based variables
crossid <- "Country"
dateid <- "Date"     
agencies <- c("Moody","SP","Fitch")
cols.agency <- c("Rating.first","Rating.last","Rating.changes.pos","Rating.changes.neg",
                 "Outlook.first","Outlook.last","Announcement","Timesince")
cols.other.bin <- c("Rating.changes.pos","Rating.changes.neg","Announcement")
cols.other.meandiff <- c("Rating.first")

###############################
# Set "default" parameters
lag_setting = 12        # number of lags for lag modification
alpha_setting = 0.99    # windsorizing
alpha_2 <- c(1-alpha_setting,alpha_setting)
realcol <- "cpi_IMF"
b_fundchange <- 1
conf.ag <- "none"
conf.pool <- "none"

###############################
# Separation of variables
colfund <- c("gnipc","ip","inf","reer","yield","debt","fiscbal","current","reserves","corrupt","ka_open")
col.candidates <- list(rating=c("rating","rating.sq"),UpDown=c("Up12","Down12"),UpDown.other=c("Up12.other","Down12.other"),
                       UpDownAll=c("UpAll12.tot","DownAll12.tot"),
                       years=c("years"),
                       ratingdiff=c("rating.diff.pos","rating.diff.neg"),N.other="N.other",gnipc=c("gnipc"),
                       ip="ip",inf="inf",reer="reer",yield="yield",debt="debt",fiscbal="fiscbal",current="current",reserves="reserves",
                       corrupt="corrupt",ka_open="ka_open",default="default",changefund="changefund")
col.full <- unique(unlist(col.candidates))

load.pool <- c("Yoprob","Yprob","Ylevel","posc.low","posc.high","X","Z","Nobs","res","pos","adj_mat")
b.constant.pool <- FALSE
changefund.placement <- "jointXZ"

###############################
# list of elements to be saved from every specification
save.elems <- c("agencies","crossid","dateid","lag_setting","alpha_setting","alpha_2","realcol","interactcol","interactcol2","savename","b_fundchange","colfund",
                "additions","addX","addZ","addXZ","baseX","baseZ","baseXZ","pos","corr","col.ind","adj_mat",
                "data.orig","Yoprob","Yprob","Ylevel","posc.low","posc.high","X","Z","Nobs","res","ptm","restable","changefund.placement",
                "p.poolingtest")

###############################
# list of specifications to be run in the for-loop below
spec.list <- list(list(savename="checkpool",additions=c("baseadd"),corr=FALSE,b_fundchange=1,interactcol=NULL))

runmodels <- 1

###############################
# Definition of variable transformations for baseline scenarios
# Transformations are performed by functions in "fun_collection.R"
# They are based on the syntax "OLD_NAME:FUNCTION:NEW_NAME". Additional qualifiers (normalize; factor; demean) can be added using "!QUALIFIER"

# Set baseline scenario used for all estimations
#   variables in probit (X) AND ordered probit (Z)
baseXZ = c("Rating.first:keep:rating!normalize=TRUE!factor=24!demean=FALSE",
           "rating:square:rating.sq",
           "Rating.changes.pos:calc_cumsum:Up12!normalize=FALSE",
           "Rating.changes.neg:calc_cumsum:Down12!normalize=FALSE",
           "Timesince:keep:years!normalize=TRUE!factor=12!demean=FALSE",
           # "years:square:years.sq",
           "Rating.first.diff.other.pos:keep:rating.diff.pos!normalize=TRUE!factor=24",
           "Rating.first.diff.other.neg:keep:rating.diff.neg!normalize=TRUE!factor=24",
           "Rating.changes.pos.total:calc_cumsum:UpAll12.tot!byCountry=FALSE!normalize=FALSE",
           "Rating.changes.neg.total:calc_cumsum:DownAll12.tot!byCountry=FALSE!normalize=FALSE",
           "Rating.changes.pos.other:calc_cumsum:Up12.other!normalize=FALSE",
           "Rating.changes.neg.other:calc_cumsum:Down12.other!normalize=FALSE",
           "N.other::N.other",
           "GDP_capita:norm_country:gnipc",
           # "gnipc:square:gnipc.sq",
           "ip.comb:calc_gr:ip",
           "Reserves:calc_windgr:reserves",
           "cpi_IMF:calc_windsor:inf",
           "reer.comb:calc_windgr:reer!alpha=alpha_2",
           "ryields.rEMBI:calc_windsor:yield",
           "debt.comb:keep:debt!normalize=TRUE!factor=100",
           "fiscal.balance.comb:keep:fiscbal!normalize=TRUE!factor=100",
           "ext.balance.comb:keep:current!normalize=TRUE!factor=100",
           "corruption_TI::corrupt",
           "ka_open::ka_open")
baseaddXZ = c("DummyDefPast::default")

#   variables only in probit
baseX = NULL
baseaddX = NULL

#   variables only in ordered probit
baseZ = NULL
baseaddZ =NULL



##################################
##################################
##################################
# Start of the main program
#   For-loop over model specifications (this case: two different pooling specifications per variable)
for (k in runmodels){
  print(spec.list[[k]])
  
  additions <- spec.list[[k]]$additions
  corr <- spec.list[[k]]$corr
  interactcol <- spec.list[[k]]$interactcol
  interactcol2 <- spec.list[[k]]$interactcol2
  b_fundchange <- spec.list[[k]]$b_fundchange
  
  estimations <- c(agencies,"full.pool",col.candidates)
  table.poolingtest <- matrix(NA,length(c(col.candidates)),12,
                              dimnames=list(names(col.candidates),
                                            c("LL.onepool","LL.nonpooled","N.onepool","N.nonpooled","stat.onepool","p.onepool",
                                              "LL.oneind","LL.fullpooled","N.oneind","N.fullpooled","stat.oneind","p.oneind")))
  
  # Run estimation over (1) individual agencies; (2) full pooled model; (3) single indicators pooled/non-pooled
  for  (pos.ag in 1:length(c(estimations,col.candidates))){
    if (pos.ag>length(estimations)){
      ag <- col.candidates[[pos.ag-length(estimations)]]
      onepool <- FALSE
      tabpos <- pos.ag-length(estimations)
    }else{
      ag <- estimations[[pos.ag]]
      onepool <- TRUE
      tabpos <- pos.ag-4
    }
    
    savename <- paste(spec.list[[k]]$savename,ag,sep=".")
    
    if (ag%in%agencies){
      # Construct data for agency-specific estimations
      
      data <- data.orig
      data[,cols.agency] <- data[,paste(cols.agency,ag,sep=".")]
      data[,"Outlook.pos"] <- 1*(data[,"Outlook.first"]>0)
      data[,"Outlook.neg"] <- 1*(data[,"Outlook.first"]<0)
      
      ag.other <- setdiff(agencies,ag)
      data[,"N.other"] <- rowSums(!is.na(data[,paste("Rating.first",ag.other,sep=".")]))
      for (cname in cols.other.bin){
        data[,cname] <- (data[,cname]>0)*1
        data[,paste(cname,"other",sep=".")] <- (rowSums(data[,paste(cname,ag.other,sep=".")],na.rm=TRUE)>0)*1
        data[,paste(cname,"total",sep=".")] <- (rowSums(data[,paste(cname,agencies,sep=".")],na.rm=TRUE)>0)*1
      }
      for (cname in cols.other.meandiff){
        temp <- rowMeans(data[,paste(cname,ag.other,sep=".")],na.rm=TRUE)-data[,cname]
        temp[is.na(temp)] <- 0
        data[,paste(cname,"diff.other",sep=".")] <-  temp
        data[,paste(cname,"diff.other.pos",sep=".")] <- temp
        data[temp<0,paste(cname,"diff.other.pos",sep=".")] <- 0
        data[,paste(cname,"diff.other.neg",sep=".")] <- temp
        data[temp>0,paste(cname,"diff.other.neg",sep=".")] <- 0
      }
      addX <- NULL
      addZ <- NULL
      addXZ <- NULL
      
      # Define dependent variable
      Yoprob <- sign(data[,"Rating.last"]-data[,"Rating.first"])
      
      data[data[,"Rating.first"]<=2 & !is.na(data[,"Rating.first"]),"Rating.first"] <- 2
      
      Yprob <- data[,"Announcement"]
      Ylevel <- data[,"Rating.first"]
      # add all required specifications together
      if (length(additions)>0) {
        for (add in additions)  {
          addX = c(addX,get(paste(add,"X",sep="")))
          addZ = c(addZ,get(paste(add,"Z",sep="")))
          addXZ = c(addXZ,get(paste(add,"XZ",sep="")))
        }
      }
      
      jointXZ <- select_data(c(baseXZ,addXZ), data)
      onlyZ <- select_data(c(baseZ,addZ),cbind(data,jointXZ))
      onlyX <- select_data(c(baseX,addX),cbind(data,jointXZ,onlyZ))
      
      # Cut data for main model and run
      pos <- which(rowSums(is.na(cbind(Yoprob,Yprob,jointXZ,onlyX,onlyZ,Ylevel)))==0)
      
      Nobs <- length(pos)
      jointXZ <- jointXZ[pos,]
      onlyX <- onlyX[pos,]
      onlyZ <- onlyZ[pos,]
      Yoprob = as.numeric(Yoprob[pos])
      Yprob = Yprob[pos]
      Ylevel <- Ylevel[pos]
      # boundlow <- Ylevel==min(Ylevel)
      # boundhigh <- Ylevel==max(Ylevel)
      posc.low <- Yoprob+2
      posc.low[(Ylevel==min(Ylevel)) & (Yoprob<1)] <- 1
      posc.high <- Yoprob+3
      posc.high[(Ylevel==max(Ylevel)) & (Yoprob>-1)] <- 4
      
      if (b_fundchange) {
        changefund <- as.matrix(calc_fundchange(cbind(jointXZ,onlyX,onlyZ),colfund,data[pos,crossid],Yprob,normalize=TRUE))
        colnames(changefund) <- "changefund"
        if (!is.null(interactcol)){
          changefund <- cbind(changefund,changefund * jointXZ[,interactcol])
          colnames(changefund) <- c("changefund",paste("changefund.X.",interactcol,sep=""))
        }
        if (changefund.placement=="jointXZ"){
          jointXZ <- cbind(jointXZ,changefund)
        }else if (changefund.placement=="onlyX"){
          onlyX <- cbind(onlyX,changefund)
        }else if (changefund.placement=="onlyZ"){
          onlyZ <- cbind(onlyZ,changefund)
        }
      }
      X = cbind(jointXZ,onlyX)
      Z = cbind(jointXZ,onlyZ)
      
      # Calculating adj_mat
      jointXZ.orig <- select_mean_sd(c(baseXZ,addXZ), data)
      onlyZ.orig <- select_mean_sd(c(baseZ,addZ),cbind(data,jointXZ.orig))
      onlyX.orig <- select_mean_sd(c(baseX,addX),cbind(data,jointXZ.orig,onlyZ.orig))
      jointXZ.orig <- jointXZ.orig[pos,]
      onlyZ.orig <- onlyZ.orig[pos,]
      onlyX.orig <- onlyX.orig[pos,]
      
      FULL <- cbind(X,Z[,setdiff(colnames(Z),colnames(X))])
      FULL.orig <- cbind(jointXZ.orig,onlyX.orig,onlyZ.orig)[,intersect(colnames(FULL),colnames(cbind(jointXZ.orig,onlyX.orig,onlyZ.orig)))]
      if (b_fundchange){
        FULL.orig <- cbind(FULL.orig,FULL[,setdiff(colnames(FULL),colnames(FULL.orig))])
        colnames(FULL.orig)[colnames(FULL.orig)==""] <- setdiff(colnames(FULL),colnames(cbind(jointXZ.orig,onlyX.orig,onlyZ.orig)))
      }
      adj_mat <- NULL
      for (col.adj in colnames(FULL)){
        check <- FALSE
        while (!check){
          pos_ms <- sample(1:dim(X)[1],2)
          check <- !(diff(FULL[pos_ms,col.adj])==0)
        }
        y <- FULL[pos_ms,col.adj]
        x <- FULL.orig[pos_ms,col.adj]
        sd <- (x[1]-x[2])/(y[1]-y[2])
        mean <- x[1]-y[1]*sd
        adj_mat <- rbind(adj_mat,cbind(sd,mean))
      }
      rownames(adj_mat) <- colnames(FULL)
      adj_mat <- adj_mat[!grepl(".sq",rownames(adj_mat)) & !grepl(".X.",rownames(adj_mat)),]
      
      conf <- conf.ag
    }else{
      print("Estimating pooled model: drawing data from individual files")
      if (ag == "full.pool"){
        print("pooled: all variables except constants")
        col.ind <- NULL
      }else{
        if (onepool){
          col.ind <- setdiff(col.full,ag)
          print(paste("pooled:",ag))
        }else{
          col.ind <- ag
          print(paste("not pooled:",ag))
        }
        
      }
      
      # Construct data from agency-specific estimations
      var.Moody <- loadpartial(paste(spec.list[[k]]$savename,"Moody.RData",sep="."),load.pool)
      var.SP <- loadpartial(paste(spec.list[[k]]$savename,"SP.RData",sep="."),load.pool)
      var.Fitch <- loadpartial(paste(spec.list[[k]]$savename,"Fitch.RData",sep="."),load.pool)
      
      # prior results for poolability tests
      if (any(is.null(var.Moody$res),is.null(var.SP$res),is.null(var.Fitch$res))){
        LL_ind <- -Inf
        Ncoeff_ind <- 0
      }else{
        LL_ind <- var.Moody$res$LL + var.SP$res$LL + var.Fitch$res$LL
        Ncoeff_ind <- length(c(var.Moody$res$opt_coeff,var.SP$res$opt_coeff,var.Fitch$res$opt_coeff))
      }
      
      
      Yoprob <- c(var.Moody$Yoprob,var.SP$Yoprob,var.Fitch$Yoprob)
      Yprob <- c(var.Moody$Yprob,var.SP$Yprob,var.Fitch$Yprob)
      Ylevel <- c(var.Moody$Ylevel,var.SP$Ylevel,var.Fitch$Ylevel)
      
      colnames.add = function(cnames,add,col.ind){
        colsq <- grep(".sq",cnames)
        colX <- grep(".X.",cnames)
        if (length(c(colsq,colX))==0){return(paste(cnames,add,sep="."))}
        
        cnames.new <- cnames
        cnames.new[-unique(c(colsq,colX))] <- paste(cnames.new[-unique(c(colsq,colX))],add,sep=".")
        cnames.new[setdiff(colsq,colX)] <- paste(unlist(strsplit(cnames.new[setdiff(colsq,colX)],".sq")),add,"sq",sep=".")
        for (k in setdiff(colX,colsq)){
          parts <- strsplit(cnames.new[k],".X.")[[1]]
          parts[parts %in% col.ind] <- paste(parts[parts %in% col.ind],add,sep=".")
          cnames.new[k] <- paste(parts,collapse=".X.")
        }
        for (k in intersect(colX,colsq)){
          parts <- strsplit(cnames.new[k],".sq.X.")[[1]]
          parts[parts %in% col.ind] <- paste(parts[parts %in% col.ind],add,sep=".")
          cnames.new[k] <- paste(parts,collapse=".sq.X.")
        }
        return(cnames.new)
      }
      
      if (!is.null(col.ind)){
        colind.X <- intersect(colnames(var.Moody$X),col.ind)
        colind.Z <- intersect(colnames(var.Moody$Z),col.ind)
      }else{
        colind.X <- NULL
        colind.Z <- NULL
      }
      
      colpool.X <- setdiff(colnames(var.Moody$X),colind.X)
      colpool.Z <- setdiff(colnames(var.Moody$Z),colind.Z)
      X.pool <- rbind(as.matrix(var.Moody$X[,colpool.X]),as.matrix(var.SP$X[,colpool.X]),as.matrix(var.Fitch$X[,colpool.X]))
      Z.pool <- rbind(as.matrix(var.Moody$Z[,colpool.Z]),as.matrix(var.SP$Z[,colpool.Z]),as.matrix(var.Fitch$Z[,colpool.Z]))
      
      if (length(colind.X)>0){
        X.ind <- as.matrix(bdiag(var.Moody$X[,colind.X],var.SP$X[,colind.X],var.Fitch$X[,colind.X]))
        colnames(X.ind) <- c(colnames.add(colind.X,"Moody",col.ind),colnames.add(colind.X,"SP",col.ind),colnames.add(colind.X,"Fitch",col.ind))
      }else{
        X.ind <- NULL
      }
      if (length(colind.Z)>0){
        Z.ind <- as.matrix(bdiag(var.Moody$Z[,colind.Z],var.SP$Z[,colind.Z],var.Fitch$Z[,colind.Z]))
        colnames(Z.ind) <- c(colnames.add(colind.Z,"Moody",col.ind),colnames.add(colind.Z,"SP",col.ind),colnames.add(colind.Z,"Fitch",col.ind))
      }else{
        Z.ind <- NULL
      }
      X <- cbind(X.pool,X.ind)
      Z <- cbind(Z.pool,Z.ind)
      
      if (max(abs(var.Moody$adj_mat-var.SP$adj_mat),abs(var.Moody$adj_mat-var.Fitch$adj_mat))<10^-8){
        adj_mat <- var.Moody$adj_mat
      }else{
        colind.adjmat <- rownames(var.Moody$adj_mat)[unique(unlist(apply(as.matrix(col.ind),1,FUN=function(x){return(grep(x,rownames(var.Moody$adj_mat)))})))]
        colpool.adjmat <- setdiff(rownames(var.Moody$adj_mat),colind.adjmat)
        adj_mat <- rbind(var.Moody$adj_mat[colpool.adjmat,],
                         var.Moody$adj_mat[colind.adjmat,],
                         var.SP$adj_mat[colind.adjmat,],
                         var.Fitch$adj_mat[colind.adjmat,])
        rownames(adj_mat) <- c(colpool.adjmat,
                               colnames.add(colind.adjmat,"Moody",col.ind),colnames.add(colind.adjmat,"SP",col.ind),colnames.add(colind.adjmat,"Fitch",col.ind))
      }
      
      pos <- c(var.Moody$pos,var.SP$pos,var.Fitch$pos)
      Nobs <- c(var.Moody$Nobs,var.SP$Nobs,var.Fitch$Nobs)
      
      if (b.constant.pool){
        posc.low <- c(var.Moody$posc.low,var.SP$posc.low,var.Fitch$posc.low)
        posc.high <- c(var.Moody$posc.high,var.SP$posc.high,var.Fitch$posc.high)
      }else{
        #positions of thresholds: if thresholds are not to be pooled, they should be c(-Inf,c_1^Moody,c_2^Moody,c_1^SP,...,Inf)
        posshifter <- function(x,shift1,shift2){
          # position of thresholds in posc.low and posc.high need to be shifted:
          # Position of Inf by 4
          # Position of polr-thresholds depending on the agency: Moody=0,SP=2,Fitch=4
          xout <- x
          xout[x>1 & x<4] <- x[x>1 & x<4]+shift1
          xout[x==4] <- x[x==4]+shift2
          return(xout)
        }
        
        posc.low <- c(posshifter(var.Moody$posc.low,0,4),posshifter(var.SP$posc.low,2,4),posshifter(var.Fitch$posc.low,4,4))
        posc.high <- c(posshifter(var.Moody$posc.high,0,4),posshifter(var.SP$posc.high,2,4),posshifter(var.Fitch$posc.high,4,4))
      }
      
      conf <- conf.pool
    }
    
    if (corr){conf="none"}
    
    rm(res)
    ptm <- proc.time()
    res <- tryCatch(
      {
        siop(Yoprob=Yoprob,Yprob=Yprob,Ylevel=Ylevel,X=X,Z=Z,posc.low.oprob=posc.low,posc.high.oprob=posc.high,corr=corr,conf=conf,method=c("BFGS","Nelder-Mead","CG"),b.constant.pool=b.constant.pool,Nobs=Nobs,agencies=agencies)
      },
      error = function(cond){
        message(cond)
        return(NULL)
      }
    )
    ptm <- proc.time()-ptm
    
    if (!is.null(res)){
      print(paste("Likelihood Model",k,ag,": ",res$LL))
      
      if (pos.ag <=4){
        if (ag == "full.pool"){
          LL_full <- res$LL
          Ncoeff_full <- length(res$opt_coeff)
        }
        p.poolingtest <- NA
        save(list=intersect(save.elems,ls()),file=paste(savename,".RData",sep=""))
      }else{
        if (onepool){
          p.poolingtest <- 1-pchisq(2*(LL_ind-res$LL),Ncoeff_ind-length(res$opt_coeff))
          table.poolingtest[tabpos,1:6] <- c(res$LL,LL_ind,length(res$opt_coeff),Ncoeff_ind,2*(LL_ind-res$LL),p.poolingtest)
        }else{
          p.poolingtest <- 1-pchisq(2*(res$LL-LL_full),length(res$opt_coeff)-Ncoeff_full)
          table.poolingtest[tabpos,7:12] <- c(res$LL,LL_full,length(res$opt_coeff),Ncoeff_full,2*(res$LL-LL_full),p.poolingtest)
        }
      }
    }else{
      pause
    }
    rm(res)
  }
}

write.csv(table.poolingtest,file="tab.poolingtest.csv")