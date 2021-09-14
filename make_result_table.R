make_result_table <- function(results,ptm,agencies=NULL,p.poolingtest="NA") {
  # MAKE_RESULT_TABLE generates a table from the results estimated by SIOP.R
  #
  # results   result list as produced by SIOP.R
  #
  # last changed: 10/9/2017
  #
  # author: mei,gsz
  
  prob_coeff <- results$prob_coeff
  if(!is.matrix(prob_coeff)){prob_coeff <- matrix(prob_coeff,nrow=1,ncol=length(prob_coeff),dimnames=list(1,names(prob_coeff)))}
  pos.inter <- grep("Intercept",colnames(prob_coeff))
  if (length(pos.inter)>1){
    prob.const <- paste("Constant",agencies,sep=".")
    oprob.const <- c(paste("Thresh -1.0",agencies,sep="."),paste("Thresh 0.1",agencies,sep="."))
  }else{
    prob.const <- "Constant"
    oprob.const <- c("Thresh -1.0","Thresh 0.1")
  }
  colnames(prob_coeff)[pos.inter] <- prob.const
  
  oprob_coeff <- results$oprob_coeff
  if(!is.matrix(oprob_coeff)){oprob_coeff <- matrix(oprob_coeff,nrow=1,ncol=length(oprob_coeff),dimnames=list(1,names(oprob_coeff)))}
  pos.thresh1 <- which(grepl("X.1.0",colnames(oprob_coeff)))
  pos.thresh2 <- which(grepl("X0.1",colnames(oprob_coeff)))
  colnames(oprob_coeff)[c(pos.thresh1,pos.thresh2)] <- oprob.const
  
  pos.p_coeffs <- setdiff(1:length(prob_coeff),pos.inter)
  pos.o_coeffs <- setdiff(1:length(prob_coeff),c(pos.thresh1,pos.thresh2))
  colnames(oprob_coeff)[pos.o_coeffs] <- unlist(strsplit(colnames(oprob_coeff)[pos.o_coeffs],".1",fixed=TRUE))
  p_coeffs <- colnames(prob_coeff)
  
  o_coeffs <- colnames(oprob_coeff)
  if (length(results$opt_coeff)-(length(prob_coeff)+length(oprob_coeff))==1){
    rho <- (exp(2*results$rho)-1)/(exp(2*results$rho)+1)
    o_coeffs <- c(o_coeffs,"corr.")
    oprob_coeff <- t(as.matrix(c(oprob_coeff,rho)))
    colnames(oprob_coeff) <- o_coeffs
  }
  breakpoint1 <- length(prob_coeff)
  breakpoint2 <- length(prob_coeff)+length(oprob_coeff)
  all_coeffs = unique(c(p_coeffs,o_coeffs))
  # print(all_coeffs)
  p.char <- rep("",breakpoint2)
  if (!is.null(results[["p"]])){
    p = results[["p"]]
    p.char[is.na(p)] <- "NA"
    p.char[p<0.1] <- "*"
    p.char[p<0.05] <- "**"
    p.char[p<0.01] <- "***"
  }else{
    p.char <- rep("NA",breakpoint2)
  }

  table = mat.or.vec(length(all_coeffs)*2,4)
  
  pi <- match(p_coeffs,all_coeffs)
  pi.not <- setdiff(1:length(all_coeffs),pi)
  table[pi*2-1,1] <- format(round(as.numeric(prob_coeff),digits=3),nsmall=3)
  if (is.null(results[["sd"]])){
    table[pi*2,1] <- "NA"
  }else{
    for (k in 1:breakpoint1){table[pi[k]*2,1] <- paste("[",format(round(results$sd[k],digits=3),nsmall=3),"]",sep="")}
  }
  table[pi*2-1,2] = p.char[1:breakpoint1]
  table[pi*2,2] = ""
  table[pi.not*2-1,1:2] <- ""
  table[pi.not*2,1:2] <- ""
  
  po <- match(o_coeffs,all_coeffs)
  po.not <- setdiff(1:length(all_coeffs),po)
  table[po*2-1,3] <- format(round(as.numeric(oprob_coeff),digits=3),nsmall=3)
  if (is.null(results[["sd"]])){
    table[po*2,1] <- "NA"
  }else{
    for (k in (breakpoint1+1):breakpoint2){table[po[k-breakpoint1]*2,3] <- paste("[",format(round(results$sd[k],digits=3),nsmall=3),"]",sep="")}
  }
  # table[po*2,3] <- ifelse(is.null(results[["sd"]]),"NA",
  #                         paste("[",format(round(results$sd[(breakpoint1+1):breakpoint2],digits=3),nsmall=3),"]",sep=""))
  table[po*2-1,4] = p.char[(breakpoint1+1):breakpoint2]
  table[po*2,4] = ""
  table[po.not*2-1,3:4] <- ""
  table[po.not*2,3:4] <- ""
  
  nms = rep("",length(all_coeffs)*2)
  nms[(1:length(all_coeffs))*2-1] = all_coeffs
  table <- rbind(table,c(format(round(results$LL,digits = 3),nsmall=3),rep("",3)))
  table <- rbind(table,c(format(round(p.poolingtest,3)),rep("",3)))
  table <- rbind(table,c(as.character(length(Yoprob)),rep("",3)))
  table <- rbind(table,c(format(round(ptm[3],3)),rep("",3)))
  table <- cbind(c(nms,"LL","p.poolingtest","N","time"),table)
  colnames(table) = c("rownames","Coeffs (probit)", "Sig. (probit)", "Coeffs (ordered)", " Sig. (ordered)")
  table <- as.data.frame(table)
  return(table)
}