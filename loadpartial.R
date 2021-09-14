loadpartial <- function(filename,varnames){
  temp <- load(filename)
  if (length(setdiff(varnames,temp))>0){
    warning("not all requested elements are in the given file")
    varnames <- intersect(varnames,temp)
  }
  out <- lapply(varnames,FUN=function(x){return(x=get(x))})
  names(out) <- varnames
  return(out)
}