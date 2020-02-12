library("ghyp")
library("Rcpp")


a=get(load("data/a.rda"))
f=get(load("data/f.rda"))


b=function(r0,N){
  r0Nalpha=r0*N^.42
  return(exp(f(r0Nalpha))*r0^1.6)
}

theta=function(r0,N,nu){
  myEstimate=a(r0)-b(r0,N)*nu^(-1.5)
  return(myEstimate)
}

estimateSNR=function(x,numPerm=1000,nu=NA){
  if(is.na(nu))
    nu_fixed=FALSE
  else
    nu_fixed=TRUE
  
  if(length(x)<1){
    return(NULL)
  }
  x=x[is.finite(x)]
  N=length(x)
  if(N<1){
    return(NULL)
  }
  x=as.numeric(x)
  if(is.na(nu) || nu<0 ){
    myfit=tryCatch(fit.tuv(x,silent=TRUE,nu=6),error=function(e){return(NA)})

    if(!class(myfit)=="mle.ghyp"){
      nu=100000  #fit.tuv fails on Gaussian data. TODO: use a test to find it out
    }else{
      nu=coef(myfit)$nu
    }
  }

  R0bar=computeR0bar(x,numPerm = numPerm)   #computes R_+ - R_- over numPerm permutations

  if(abs(R0bar)<1){
    R0bar=0
    return(list(nu=nu,SNR=0,R0bar=0,N=length(x)))
  }

  N=length(x)   # note: this is the length of x[is.finite(x)], which removes NAs and Inf
  SNR_raw=theta(abs(R0bar)/N,N,nu)                #splines have been calibrated for positive SNRs, hence abs(R0bar)

  SNR_raw=SNR_raw * (sign(SNR_raw)>0)                   #if(spline(abs(x)))<0, then set SNR to 0 (this may occur when SNR is very small)
  SNR=SNR_raw*sign(R0bar)                               #then apply sign of R0bar

  if(nu_fixed){
    SNR=SNR*(1-(exp(3.7726-6.2661*log(nu))-exp(-0.36538*nu-1.58686)))*(1-0.009248)
  }else{
    SNR=SNR*((1-(exp(4.4635-6.9240*log(nu))-exp(-0.18441*nu-3.18909)))*(1+exp(-log(nu)*3.2+0.9)-0.009248))
  }

  return(list(SNR=SNR,nu=nu,R0bar=R0bar,N=length(x)))
}

.onUnload <- function (libpath) {
  library.dynam.unload("sharpeRratio", libpath)
}
