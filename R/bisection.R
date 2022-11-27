#'@title bisection for root finding. Vectorized.
#'@description modified from https://stat.ethz.ch/pipermail/r-help/2012-November/340295.html
#'@param f If auto_adjust_interval is TRUE, f must be an increasing function.
#'@param auto_adjust_interval if signs of f(lower) and f(upper) are not opposite, search for new interval. If TRUE ,f must be an increasing function

bisection= function(f, lower, upper, ...,
                    maxiter =1000,
                    tol = 1e-8,
                    search_step = 1,
                    auto_adjust_interval = TRUE) {

  if(auto_adjust_interval){
    # make sure all f(lower)<0 and all f(upper)>0
    if(length(search_step)==1){
      search_step = rep(search_step,length(lower))
    }
    l = which(f(lower,...) > 0)
    u = which(f(upper,...) < 0)
    counter = 0
    while((sum(l)+sum(u))!=0){
      if(sum(l)!=0){
        lower[l] = lower[l] - search_step[l]
        fl = f(lower,...)
        if(any(is.nan(fl))){
          stop('NaN in f(lower)')
        }
        l = which(fl > 0)
        #print(l)
      }
      if(sum(u)!=0){
        upper[u] = upper[u] + search_step[u]
        fu = f(upper,...)
        if(any(is.nan(fu))){
          stop('NaN in f(upper)')
        }
        u = which(fu < 0)
        #print(u)
      }
      counter = counter + 1
      if(counter>maxiter){
        stop('Cannot find searching interval, consider increasing maxiter')
      }
    }
  }



  fl=f(lower, ...)
  fu=f(upper, ...)

  if(any(fl*fu>0)){
    print(which(fl*fu > 0))
    stop('Need to adjust lower and upper such that f(lower) and f(upper) have opposite signs')
  }

  for(n in 1:maxiter) {
    mid=(lower+upper)/2
    fmid=f(mid, ...)
    if(all(abs(fmid) < tol)){
      break
    }
    #samesign = ((fmid<0)&(fl<0))|((fmid>=0)&(fl>=0))
    samesign = ((fmid<0)&(fl<0))
    lower = ifelse(samesign, mid, lower)
    fl = ifelse(samesign, fmid, fl)
    upper = ifelse(!samesign, mid, upper)
    fu = ifelse(!samesign, fmid, fu)
    if(n == maxiter){
      warning('Not converged yet.')
    }
  }
  return(mid)
}
