#'@title Optimize vga poisson, using bisection for v
#'@param x,s data and scale factor
#'@param beta,sigma2 prior mean and variance
#'@export
vga_pois_solver = function(init_val,x,s,beta,sigma2,maxiter=1000,tol=1e-8){

  n = length(x)
  if(length(sigma2)==1){
    upper = rep(sigma2,n)
  }else if(length(sigma2)==n){
    upper = sigma2
  }else{
    stop('check length of sigma2')
  }
  # bisection
  v = try(bisection(h_v,
                lower = rep(0,n),upper = upper,
                x=x,s=s,beta=beta,sigma2=sigma2,
                auto_adjust_interval = FALSE,
                maxiter=maxiter,tol=tol),silent = TRUE)
  # if not working, try optim
  if(class(v)=='try-error'){
    return(vga_optimize_m(init_val,x,s,beta,sigma2))
  }else{
    m = sigma2*x + beta + 1 - sigma2/v
    return(list(m=m,v=v))
  }

}

h_v = function(v,x,s,beta,sigma2){
  val = (-s*exp(sigma2*x+beta+1-sigma2/v + v/2) - 1/sigma2 + 1/v)
  return(-val)
}
