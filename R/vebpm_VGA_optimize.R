
#'@title optimze vga for posterior mean and var
#'@importFrom nloptr lbfgs
vga_optimize = function(init_val,x,s,beta,sigma2,method='lbfgs'){
  if(method!="lbfgs"){
    stop('only lbfgs is implemented')
  }
  n = length(x)
  opt = try(optim(init_val,
              fn = pois_mean_GG_opt_obj,
              gr = pois_mean_GG_opt_obj_gradient,
              x=x,
              s=s,
              beta=beta,
              sigma2=sigma2,
              n=n,
              #const=const,
              method = 'L-BFGS-B'),silent = TRUE)
  if(class(opt)=='try-error'){
    opt = optim(c(rep(0,n),rep(1,n)),
                    fn = pois_mean_GG_opt_obj,
                    gr = pois_mean_GG_opt_obj_gradient,
                    x=x,
                    s=s,
                    beta=beta,
                    sigma2=sigma2,
                    n=n,
                    #const=const,
                    method = 'L-BFGS-B',
                control = list(factr = 1e5))
  }
  return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
  # if(class(opt)=='try-error'){
  #   opt = optim(init_val,
  #               fn = pois_mean_GG_opt_obj0,
  #               gr = pois_mean_GG_opt_obj_gradient0,
  #               x=x,
  #               s=s,
  #               beta=beta,
  #               sigma2=sigma2,
  #               n=n,
  #               #const=const,
  #               method = 'L-BFGS-B',
  #               lower=c(rep(-Inf,n),rep(1e-10,n)))
  #   return(list(m=opt$par[1:n],v=(opt$par[(n+1):(2*n)]),opt=opt))
  # }else{
  #   return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
  # }


}

#'@title calculate VGA Poisson objective function, variance = exp(.)
pois_mean_GG_opt_obj = function(theta,x,s,beta,sigma2,n){
  #n = length(x)
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  val = -sum(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2)
  return(val)
}

#'@title calculate VGA Poisson objective function gradient, variance = exp(.)
pois_mean_GG_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
  #n = length(x)
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}

#'@title calculate VGA Poisson objective function
#'@param theta (m,v)
pois_mean_GG_opt_obj0 = function(theta,x,s,beta,sigma2,n){
  #n = length(x)
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  val = -sum(x*m-s*exp(m+v/2)-(m^2+v-2*m*beta)/2/sigma2+log(v)/2)
  return(val)
}

#'@title calculate VGA Poisson objective function gradient
#'@param theta (m,v)
pois_mean_GG_opt_obj_gradient0 = function(theta,x,s,beta,sigma2,n){
  #n = length(x)
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+v/2)-m/sigma2+beta/sigma2)
  g2 = -(-1/2*s*exp(m+v/2) - 1/2/sigma2 + 1/2/v)
  return(c(g1,g2))
}
