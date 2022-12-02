#'
#' #'@title optimize vga for posterior mean and var.
#' #'@param init_val initial value for (m,log(v))
#' #'@param x,s data and scale factor
#' #'@param beta,sigma2 prior mean and variance
#' vga_optimize = function(init_val,x,s,beta,sigma2,method='lbfgs'){
#'   if(method!="lbfgs"){
#'     stop('only lbfgs is implemented')
#'   }
#'   n = length(x)
#'   opt = try(optim(init_val,
#'                   fn = pois_mean_GG_opt_obj,
#'                   gr = pois_mean_GG_opt_obj_gradient,
#'                   x=x,
#'                   s=s,
#'                   beta=beta,
#'                   sigma2=sigma2,
#'                   n=n,
#'                   #const=const,
#'                   method = 'L-BFGS-B'),silent = TRUE)
#'
#'   if(class(opt)=='try-error'){
#'     warning(paste(opt[1],'; returning the initialization values'))
#'     return(list(m=init_val[1:n],v=exp(init_val[(n+1):(2*n)])))
#'   }else{
#'     return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
#'   }
#' }
#'
#' #'@title calculate VGA Poisson objective function, variance = exp(.)
#' pois_mean_GG_opt_obj = function(theta,x,s,beta,sigma2,n){
#'   #n = length(x)
#'   m = theta[1:n]
#'   v = theta[(n+1):(2*n)]
#'   val = -sum(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2)
#'   return(val)
#' }
#'
#' #'@title calculate VGA Poisson objective function gradient, variance = exp(.)
#' pois_mean_GG_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
#'   #n = length(x)
#'   m = theta[1:n]
#'   v = theta[(n+1):(2*n)]
#'   g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
#'   g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
#'   return(c(g1,g2))
#' }
#'
#' #'@title optimize vga for posterior mean then set v = v(m)
#' vga_optimize_m = function(init_val,x,s,beta,sigma2,method='lbfgs'){
#'   if(method!="lbfgs"){
#'     stop('only lbfgs is implemented')
#'   }
#'   n = length(x)
#'   opt = try(optim(init_val,
#'                   fn = pois_mean_GG_opt_obj_m,
#'                   gr = pois_mean_GG_opt_obj_m_gradient,
#'                   x=x,
#'                   s=s,
#'                   beta=beta,
#'                   sigma2=sigma2,
#'                   method = 'L-BFGS-B'),silent = TRUE)
#'
#'   if(class(opt)=='try-error'){
#'     warning(paste(opt[1],'; returning the initialization values'))
#'     return(list(m=init_val,v=v_m(init_val,x,beta,sigma2)))
#'   }else{
#'     return(list(m=opt$par,v=v_m(opt$par,x,beta,sigma2),opt=opt))
#'   }
#' }
#'
#' v_m = function(m,x,beta,sigma2){
#'   sigma2/(sigma2*x-m+beta + 1)
#' }
#'
#' v_dm = function(m,x,beta,sigma2){
#'   sigma2/(sigma2*x-m+beta + 1)^2
#' }
#' v_dm2 = function(m,x,beta,sigma2){
#'   (2*sigma2)/(sigma2*x-m+beta + 1)^3
#' }
#'
#' #'@title calculate VGA Poisson objective function of posterior mean
#' pois_mean_GG_opt_obj_m = function(m,x,s,beta,sigma2){
#'   v = v_m(m,x,beta,sigma2)
#'   val = -sum(x*m-s*exp(m+v/2)-(m^2+v-2*m*beta)/2/sigma2+log(v)/2)
#'   return(val)
#' }
#'
#' #'@title calculate gradient of VGA Poisson objective function of posterior mean
#' pois_mean_GG_opt_obj_m_gradient = function(m,x,s,beta,sigma2){
#'   v = v_m(m,x,beta,sigma2)
#'   v_g = v_dm(m,x,beta,sigma2)
#'   val = -(x-s*exp(m+v/2)*(1+v_g/2)-m/sigma2 - v_g/2/sigma2 + 1/sigma2 + v_g/v/2)
#'   return(val)
#' }
#' #'@title calculate hessian of VGA Poisson objective function of posterior mean
#' pois_mean_GG_opt_obj_m_hess = function(m,x,s,beta,sigma2){
#'   v = v_m(m,x,beta,sigma2)
#'   v_h = v_dm2(m,x,beta,sigma2)
#'   v_g = v_dm(m,x,beta,sigma2)
#'   val = -(-s*exp(m+v/2)*(1+v_g/2)-s*exp(m+v/2)*(v_h/2)-1/sigma2 - v_h/2/sigma2 + v_h/v/2)
#'   return(val)
#' }


#' #'@title calculate VGA Poisson objective function
#' #'@param theta (m,v)
#' pois_mean_GG_opt_obj0 = function(theta,x,s,beta,sigma2,n){
#'   #n = length(x)
#'   m = theta[1:n]
#'   v = theta[(n+1):(2*n)]
#'   val = -sum(x*m-s*exp(m+v/2)-(m^2+v-2*m*beta)/2/sigma2+log(v)/2)
#'   return(val)
#' }
#'
#' #'@title calculate VGA Poisson objective function gradient
#' #'@param theta (m,v)
#' pois_mean_GG_opt_obj_gradient0 = function(theta,x,s,beta,sigma2,n){
#'   #n = length(x)
#'   m = theta[1:n]
#'   v = theta[(n+1):(2*n)]
#'   g1 = -(x-s*exp(m+v/2)-m/sigma2+beta/sigma2)
#'   g2 = -(-1/2*s*exp(m+v/2) - 1/2/sigma2 + 1/2/v)
#'   return(c(g1,g2))
#' }
