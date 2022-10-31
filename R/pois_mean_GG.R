#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param prior_mean prior mean
#'@param prior_var prior variance
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posterior:}{posteriorMean_latent is the posterior mean of mu, posteriorMean_mean is the posterior of exp(mu)}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'pois_mean_GG(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}
#'@export


pois_mean_GG = function(x,
                        s = NULL,
                        prior_mean = NULL,
                        prior_var=NULL,
                        optim_method = 'L-BFGS-B',
                        maxiter = 1000,
                        tol = 1e-5){

  # init the posterior mean and variance?
  n = length(x)
  m = log(x+0.1)
  v = rep(1/n,n)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  #
  if(is.null(prior_mean) | is.null(prior_var)){

    if(is.null(prior_mean)){
      est_beta = TRUE
    }else{
      est_beta = FALSE
      beta = prior_mean
    }
    if(is.null(prior_var)){
      est_sigma2=TRUE
    }else{
      est_sigma2 = FALSE
      sigma2=prior_var
    }

    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    for(iter in 1:maxiter){
      if(est_beta){
        beta = mean(m)
      }
      if(est_sigma2){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
      }
      # for(i in 1:n){
      #   temp = pois_mean_GG1(x[i],s[i],beta,sigma2,optim_method,m[i],v[i])
      #   m[i] = temp$m
      #   v[i] = temp$v
      # }
      opt = optim(c(m,log(v)),
                  fn = pois_mean_GG_opt_obj,
                  gr = pois_mean_GG_opt_obj_gradient,
                  x=x,
                  s=s,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  method = optim_method)
      m = opt$par[1:n]
      v = exp(opt$par[(n+1):(2*n)])
      obj[iter+1] = pois_mean_GG_obj(x,s,beta,sigma2,m,v)
      if((obj[iter+1] - obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }

  }else{
    beta = prior_mean
    sigma2 = prior_var
    # for(i in 1:n){
    #   temp = pois_mean_GG1(x[i],s[i],prior_mean,prior_var,optim_method,m[i],v[i])
    #   m[i] = temp$m
    #   v[i] = temp$v
    # }
    opt = optim(c(m,log(v)),
                fn = pois_mean_GG_opt_obj,
                gr = pois_mean_GG_opt_obj_gradient,
                x=x,
                s=s,
                beta=beta,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    m = opt$par[1:n]
    v = exp(opt$par[(n+1):(2*n)])
    obj = pois_mean_GG_obj(x,s,prior_mean,prior_var,m,v)

  }

  return(list(posterior = list(posteriorMean_log_mean = m,
                               posteriorVar_log_mean = v,
                               posteriorMean_mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              obj_value=obj))

  #return(list(posteriorMean=m,priorMean=beta,priorVar=sigma2,posteriorVar=v,obj_value=obj))

}

#'@title calculate VGA Poisson objective function
#'@export
pois_mean_GG_opt_obj = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2))
}

#'@title calculate VGA Poisson objective function gradient
#'@export
pois_mean_GG_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}


pois_mean_GG_obj = function(x,s,beta,sigma2,m,v){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2))
}

#' #'@param x a data point
#' #'@param beta prior mean
#' #'@param sigma2 prior variance
#' #'@param optim_method optimization method in `optim` function
#' pois_mean_GG1 = function(x,s,
#'                             beta,
#'                             sigma2,
#'                             optim_method = 'BFGS',
#'                             m_init  = NULL,
#'                             v_init = NULL){
#'   # init m, v
#'   if(is.null(m_init)){
#'     m = 0
#'   }else{
#'     m = m_init
#'   }
#'
#'   if(is.null(v_init)){
#'     v = 1
#'   }else{
#'     v = v_init
#'   }
#'
#'   opt = optim(c(m,log(v)),
#'               fn = pois_mean_GG1_obj,
#'               gr = pois_mean_GG1_obj_gradient,
#'               x=x,
#'               s=s,
#'               beta=beta,
#'               sigma2=sigma2,
#'               method = optim_method)
#'
#'   return(list(m=opt$par[1],v=exp(opt$par[2]),obj=-opt$value))
#' }
#'
#' #'calculate objective function
#' pois_mean_GG1_obj = function(theta,x,s,beta,sigma2){
#'   return(-(x*theta[1]-s*exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2/sigma2+log(exp(theta[2]))/2))
#' }
#'
#' #'calculate gradient vector
#' pois_mean_GG1_obj_gradient = function(theta,x,s,beta,sigma2){
#'   g1 = -(x-s*exp(theta[1]+exp(theta[2])/2)-theta[1]/sigma2+beta/sigma2)
#'   g2 = -(-exp(theta[2])/2*s*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2/sigma2 + 1/2)
#'   return(c(g1,g2))
#' }



