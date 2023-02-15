#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param prior_mean prior mean
#'@param prior_var prior variance
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posterior:}{mean_log/var_log is the posterior mean/var of mu, mean is the posterior of exp(mu)}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'ebpm_normal(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}
#'@export


ebpm_normal = function(x,
                        s = NULL,
                       g_init = NULL,
                       fix_g = FALSE,
                       q_init = NULL,
                        maxiter = 1000,
                        tol = 1e-5,
                        vga_tol=1e-5,
                        conv_type='elbo'){

  # init the posterior mean and variance?
  n = length(x)

  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  if(is.null(q_init)){
    m = log(x/s+1)
    v = rep(1/n,n)
  }else{
    if(is.null(q_init$m_init)){
      m = log(1+x/s)
    }else{
      m = q_init$m_init
    }
    if(is.null(q_init$v_init)){
      v = rep(1/n,n)
    }else{
      v = q_init$v_init
    }
  }

  const = sum((x-1)*log(s)) - sum(lfactorial(x))
  #
  t_start = Sys.time()

  if(is.null(g_init)){
    prior_mean = NULL
    prior_var = NULL
  }else{
    prior_mean = g_init$mean
    prior_var = g_init$var
  }

  if(length(fix_g)==1){
    est_prior_mean = !fix_g
    est_prior_var = !fix_g
  }else if(length(fix_g)==2){
    est_prior_mean = !fix_g[1]
    est_prior_var = !fix_g[2]
  }else{
    stop('fix_g can be either length 1 or 2')
  }


  if(est_prior_mean | est_prior_var){

    if(is.null(prior_mean)){
      est_prior_mean = TRUE
      beta = mean(m)
    }else{
      beta = prior_mean
    }
    if(is.null(prior_var)){
      est_prior_var=TRUE
      sigma2 = mean(m^2+v-2*m*beta+beta^2)
    }else{
      sigma2=prior_var
    }

    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    for(iter in 1:maxiter){
      sigma2_old = sigma2
      m = vga_pois_solver(m,x,s,beta,sigma2,tol=vga_tol)
      v =  m$v
      m = m$m

      if(est_prior_mean){
        beta = mean(m)
      }
      if(est_prior_var){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
      }

      if(conv_type=='elbo'){
        obj[iter+1] = ebpm_normal_obj(x,s,beta,sigma2,m,v,const)
        if((obj[iter+1] - obj[iter])/n <tol){
          obj = obj[1:(iter+1)]
          if((obj[iter+1]-obj[iter])<0){
            warning('An iteration decreases ELBO. This is likely due to numerical issues.')
          }
          break
        }
      }
      if(conv_type=='sigma2abs'){
        obj[iter+1] = abs(sigma2-sigma2_old)
        if(obj[iter+1] <tol){
          obj = obj[1:(iter+1)]
          break
        }
      }

    }

  }else{
    beta = prior_mean
    sigma2 = prior_var
    m = vga_pois_solver(m,x,s,beta,sigma2,tol=vga_tol)
    v = m$v
    m = m$m
    obj = ebpm_normal_obj(x,s,prior_mean,prior_var,m,v,const)

  }
  t_end = Sys.time()

  return(list(posterior = list(mean_log = m,
                               var_log = v,
                               mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              elbo=ebpm_normal_obj(x,s,beta,sigma2,m,v,const),
              obj_trace = obj,
              run_time = difftime(t_end,t_start,units='secs')))

}


ebpm_normal_obj = function(x,s,beta,sigma2,m,v,const){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2)+const)
}

#' #'@param x a data point
#' #'@param beta prior mean
#' #'@param sigma2 prior variance
#' #'@param optim_method optimization method in `optim` function
#' ebpm_normal1 = function(x,s,
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
#'               fn = ebpm_normal1_obj,
#'               gr = ebpm_normal1_obj_gradient,
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
#' ebpm_normal1_obj = function(theta,x,s,beta,sigma2){
#'   return(-(x*theta[1]-s*exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2/sigma2+log(exp(theta[2]))/2))
#' }
#'
#' #'calculate gradient vector
#' ebpm_normal1_obj_gradient = function(theta,x,s,beta,sigma2){
#'   g1 = -(x-s*exp(theta[1]+exp(theta[2])/2)-theta[1]/sigma2+beta/sigma2)
#'   g2 = -(-exp(theta[2])/2*s*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2/sigma2 + 1/2)
#'   return(c(g1,g2))
#' }



