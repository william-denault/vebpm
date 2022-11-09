#'@title Solve Negative Binomial mean problem by Jaakkola Jordan bound
#'@param x data vector
#'@param r NB binomial parameter r
#'@param ebnm_params a list of `ebnm` parameters
#'@param est_r whether estimate r. For Poisson data, set it to False.
#'@param update_r_every if estimate r, how often to update it
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posterior:}{posteriorMean_latent is the posterior mean of mu, posteriorMean_mean is for Poisson mean, posteriorMean_log_mean is for Poisson log mean}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'  \item{fit:}{fitted object from `ebnm`}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'nb_mean_lower_bound(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\log(1+\exp(\mu_i))),}
#'\deqn{\mu_i\sim g(\cdot).}
#'
#'@import ebnm
#'@export


nb_mean_lower_bound = function(x,
                               r = 1e2,
                               ebnm_params = NULL,
                               tol=1e-5,
                               maxiter=1000,
                               est_r = FALSE,
                               update_r_every=1,
                               m_init = NULL,
                               v_init = NULL){
  n = length(x)

  obj = rep(0,maxiter+1)
  obj[1] = -Inf
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }

  if(is.null(m_init)){
    m = log(x/r+1)
  }
  if(is.null(v_init)){
    v = rep(1/n,n)
  }

  if(is.null(r)){
    r = round(mean(x)*2)
    est_r = TRUE
  }
  r_trace = r
  for(iter in 1:maxiter){

    # update xi
    xi = sqrt(m^2+v)

    # update g,q by performing ebnm
    pseudo_s = sqrt(2*xi/(x+r)/tanh(xi/2))
    pseudo_x = (x-r)*pseudo_s^2/2
    res = ebnm(pseudo_x,pseudo_s,
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    m = res$posterior$mean
    v = res$posterior$sd^2
    # get Elog(g/q)
    H_mu = res$log_likelihood + sum(log(2*pi*pseudo_s^2)/2)+sum((pseudo_x^2-2*pseudo_x*m+m^2+v)/pseudo_s^2/2)

    if(est_r){
      if(iter%%update_r_every==0){
        r = update_r_lb(x,m,v,xi,r)
        r_trace = c(r_trace,r)
      }
    }
    # calc objective function
    obj[iter+1] = neg_binom_mean_lb_obj(x,m,v,r,xi,H_mu)
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }
  }
  if(iter==maxiter){
    message('Not converged - Increase number of iterations.')
  }
  p = 1/(1+exp(-m))




  return(list(posterior = list(mean_latent = m,
                               var_latent = v,
                               mean = r*S_exp(pseudo_x,pseudo_s,res$fitted_g$pi,res$fitted_g$mean[1],res$fitted_g[[3]]),
                               mean_log = log(r)+m),
              fitted_g = res$fitted_g,
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = res))

  # return(list(posteriorMean=m,
  #             posteriorVar=v,
  #             r=r,
  #             obj_value=obj,
  #             ebnm_res=res,
  #             p=1/(1+exp(-m)),
  #             poisson_mean_est = r*exp(m),
  #             poisson_log_mean_est = log(r) + m,
  #             r_trace=r_trace))
}

neg_binom_mean_lb_obj = function(x,m,v,r,xi,H_mu){
  n = length(x)
  #xi = sqrt(m^2+v)
  temp = sum((x-r)/2*m + (x+r)*(-xi/2-log(1+exp(-xi))-1/4/xi*tanh(xi/2)*(m^2+v-xi^2)) + lgamma(x+r)) - n*lgamma(r)
  return(temp)
}

update_r_lb = function(x,m,v,xi,r_init){

  delta = sum(-xi/2-log(1+exp(-xi))-1/4/xi*tanh(xi/2)*(m^2+v-xi^2)-m/2)

  res = optim(log(r_init),
              Fr_lb,
              Fr_lb_d1,
              x=x,
              m=m,
              v=v,
              delta=delta,
              method='BFGS')
  return(exp(res$par))
}

# r is log(r) in these functions.
Fr_lb = function(r,x,m,v,delta){
  n = length(x)
  temp = (exp(r)*delta+sum(lgamma(x+exp(r)))-n*lgamma(exp(r)))
  return(-temp)
}
Fr_lb_d1 = function(r,x,m,v,delta){
  n = length(x)
  temp = (delta+sum(digamma(x+exp(r)))-n*digamma(exp(r)))
  return(-temp)
}



