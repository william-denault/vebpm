#'@title Solve Negative Binomial mean problem by polya-gamma augmentation
#'@param x data vector
#'@param r NB binomial parameter r
#'@param ebnm_params a list of ebnm parameters
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
#'n = 100
#'mu = rnorm(n)
#'p = exp(mu)/(1+exp(mu))
#'r = 10
#'x = rnbinom(n,r,1-p)
#'nb_mean_polya_gamma(x,r=r)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\log(1+\exp(\mu_i))),}
#'\deqn{\mu_i\sim g(\cdot).}
#'
#'@import ebnm
#'@export

nb_mean_polya_gamma = function(x,
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
  t_start = Sys.time()
  for(iter in 1:maxiter){

    # update Ew
    Ew = (x+r)/(2*sqrt(m^2+v))*tanh(sqrt(m^2+v)/2)

    # update g,q by performing ebnm
    pseudo_x = (x-r)/2/Ew
    pseudo_s = sqrt(1/Ew)
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
       r = update_r(x,m,v,r)
       r_trace = c(r_trace,r)
      }
    }
    # calc objective function
    obj[iter+1] = neg_binom_mean_pg_obj(x,m,v,r,H_mu)
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }
  }
  # if(iter==maxiter){
  #   message('Not converged - Increase number of iterations.')
  # }

  t_end = Sys.time()
  return(list(posterior = list(mean_latent = m,
                               var_latent = v,
                               mean = r*S_exp(pseudo_x,pseudo_s,res$fitted_g$pi,res$fitted_g$mean[1],res$fitted_g[[3]]),
                               mean_log = log(r)+m),
              fitted_g = res$fitted_g,
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = list(ebnm_res = res,r=r),
              run_time = difftime(t_end,t_start,units='secs')))
  #
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

neg_binom_mean_pg_obj = function(x,m,v,r,H_mu){
  #return(sum(-Ew*(m^2+v)/2+(x-r)*m/2)+H_w+H_mu + sum(lgamma(x+r)-lgamma(r)-(x+r)*log(2)))
  Ew = (x+r)/(2*sqrt(m^2+v))*tanh(sqrt(m^2+v)/2)
  H_w = sum(-(x+r)*log(cosh(sqrt(m^2+v)/2))) + sum((m^2+v)/2*Ew)
  return(sum(-Ew*(m^2+v)/2+(x-r)*m/2)+H_w+H_mu + sum(lgamma(x+r)-lgamma(r)-(x+r)*log(2)))
}

update_r = function(x,m,v,r_init){
  delta = -sum(m)/2-sum(log(cosh(sqrt(m^2+v)/2)))-n*log(2)
  res = optim(log(r_init),
              Fr,
              Fr_d1,
              x=x,
              m=m,
              v=v,
              delta=delta,
              method='BFGS')
  return(exp(res$par))
}

# r is log(r) in these functions.
Fr = function(r,x,m,v,delta){
  n = length(x)
  temp = (exp(r)*delta+sum(lgamma(x+exp(r)))-n*lgamma(exp(r)))
  return(-temp)
}
Fr_d1 = function(r,x,m,v,delta){
  n = length(x)
  temp = (delta+sum(digamma(x+exp(r)))-n*digamma(exp(r)))
  return(-temp)
}


