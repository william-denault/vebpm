#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Poisson mean problem, log(1+exp(x)) link function.
#'@param x data vector
#'@param maxiter max number of iterations
#'@param ebnm_params a list of `ebnm` parameters
#'@param tol tolerance for stopping the updates
#'@param kapa default to be 0.25 + 0.17*max(x)
#'@return a list of
#'  \item{posterior:}{posteriorMean_latent is the posterior mean of mu, posteriorMean_mean is for Poisson mean}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'  \item{fit:}{fitted object from `ebnm`}
#'@examples
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'pois_mean_log1exp(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\log(1+\exp(\mu_i))),}
#'\deqn{\mu_i\sim g(\cdot).}
#'@import ebnm
#'@export
#'

pois_mean_log1exp = function(x,ebnm_params = NULL,tol=1e-5,maxiter=1e3,kapa = NULL){
  n = length(x)
  if(is.null(kapa)){
    kapa = 0.25+0.17*max(x)
  }

  pseudo_s = sqrt(1/kapa)
  const = sum(lfactorial(x))

  mu_tilde = x
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
  for(iter in 1:maxiter){
    # update g,q by performing ebnm
    pseudo_x = mu_tilde-nll_d1(x,mu_tilde)/kapa
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
    H = res$log_likelihood + sum(log(2*pi*pseudo_s^2)/2)+sum((pseudo_x^2-2*pseudo_x*m+m^2+v)/pseudo_s^2/2)

    # update mu_tilde
    mu_tilde = m

    # calc objective function
    obj[iter+1] = -sum(nll(x,mu_tilde)+nll_d1(x,mu_tilde)*(m-mu_tilde)+kapa/2*(m^2+v+mu_tilde^2-2*m*mu_tilde))+H - const
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }

  }
  return(list(posterior = list(mean_latent = m,
                               var_latent = v,
                               mean = log(1+exp(m)),
                               mean_log = log(log(1+exp(m)))),
              fitted_g = res$fitted_g,
              obj_value=obj,
              fit = res))

  #return(list(posteriorMean=m,posteriorVar=v,obj_value=obj,ebnm_res=res,kappa=kapa))
}

#'@title calc log(1+exp(x))
#'@description from package `qgam`
log1pexp = function (x){
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), right = TRUE,
                   include.lowest = TRUE)
  kk <- which(indx == 1)
  if (length(kk)) {
    x[kk] <- exp(x[kk])
  }
  kk <- which(indx == 2)
  if (length(kk)) {
    x[kk] <- log1p(exp(x[kk]))
  }
  kk <- which(indx == 3)
  if (length(kk)) {
    x[kk] <- x[kk] + exp(-x[kk])
  }
  return(x)
}

# log1exp = function(x){
#   # log(1+exp(40)) == 40 is TRUE
#   # log(1+exp(-40)) == 0 is TRUE
#   if(x<40){
#     if(x>-40){
#       return(log(1+exp(x)))
#     }else{
#       return(0)
#     }
#   }else{
#     return(x)
#   }
# }

nll = function(x,mu){
  return(log1pexp(mu)-x*log(log1pexp(mu)))
}

nll_d1 = function(x,mu){
  n = (log1pexp(mu)-x)
  d = (1+exp(-mu))*log1pexp(mu)
  return(n/d)
}


