#'@title Solve Poisson mean problem by a splitting method
#'@param x data vector
#'@param s scaling vecto
#'@param sigma2 prior variance
#'@param ebnm_params a list of `ebnm` parameters
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
#'pois_mean_split(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'@import ebnm
#'@export
#'

pois_mean_split = function(x,s=NULL,
                           sigma2 = NULL,
                           tol=1e-5,
                           maxiter=1e3,
                           ebnm_params=NULL,
                           optim_method ='L-BFGS-B'){
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
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  # const in objective function
  const = sum((x-1)*log(s)) - sum(lfactorial(x))

  b_pm = rep(0,n)
  #b_pv = rep(1/n,n)
  mu_pm = log(1+x)
  mu_pv = rep(1/n,n)
  if(is.null(sigma2)){
    sigma2 = var(log(x+1))
  }
  for (iter in 1:maxiter) {

    # # VGA
    # for(i in 1:n){
    #   temp = pois_mean_GG1(x[i],s[i],b_pm[i],sigma2,optim_method,mu_pm[i],mu_pv[i])
    #   mu_pm[i] = temp$m
    #   mu_pv[i] = temp$v
    # }
    opt = vga_optimize(c(mu_pm,mu_pv),x,s,b_pm,sigma2)
    mu_pm = opt$m
    mu_pv = opt$v

    # opt = optim(c(mu_pm,log(mu_pv)),
    #             fn = pois_mean_GG_opt_obj,
    #             gr = pois_mean_GG_opt_obj_gradient,
    #             x=x,
    #             s=s,
    #             beta=b_pm,
    #             sigma2=sigma2,
    #             n=n,
    #             #const=const,
    #             method = optim_method)
    # mu_pm = opt$par[1:n]
    # mu_pv = exp(opt$par[(n+1):(2*n)])

    # EBNM
    res = ebnm(mu_pm,sqrt(sigma2),
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    b_pm = res$posterior$mean
    b_pv = res$posterior$sd^2
    H = res$log_likelihood + n*(log(2*pi*sigma2)/2)+sum((mu_pm^2-2*mu_pm*b_pm+b_pm^2+b_pv)/sigma2/2)

    # Update sigma2
    sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)

    # ELBO
    obj[iter+1] = sum(x*mu_pm-s*exp(mu_pm+mu_pv/2)) +const - n/2*log(2*pi*sigma2) - sum(mu_pm^2 + mu_pv + b_pm^2 + b_pv - 2*mu_pm*b_pm)/2/sigma2 + H + sum(log(2*pi*mu_pv))/2 - n/2
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }

  }
  return(list(posterior = list(mean_log = mu_pm,
                               mean_b = b_pm,
                               #osteriorVar_log_mean = mu_pv,
                               #posteriorVar_latent_b = b_pv,
                               mean = exp(mu_pm + mu_pv/2)),
              fitted_g = list(sigma2=sigma2,g_b = res$fitted_g),
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = res))

  # return(list(posteriorMean = mu_pm,
  #             posteriorVar = mu_pv,
  #             sigma2=sigma2,
  #             ebnm_res=res,
  #             obj_value=obj))

}


