#'@title Solve Poisson mean problem by a splitting method
#'@param x data vector
#'@param s scaling vector
#'@param sigma2 prior variance
#'@param est_sigma2 whether fix sigma2 at input or update it in iterations
#'@param ebnm_params a list of `ebnm` parameters
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@param mu_pm_init init value of E_q(mu), default to log(1+x)
#'@return a list of
#'  \item{posterior:}{posterior mean of mu, posterior of exp(mu), etc.}
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
#'\deqn{b_i\sim g(.).}
#'@import ebnm
#'@export
#'

pois_mean_split = function(x,s=NULL,
                           sigma2 = NULL,
                           est_sigma2 = TRUE,
                           tol=1e-5,
                           maxiter=1e3,
                           ebnm_params=NULL,
                           mu_pm_init = NULL){
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

  # #b_pv = rep(1/n,n)
  # if(is.null(mu_pm_init)){
  #   mu_pm_init = log(1/s+x/s)
  # }
  # mu_pm = mu_pm_init
  #
  # #mu_pv = rep(1/n,n)
  # if(is.null(sigma2)){
  #   sigma2 = var(mu_pm)
  #   est_sigma2 = TRUE
  # }
  if(is.null(sigma2)){
    est_sigma2 = TRUE
    init_res = ebpm_normal(x,s,g_init=list(mean=log(sum(x)/sum(s)),var=NULL),fix_g=c(T,F),maxiter=maxiter,tol=tol)
    mu_pm = init_res$posterior$mean_log
    sigma2 = init_res$fitted_g$var
  }else{
    mu_pm = log(1/s+x/s)
  }
  if(!is.null(mu_pm_init)){
    mu_pm = mu_pm_init
  }
  #sigma2 = min(sigma2,var(mu_pm)/4)
  sigma2_trace = c()
  t_start = Sys.time()
  for (iter in 1:maxiter) {

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

    #print(paste('after ebnm, obj=',splitting_obj(x,s,mu_pm,mu_pv,b_pm,b_pv,sigma2,H,const,n)))

    opt = vga_pois_solver(mu_pm,x,s,b_pm,sigma2)
    #opt = vga_optimize(c(mu_pm,log(mu_pv)),x,s,b_pm,sigma2)
    mu_pm = opt$m
    mu_pv = opt$v

    #print(paste('after vga, obj=',splitting_obj(x,s,mu_pm,mu_pv,b_pm,b_pv,sigma2,H,const,n)))

    # Update sigma2
    if(est_sigma2){
      sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)
      sigma2_trace[iter] = sigma2
    }

    #print(paste('after sigma2, obj=',splitting_obj(x,s,mu_pm,mu_pv,b_pm,b_pv,sigma2,H,const,n)))

    # ELBO
    obj[iter+1] = splitting_obj(x,s,mu_pm,mu_pv,b_pm,b_pv,sigma2,H,const,n)
    if((obj[iter+1]-obj[iter])/n<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }

  }
  t_end = Sys.time()
  if(class(res$fitted_g)=="normalmix"){
    mean_exp_b = S_exp(mu_pm,sqrt(sigma2),w=res$fitted_g$pi,mu = res$fitted_g$mean,grid = res$fitted_g$sd)
  }else{
    mean_exp_b = exp(b_pm)
  }
  return(list(posterior = list(mean_log = mu_pm,
                               mean_b = b_pm,
                               mean_exp_b = mean_exp_b,
                               var_log = mu_pv,
                               var_b = b_pv,
                               mean = exp(mu_pm + mu_pv/2)),
              fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace,g_b = res$fitted_g),
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = res,
              run_time = difftime(t_end,t_start,units='secs')))

}

splitting_obj = function(x,s,mu_pm,mu_pv,b_pm,b_pv,sigma2,H,const,n){
  sum(x*mu_pm-s*exp(mu_pm+mu_pv/2)) +const - n/2*log(2*pi*sigma2) - sum(mu_pm^2 + mu_pv + b_pm^2 + b_pv - 2*mu_pm*b_pm)/2/sigma2 + H + sum(log(2*pi*mu_pv))/2 - n/2
}

# calculate the posterior expectation of exp(x) with a mixture of uniform prior.
# S_exp_unif = function(fitted_g){
#   sum(fitted_g$pi * (exp(fitted_g$b)-exp(fitted_g$a)))
# }



