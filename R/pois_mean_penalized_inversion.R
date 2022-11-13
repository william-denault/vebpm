#'@title Solve poisson penalized problem via inversion method
#'@description Formulate Poisson mean problem as likelihood + penalty
#'@param x data vector
#'@param w prior weights
#'@param prior_mean prior mean
#'@param mixsd prior sd
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posterior:}{posteriorMean_log_mean is the posterior mean of mu, posteriorMean_mean is the posterior of exp(mu)}
#'  \item{fitted_g:}{estimated prior}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'pois_mean_penalized_inversion(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim \sum_k w_k N(\beta,\sigma_k^2).}
#'@export
#'
pois_mean_penalized_inversion = function(x,
                                         w = NULL,
                                         prior_mean = NULL,
                                         mixsd=NULL,
                                         point_mass = TRUE,
                                         optim_method = 'L-BFGS-B',
                                         tol=1e-8,
                                         maxiter = 1000,
                                         verbose=FALSE,
                                         theta_init = NULL){
  ##########
  # S_inv is too slow
  # takes almost all time
  n = length(x)
  # init prior mean,
  if(is.null(prior_mean)){
    beta = log(sum(x)/n)
  }else{
    beta = prior_mean
  }
  if(is.null(mixsd)){
    ## how to choose grid in this case?
    #mixsd = ebnm:::default_smn_scale(log(x+1),sqrt(1/(x+1)),mode=beta)
    s = 1
    #mixsd = ashr:::autoselect.mixsd(data=list(x = log(0.1/s+x/s),s = sqrt(1/(0.1/s+x/s)),lik=list(name='normal')),sqrt(2),mode=0,grange=c(-Inf,Inf),mixcompdist = 'normal')
    mixsd = select_mixsd(x,1)
    #mixsd = c(1e-10,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
    #mixsd = c(0,1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
  }
  if(point_mass){
    if(mixsd[1]!=0){
      mixsd = c(0,mixsd)
    }
  }
  K = length(mixsd)



  if(is.null(w)){
    w = rep(1/K,K)
  }

  const = lfactorial(x)

  t_start = Sys.time()
  if(is.null(theta_init)){
    theta_init = log(x+1)
  }
  fit = optim(c(theta_init,w,beta),f_obj,f_obj_grad,
              method = optim_method,y=x,grid=mixsd,const=const,
              control=list(trace=verbose,maxit=maxiter,factr=tol/.Machine$double.eps))
  t_end = Sys.time()

  m = fit$par[1:n]
  w_hat = softmax(fit$par[(n+1):(n+K)])
  mu_hat = fit$par[n+K+1]
  s_hat = sqrt(exp(-m))
  z_hat = S_inv(m,s_hat,w_hat,mu_hat,mixsd)
  return(list(posterior = list(mean_log = m,
                               var_log = PV(z_hat,s_hat,w_hat,mu_hat,mixsd),
                               mean = S_exp(z_hat,s_hat,w_hat,mu_hat,mixsd)),
              fitted_g = list(weight = w_hat,mean = mu_hat,sd = mixsd),
              elbo=-fit$value,
              fit =list(z=z_hat,s=s_hat,optim_fit = fit),
              run_time = difftime(t_end,t_start,units='secs')))

  #return(list(posteriorMean = fit$par[1:n], w = softmax(fit$par[(n+1):(n+K)]), priorMean = fit$par[n+K+1], optim_fit = fit,sigma2k=sigma2k))


}

#'@title objective function for inversion method
#'@param params (theta,w,mu)
f_obj = function(params,y,grid,const){
  n = length(y)
  K = length(grid)
  theta = params[1:n]
  a = params[(n+1):(n+K)]
  w = softmax(a)
  mu = params[n+K+1]
  s = sqrt(exp(-theta))
  z = S_inv(theta,s,w,mu,grid)
  val = sum(-y*theta+exp(theta)+const-l_nm(z,s,w,mu,grid)-(theta-z)^2/2/s^2-log(2*pi*s^2)/2)
  #print(theta)
  #print(l_nm(z,s,w,mu,grid))
  #z = S_inv(theta,s,w,mu,grid,z_range)
  #print(val)
  return(val)
}

f_obj_grad=function(params,y,grid,const){
  n = length(y)
  K = length(grid)
  theta = params[1:n]
  a = params[(n+1):(n+K)]
  w = softmax(a)
  mu = params[n+K+1]
  s = sqrt(exp(-theta))
  z = S_inv(theta,s,w,mu,grid)
  #z = S_inv(theta,s,w,mu,grid,z_range)
  grad_theta = exp(theta)-y-l_nm_d1_theta(z,theta,s,w,mu,grid) - (2*s^2*(theta-z)*(1-z_d1_theta(z,theta,s,w,mu,grid))-(-exp(-theta))*(theta-z)^2)/2/s^4 - (-exp(-theta))/2/s^2
  grad_g = colSums(-l_nm_d1_g(z,theta,s,a,mu,grid) - 2*(z-theta)*z_d1_g(z,theta,s,a,mu,grid)/2/s^2)
  #print(z_d1_theta(z,theta,s,w,mu,grid))
  #print(grad_g)
  return(c(grad_theta,grad_g))
}


# f_obj_known_g = function(theta,y,w,mu,grid,z_range){
#   s = sqrt(exp(-theta))
#   z = S_inv(theta,s,w,mu,grid,z_range)
#   return(sum(-y*theta+exp(theta)-l_nm(z,s,w,mu,grid)-(theta-z)^2/2/s^2-log(2*pi*s^2)/2))
# }
#
# f_obj_grad_known_g = function(theta,y,w,mu,grid,z_range){
#   s=sqrt(exp(-theta))
#   z = S_inv(theta,s,w,mu,grid,z_range)
#   exp(theta)-y-l_nm_d1_theta(z,theta,s,w,mu,grid) - (2*s^2*(theta-z)*(1-z_d1_theta(z,theta,s,w,mu,grid))-(-exp(-theta))*(theta-z)^2)/2/s^4 - (-exp(-theta))/2/s^2
# }




