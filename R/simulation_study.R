
#'@title generate simulation data, gamma
#'@details The DGP is
#'\deqn{x_i\sim Poisson(\lambda_i),}
#'\deqn{\lambda_i\sim w\times Gamma(shape0,rate0) + (1-w)\times Gamma(shape1,rate1).}
#'@export

gen_data_gamma = function(n=1e3,n_simu=100,
                          prior='gamma',
                          shape0 = 2,
                          rate0=2,
                          shape1 = 10,
                          rate1=2,
                          #t_df = 3,
                          #t_delta=2,
                          #quasi_sparse=0.1,
                          #exp_rate = 0.5,
                          w = 0.8,
                          seed=12345){

  set.seed(seed)

  Mean = matrix(nrow=n_simu,ncol=n)
  X = matrix(nrow=n_simu,ncol=n)
  if(prior=='gamma'){
    for(i in 1:n_simu){
      d0 = rgamma(n,shape0,rate0)
      d1 = rgamma(n,shape1,rate1)
      idx = rbinom(n,1,1-w)
      l = idx*d1 + (1-idx)*d0
      Mean[i,] = l
      X[i,] = rpois(n,l)
    }
    return(list(Mean=Mean,X=X,log_Mean = log(Mean),prior=prior,shape0=shape0,
                rate0=rate0,shape1=shape1,rate1=rate1,w=w,seed=seed))
  }
}

#'@title generate simulation data, t
#'@details The DGP is
#'\deqn{x_i\sim Poisson(\lambda_i),}
#'\deqn{\lambda_i\sim w\times \delta_{c} + (1-w)\times (|t_{df}|+delta_{2c}).}
#'@export

gen_data_t = function(n=1e3,n_simu=100,
                          prior='t',
                          t_df = 3,
                          t_delta=2,
                          #exp_rate = 0.5,
                          w = 0.8,
                          seed=12345){

  set.seed(seed)

  Mean = matrix(nrow=n_simu,ncol=n)
  X = matrix(nrow=n_simu,ncol=n)
  if(prior=='t'){
    for(i in 1:n_simu){
      d1 = abs(rt(n,t_df))+t_delta*2
      idx = rbinom(n,1,1-w)
      l = idx*d1 + (1-idx)*t_delta
      Mean[i,] = l
      X[i,] = rpois(n,l)
    }
    return(list(Mean=Mean,X=X,log_Mean = log(Mean),prior=prior,t_df=t_df,t_delta=t_delta,w=w,seed=seed))
  }
}

#'@title generate simulation data, exponential
#'@details The DGP is
#'\deqn{x_i\sim Poisson(\lambda_i),}
#'\deqn{\lambda_i\sim Exp(rate).}
#'@export

gen_data_exp = function(n=1e3,n_simu=100,
                      prior='exponential',
                      exp_rate = 0.5,
                      seed=12345){

  set.seed(seed)

  Mean = matrix(nrow=n_simu,ncol=n)
  X = matrix(nrow=n_simu,ncol=n)

  if(prior=='exponential'){
    for(i in 1:n_simu){
      l = rexp(n,exp_rate)
      Mean[i,] = l
      X[i,] = rpois(n,l)
    }
    return(list(Mean=Mean,X=X,log_Mean = log(Mean),prior=prior,exp_rate=exp_rate,seed=seed))
  }
}




#'@title generate simulation data, log-link
#'@details The DGP is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim w\times N(m,\sigma^2_0) + (1-w)\times N(m,\sigma^2_1).}
#'@export
gen_data_log_link = function(n=1e3,n_simu=100,w=0.8,
                             prior_mean=2,prior_var0 = 0,
                             prior_var1 = 2,seed=12345){
  set.seed(seed)

  Mu = matrix(nrow=n_simu,ncol=n)
  X = matrix(nrow=n_simu,ncol=n)
  for(i in 1:n_simu){
    d0 = rnorm(n,prior_mean,sqrt(prior_var0))
    d1 = rnorm(n,prior_mean,sqrt(prior_var1))
    idx = rbinom(n,1,1-w)
    l = idx*d1 + (1-idx)*d0
    Mu[i,] = l
    X[i,] = rpois(n,exp(l))
  }
  return(list(Mean = exp(Mu),log_Mean=Mu,X=X,prior_mean=prior_mean,
              prior_var0 = prior_var0, prior_var1 = prior_var1,w=w,seed=seed))
}




## get scores






