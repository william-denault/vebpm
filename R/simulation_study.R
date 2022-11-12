
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


#'@title run the algorithms
#'@import parallel
#'@importFrom ebpm ebpm_gamma
#'@importFrom ebpm ebpm_exponential_mixture
#'@importFrom ashr ash_pois
#'@export
simu_study_poisson_mean = function(sim_data,
                                   ebnm_params = list(prior_family='normal_scale_mixture'),
                                   tol=1e-6,maxiter=2e3,n_cores = 10,
                                   method_list = c('GG','GMGM',
                                                   'GMGM_pointmass',
                                                   'nb_pg',
                                                   'log1exp','split','split_mixture',
                                                   'penalty_compound',
                                                   'penalty_inversion',
                                                   'ash_pois_identity',
                                                   'ebpm_gamma',
                                                   'ebpm_exp_mixture'),
                                   save_data = TRUE){
  X = sim_data$X
  #Mean,log_Mean
  n_simu = nrow(X)
  n = ncol(X)
  n_method = length(method_list)
  res = mclapply(1:n_simu,function(i){


    fitted_model <- vector("list", n_method)
    names(fitted_model) <- method_list

    if('GG'%in%method_list){
      res_GG = try(pois_mean_GG(X[i,],tol=tol,maxiter = maxiter))
      fitted_model$GG = res_GG
      #MSE_mean$GG = mse(res_GG$posterior$posteriorMean_mean)
    }
    # if('GMG'%in%method_list){
    #   res_GMG = try(pois_mean_GMG(X[i,],tol=tol,maxiter = maxiter))
    #   fitted_model$GMG = res_GMG
    # }

    if('GMGM'%in%method_list){
      res_GMGM = try(pois_mean_GMGM(X[i,],tol=tol,maxiter = maxiter,point_mass = F))
      fitted_model$GMGM = res_GMGM
    }

    if('GMGM_pointmass'%in%method_list){
      res_GMGM_pointmass = try(pois_mean_GMGM(X[i,],tol=tol,maxiter = maxiter,point_mass = T))
      fitted_model$GMGM_pointmass = res_GMGM_pointmass
    }


    # if('nb_lb'%in%method_list){
    #   res_nb_lb = try(nb_mean_lower_bound(X[i,],r=1000,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
    #   fitted_model$nb_lb = res_nb_lb
    # }

    if('nb_pg'%in%method_list){
      res_nb_pg = try(nb_mean_polya_gamma(X[i,],r=max(100,median(X[i,])),tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
      fitted_model$nb_pg = res_nb_pg
    }

    if('log1exp'%in%method_list){
      res_log1exp = try(pois_mean_log1exp(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
      fitted_model$log1exp = res_log1exp
    }

    if('split'%in%method_list){
      res_split = try(pois_mean_split(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
      fitted_model$split = res_split
    }

    if('split_mixture'%in%method_list){
      res_split_mixture = try(pois_mean_split_mixture(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
      fitted_model$split_mixture = res_split_mixture
    }

    if('penalty_compound'%in%method_list){
      res_compound = try(pois_mean_penalized_compound(X[i,],tol=tol,maxiter=maxiter))
      fitted_model$penalty_compound = res_compound
    }

    if('penalty_inversion'%in%method_list){
      res_inversion = try(pois_mean_penalized_inversion(X[i,],tol=tol,maxiter=maxiter))
      fitted_model$penalty_inversion = res_inversion
    }

    if('ash_pois_identity'%in%method_list){
      t_start = Sys.time()
      res_ash_pois_identity = try(ash_pois(X[i,],link='identity'))
      t_end=Sys.time()
      res_ash_pois_identity$posterior = list(mean=res_ash_pois_identity$result$PosteriorMean,mean_log=log(res_ash_pois_identity$result$PosteriorMean))
      res_ash_pois_identity$run_time = difftime(t_end,t_start,units='secs')
      res_ash_pois_identity$elbo = res_ash_pois_identity$loglik
      fitted_model$ash_pois_identity = res_ash_pois_identity
    }

    # if('ash_pois_log'%in%method_list){
    #   res_ash_pois_log = try(ash_pois(X[i,],link='log'))
    #   res_ash_pois_log$posterior = list(mean=res_ash_pois_log$result$PosteriorMean,mean_log=log(res_ash_pois_log$result$PosteriorMean))
    #   fitted_model$ash_pois_log = res_ash_pois_log
    # }

    if('ebpm_gamma'%in%method_list){
      t_start = Sys.time()
      res_ebpm_gamma = try(ebpm_gamma(X[i,]))
      t_end=Sys.time()
      res_ebpm_gamma$run_time = difftime(t_end,t_start,units='secs')
      res_ebpm_gamma$elbo = res_ebpm_gamma$log_likelihood
      fitted_model$ebpm_gamma = res_ebpm_gamma
    }

    if('ebpm_exp_mixture'%in%method_list){
      t_start = Sys.time()
      res_ebpm_exp_mixture = try(ebpm_exponential_mixture(X[i,]))
      t_end=Sys.time()
      res_ebpm_exp_mixture$run_time = difftime(t_end,t_start,units='secs')
      res_ebpm_exp_mixture$elbo = res_ebpm_exp_mixture$log_likelihood
      fitted_model$ebpm_exp_mixture = res_ebpm_exp_mixture
    }

    MSE_mean = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
        mse(x$posterior$mean,sim_data$Mean[i,])
      }else{
        NA
      }

    }))

    MAE_mean = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
      mae(x$posterior$mean,sim_data$Mean[i,])
      }else{
        NA
      }
    }))

    MSE_log_mean = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
      mse(x$posterior$mean_log,sim_data$log_Mean[i,])
      }else{
        NA
      }
    }))

    MAE_log_mean = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
      mae(x$posterior$mean_log,sim_data$log_Mean[i,])
      }else{
        NA
      }
    }))

    run_times = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
        x$run_time
      }else{
        NA
      }
    }))

    elbos = simplify2array(lapply(fitted_model,function(x){
      if(class(x)!='try-error'){
        x$elbo
      }else{
        NA
      }
    }))


    return(list(fitted_model = fitted_model,
                MSE_mean=MSE_mean,
                MAE_mean=MAE_mean,
                MSE_log_mean=MSE_log_mean,
                MAE_log_mean=MAE_log_mean,
                run_times=run_times,
                elbos=elbos))

  },mc.cores = n_cores)
  if(save_data){
    return(list(sim_data=sim_data,output = res))
  }else{
    return(res)
  }

}

## get scores






