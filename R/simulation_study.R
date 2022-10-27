
#'@title generate simulation data
#'@export

gen_data_identity = function(n=1e3,n_simu=100,
                          prior='gamma',
                          gamma_shape = 10,
                          gamma_rate=2,
                          t_df = 3,
                          t_delta=2,
                          quasi_sparse=0.1,
                          w = 0.8,
                          seed=12345){

  set.seed(seed)

  Mean = matrix(nrow=n_simu,ncol=n)
  X = matrix(nrow=n_simu,ncol=n)
  if(prior=='gamma'){
    for(i in 1:n_simu){
      d1 = rgamma(n,gamma_shape,gamma_rate)
      idx = rbinom(n,1,1-w)
      l = idx*d1 + (1-idx)*quasi_sparse
      Mean[i,] = l
      X[i,] = rpois(n,l)
    }
    return(list(Mean=Mean,X=X,log_Mean = log(Mean),prior=prior,gamma_shape=gamma_shape,gamma_rate=gamma_rate,quasi_sparse=quasi_sparse,w=w,seed=seed))
  }

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

#'@title generate simulation data
#'@export
gen_data_exp = function(n=1e3,n_simu=100,w=0.8,prior_mean=2,prior_var0 = 0, prior_var1 = 3,seed=12345){
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
  return(list(Mean = exp(Mu),log_Mean=Mu,X=X,prior_mean=prior_mean,prior_var0 = prior_var0, prior_var1 = prior_var1,w=w,seed=seed))
}


#'@title run the algorithms
#'@export
simu_study_poisson_mean = function(sim_data,ebnm_params = list(prior_family='normal_scale_mixture'),tol=1e-5,maxiter=2e3,n_cores = 10){
  X = sim_data$X
  #Mean,log_Mean
  n_simu = nrow(X)
  n = ncol(X)
  res = mclapply(1:n_simu,function(i){


    res_GG = try(pois_mean_GG(X[i,],tol=tol,maxiter = maxiter))
    res_GMG = try(pois_mean_GMG(X[i,],tol=tol,maxiter = maxiter))
    res_GMGM = try(pois_mean_GMGM(X[i,],tol=tol,maxiter = maxiter))

    res_nb_lb = try(nb_mean_lower_bound(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
    res_nb_pg = try(nb_mean_polya_gamma(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
    res_log1exp = try(pois_mean_log1exp(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
    res_split = try(pois_mean_split(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
    res_split_mixture = try(pois_mean_split_mixture(X[i,],tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))

    res_compound = try(pois_mean_penalized_compound(X[i,]))
    res_inversion = try(pois_mean_penalized_inversion(X[i,]))

    return(list(fitted_model = list(res_GG=res_GG,
                                    res_GMG=res_GMG,
                                    res_GMGM=res_GMGM,
                                    res_nb_lb=res_nb_lb,
                                    res_nb_pg=res_nb_pg,
                                    res_log1exp=res_log1exp,
                                    res_split=res_split,
                                    res_split_mixture=res_split_mixture,
                                    res_compound=res_compound,
                                    res_inversion=res_inversion
                                    )))

  },mc.cores = n_cores)
  res
}

## get scores






