#set.seed(1) error in this seed for penalized-inversion method
set.seed(1234)
x = rpois(10,5)
tol = 1e-3
maxiter = 100
ebnm_params=list(prior_family='normal_scale_mixture')

res_GG = (pois_mean_GG(x,tol=tol,maxiter = maxiter))
#res_GMG = (pois_mean_GMG(x,tol=tol,maxiter = maxiter))
res_GMGM = (pois_mean_GMGM(x,tol=tol,maxiter = maxiter))

res_nb_lb = (nb_mean_lower_bound(x,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
res_nb_pg = (nb_mean_polya_gamma(x,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
res_log1exp = (pois_mean_log1exp(x,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
res_split = (pois_mean_split(x,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))
res_split_mixture = (pois_mean_split_mixture(x,tol=tol,maxiter = maxiter,ebnm_params = ebnm_params))

res_compound = (pois_mean_penalized_compound(x))
res_inversion = pois_mean_penalized_inversion(x,tol=1e-4)



#simdata = gen_data_gamma(n=100,n_simu=2)
#out= simu_study_poisson_mean(simdata,n_cores = 1)
