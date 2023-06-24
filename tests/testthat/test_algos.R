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


set.seed(12345)
n = 1000
#pi0=0.9
#mu = 5
#b = c(rep(mu,n*pi0),rnorm(n-n*pi0,mu,1))
#b = c(rep(0,n*pi0),rep(mu,n-n*pi0))
#x = rpois(n,exp(b))
x = rpois(n,10)
fit = pois_mean_split(x,ebnm_params = list(prior_family='normal'),mu_pm_init = NULL,sigma2 = 1,tol=1e-100)
fit$fitted_g$sigma2_trace
fit$fitted_g$g_b

fit$posterior$mean_log[1:5]
fit$fit$posterior$mean[1:5]

sqrt(fit$posterior$var_log[1:5])
fit$fit$posterior$sd[1:5]

plot(x,col='grey80')
lines(fit$posterior$mean_exp_b)
lines(fit$posterior$mean)

plot(b,col='grey80')
lines(fit$posterior$mean_b)

fit_true = pois_mean_GMGM(x,mixsd=fit$fitted_g$g_b$sd[-1],tol=1e-10)

plot(fit$posterior$mean_log,fit$posterior$mean_b)
abline(a=0,b=1)

hist(fit$posterior$var_b)
hist(fit$posterior$var_log)


fit2 = pois_mean_split(x,ebnm_params = list(prior_family='point_normal',g_init=ashr::normalmix(pi=c(0.8,0.2),mean=c(5,5),sd=c(0,1)),fix_g=TRUE),
                       mu_pm_init = b,sigma2 = 1,tol=1e-10)
fit2$fitted_g$sigma2_trace
fit2$fitted_g$g_b

plot(x,col='grey80')
lines(fit2$posterior$mean_exp_b)
lines(fit2$posterior$mean,col='grey80')










