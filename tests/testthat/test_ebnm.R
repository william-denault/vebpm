set.seed(12345)
n = 200
w0 = 0.9
mu = c(rep(0,round(n*w0)),rep(10,n-round(n*w0)))
w_true = c(w0,1-w0)
grid_true = c(0.01,7)
s = rep(1,n)
y = rnorm(n,mu,s)
library(ashr)
grids = ebnm:::get_ashr_grid(y,s)

fit.ash = ash(y,s,mixcompdist = 'normal',pointmass=FALSE,prior='uniform',mixsd=grids)
plot(y,main='ash fit',col='grey80')
lines(mu,col='grey60')
lines(fit.ash$result$PosteriorMean)
legend('topleft',c('data','true mean','ash posteriorMean'),
       pch=c(1,NA,NA),lty=c(NA,1,1),col=c('grey80','grey60',1))
fit.ash$loglik

fit_inv = ebnm_penalized_inversion(y,s,g_init = list(sd=grids))
plot(y,main='inversion',col='grey80',pch=20)
lines(mu,col='grey60')
lines(fit_inv$posterior$mean)

fit_compound = ebnm_penalized_inversion(y,s,g_init = list(sd=grids))
fit_compound$log_likelihood
plot(y,main='compound',col='grey80',pch=20)
lines(mu,col='grey60')
lines(fit_compound$posterior$mean)
