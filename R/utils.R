#'@title MSE
#'@export
mse = function(x,y){mean((x-y)^2)}

#'@title RMSE
#'@export
rmse = function(x,y){sqrt(mean((x-y)^2))}

#'@title MAE
#'@export
mae = function(x,y){mean(abs(x-y))}


#'@title select grids
#'@param x poisson data
#'@param s scaling
#'@export
select_mixsd = function (x,s,mult=sqrt(2),mode=0){
  betahat = log(1/s+x/s) - mode
  sebetahat = sqrt(1/(1/s+x/s))
  sigmaamin = min(sebetahat)/10
  sigmaamin = max(sigmaamin,3.2e-3)
  sigmaamax = max(x/s)
  sigmaamax = log(sigmaamax + 3*sqrt(sigmaamax))/2
  npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
  return(mult^((-npoint):0) * sigmaamax)
}
