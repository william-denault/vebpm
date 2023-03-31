#'@title ebnm penalized compound
#'@param x,s observations, and sd
#'@param g_init normalmix object, list(pi,mean,sd)
#'@param fix_g FALSE
#'@param mode currently only supports mode = 0
#'@param scale default to 'estimate'
#'@param control NULL
#'@return posterior, fitted_g, log_likelihood
#'@export
ebnm_penalized_compound = function(x,s,mode = 0,
                                   scale = 'estimate',
                                   g_init=NULL,
                                   fix_g=FALSE,
                                   control=list(opt_method = 'L-BFGS-B')){
  n = length(x)
  if(length(s)==1){
    s = rep(s,n)
  }
  if(is.null(g_init$sd)){
    grid = ebnm:::get_ashr_grid(x,s,mode)
  }else{
    grid = g_init$sd
  }
  K = length(grid)
  if(is.null(g_init$pi)){
    w_init = rep(1/K,K)
  }else{
    w_init = g_init$pi
  }
  z_init = x
  out = optim(c(z_init,w_init),
              fn=f_obj_ebnm_compound,
              method=control$opt_method,
              y=y,grid=grid,s=s)
  z = out$par[1:n]
  w = softmax(out$par[-(1:n)])
  posteriorMean = S(z,s,w,0,grid)
  return(list(posterior = list(mean=posteriorMean),
              fitted_g=list(pi=w,sd=grid),
              log_likelihood = -out$value,
              opt = out))
}


#'objective function
#'@param theta (mu_bar,w)
#'@param grid prior sds
f_obj_ebnm_compound = function(theta,y,s,grid){
  n = length(y)
  w = softmax(theta[-(1:n)])
  z = theta[1:n]
  res = sum((y-z-s^2*l_nm_d1_z(z,s,w,0,grid))^2/2/s^2 - l_nm(z,s,w,0,grid) - s^2*(l_nm_d1_z(z,s,w,0,grid))^2/2)
  return(res)
}

#' #'objective function
#' #'@param theta (mu_bar,w)
#' #'@param grid prior sds
#' f_obj_grad_ebnm_compound = function(theta,y,s,grid){
#'   n = length(y)
#'   a = theta[-(1:n)]
#'   w = softmax(a)
#'   z = theta[1:n]
#'
#'   grad_z = (1+s^2*l_nm_d2_z(z,s,w,grid))*(z-y)/s^2
#'   grad_a = colSums((s^2*l_nm_d1_z(z,s,w,grid)-y+z)*l_nm_d2_za(z,s,a,grid) - l_nm_d1_a(z,s,a,grid) - s^2*l_nm_d1_z(z,s,w,grid)*l_nm_d2_za(z,s,a,grid))
#'
#'   return(c(grad_z,c(grad_a)))
#' }
