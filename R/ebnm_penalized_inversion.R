#'@title ebnm penalized version
#'@param x,s observations, and sd
#'@param g_init normalmix object, list(pi,mean,sd)
#'@param fix_g FALSE
#'@param mode currently only supports mode = 0
#'@param scale default to 'estimate'
#'@param control NULL
#'@return posterior, fitted_g, log_likelihood
#'@export
ebnm_penalized_inversion = function(x,s,
                                    mode = 0,
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

  theta_init = rep(0,n)

  params = c(theta_init,w_init)
  out = optim(params,
              fn=f_obj_ebnm_inversion,
              gr=f_obj_grad_ebnm_inversion,
              y=x,
              s=s,
              grid=grid,
              method = control$opt_method)

  return(list(posterior=list(mean = out$par[1:n]),
              fitted_g = list(pi = softmax(out$par[-(1:n)]),sd=grid),
              log_likelihood = -out$value,
              opt = out))
}

#'objective function ebnm_compound
#'@param theta (theta,w)
#'@param grid prior sds
f_obj_ebnm_inversion = function(params,y,s,grid){
  n = length(y)
  w = softmax(params[-(1:n)])
  theta = params[1:n]
  z = S_inv(theta,s,w,0,grid)
  return(sum((y-theta)^2/2/s^2 - l_nm(z,s,w,0,grid)-(z-theta)^2/2/s^2))
}


#'objective function gradient ebnm_compound
#'@param theta (theta,w)
#'@param grid prior sds
f_obj_grad_ebnm_inversion = function(params,y,s,grid){
  n = length(y)
  a = params[-(1:n)]
  w = softmax(a)
  theta = params[1:n]
  z = S_inv(theta,s,w,0,grid)
  grad_theta = (z-y)/s^2
  grad_a = -colSums(l_nm_d1_a(z,s,a,0,grid))
  return(c(grad_theta,c(grad_a)))
}
