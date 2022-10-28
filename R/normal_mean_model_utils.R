#'log marginal likelihood of normal mean model
#'@param x data vector of length n
#'@param s standard error
#'@param w prior weights
#'@param mu prior mean
#'@param grid grid of sd in prior
#'@return a vector of length n
l_nm = function(x,s,w,mu,grid){
  return(log(f_nm(x,s,w,mu,grid)))
}

#'@return a n by K matrix of normal density
nm_density = function(x,s,mu,grid){
  K = length(grid)
  n = length(x)
  if(length(s)==1){
    s = rep(s,n)
  }
  sdmat = sqrt(outer(s^2,grid^2,FUN="+"))
  return(dnorm(outer(x,rep(1,K),FUN="*"),mu,sdmat))
}

#'@return a vector of likelihood, of length n
f_nm = function(x,s,w,mu,grid){
  temp = c(nm_density(x,s,mu,grid)%*%w)
  #return(temp)
  return(pmax(temp,1e-5))
}

# n = 100000
# x = 1:n
# t1 = Sys.time();a = matrix(x,nrow=n,ncol=200,byrow=FALSE);Sys.time()-t1
# t1 = Sys.time();a = outer(x,rep(1,200));Sys.time()-t1

# l_nm_d1 = function(x,w,grid){
#   s = sqrt(exp(-x))
#   f = sum(w*dnorm(x,0,sqrt(grid^2+s^2)))
#   f_d1 = sum(w*dnorm(x,0,sqrt(grid^2+s^2))*x/(grid^2+s^2))
#   f_d1/f
# }

#'@return a vector of gradient df/dz, length n
f_nm_d1_z = function(x,s,w,mu,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c(-(nm_density(x,s,mu,grid)/vmat*(x-mu))%*%w))
}

#'@return a vector of second derivative d^2f/dz^2, length n
f_nm_d2_z = function(x,s,w,mu,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,mu,grid)*((x-mu)^2/vmat^2-1/vmat))%*%w))
}

#'@return a vector of third derivative d^3f/dz^3, length n
f_nm_d3_z = function(x,s,w,mu,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,mu,grid)*(3*(x-mu)/vmat^2-(x-mu)^3/vmat^3))%*%w))
}

#'@return a vector of derivative dl_nm/dz, length n
l_nm_d1_z = function(x,s,w,mu,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d1_z(x,s,w,mu,grid)/temp)
  }
}

#'@return a vector of second derivative d^2l_nm/dz^2, length n
l_nm_d2_z = function(x,s,w,mu,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d2_z(x,s,w,mu,grid)/temp - (f_nm_d1_z(x,s,w,mu,grid)/temp)^2)
  }
}

#'@return a vector of third derivative d^3l_nm/dz^3, length n
l_nm_d3_z = function(x,s,w,mu,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d3_z(x,s,w,mu,grid)/temp - 3*f_nm_d1_z(x,s,w,mu,grid)*f_nm_d2_z(x,s,w,mu,grid)/temp^2 +2*(f_nm_d1_z(x,s,w,mu,grid)/temp)^3)
  }

}

#'@return a vector of derivative df_nm/ds2, length n
f_nm_d1_s2 = function(x,s,w,mu,grid){
  K = length(grid)
  xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,mu,grid)/vmat^2*((xmat-mu)^2-vmat))%*%w/2))
}

#'@return a vector of second derivative d^2f_nm/dzds2, length n
f_nm_d2_zs2 = function(x,s,w,mu,grid){
  K = length(grid)
  xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  return(c(((xmat-mu)*(3*vmat-(xmat-mu)^2)/vmat^3.5*exp(-(xmat-mu)^2/2/vmat))%*%w/2/sqrt(2*pi)))
}

#'@return a vector of derivative dl_nm/ds2, length n
l_nm_d1_s2 = function(x,s,w,mu,grid){
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d1_s2(x,s,w,mu,grid)/temp)
  }
}

#'@return a vector of second derivative d^2l_nm/dzds2, length n
l_nm_d2_zs2 = function(x,s,w,mu,grid){
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d2_zs2(x,s,w,mu,grid)/temp-f_nm_d1_s2(x,s,w,mu,grid)*f_nm_d1_z(x,s,w,mu,grid)/temp^2)
  }
}


#'@return a matrix of gradient, size n* K
f_nm_d1_a = function(x,s,a,mu,grid){
  n = length(x)
  dens_mat = nm_density(x,s,mu,grid)
  return((dens_mat*sum(exp(a)) - c(dens_mat%*%exp(a)))*outer(rep(1,n),exp(a))/sum(exp(a))^2)
}

#'@return a matrix of gradient, size n*K
f_nm_d2_za = function(x,s,a,mu,grid){
  K = length(grid)
  #xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  dens_mat = nm_density(x,s,mu,grid)
  lhs = c((dens_mat/vmat)%*%exp(a))
  rhs = (dens_mat/vmat)*sum(exp(a))
  return(outer(x-mu,exp(a))*(lhs-rhs)/sum(exp(a))^2)

}

#'@return a matrix of gradient, size n*K
l_nm_d1_a = function(x,s,a,mu,grid){
  w = softmax(a)
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d1_a(x,s,a,mu,grid)/temp)
  }

}

#'@return a matrix of gradient, size n*K
l_nm_d2_za = function(x,s,a,mu,grid){
  w = softmax(a)
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d2_za(x,s,a,mu,grid)/temp - f_nm_d1_a(x,s,a,mu,grid)*f_nm_d1_z(x,s,w,mu,grid)/temp^2)
  }
}

#'@return gradient w.r.t prior mean, a scalar
f_nm_d1_mu = function(x,s,w,mu,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return((c((nm_density(x,s,mu,grid)/vmat*(x-mu))%*%w)))
}

#' #'@return gradient w.r.t prior mean, a vector, do not sum over observations
#' f_nm_d1_mu_1by1 = function(x,s,w,mu,grid){
#'   vmat = outer(s^2,grid^2,FUN="+")
#'   return((c((nm_density(x,s,mu,grid)/vmat*(x-mu))%*%w)))
#' }

#'@return a vector gradient
f_nm_d2_zmu = function(x,s,w,mu,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,mu,grid)*(-(x-mu)^2/vmat^2+1/vmat))%*%w))
}

l_nm_d1_mu = function(x,s,w,mu,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d1_mu(x,s,w,mu,grid)/temp)
  }

}

l_nm_d2_zmu = function(x,s,w,mu,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,mu,grid)
  if(any(temp==0)){
    return(rep(0,length(x)))
  }else{
    return(f_nm_d2_zmu(x,s,w,mu,grid)/temp-f_nm_d1_z(x,s,w,mu,grid)*f_nm_d1_mu(x,s,w,mu,grid)/temp^2)
  }
}

#'@return derivative of l_nm(z(theta);g,s^2(theta)) w.r.t theta
l_nm_d1_theta = function(z,theta,s,w,mu,grid){
  l_nm_d1_z(z,s,w,mu,grid)*z_d1_theta(z,theta,s,w,mu,grid) + l_nm_d1_s2(z,s,w,mu,grid)*(-exp(-theta))
}

z_d1_theta = function(z,theta,s,w,mu,grid){
  numerator = 1-(-exp(-theta))*l_nm_d1_z(z,s,w,mu,grid) - exp(-theta)*(-exp(-theta))*l_nm_d2_zs2(z,s,w,mu,grid)
  denominator = 1 + exp(-theta)*l_nm_d2_z(z,s,w,mu,grid)
  return(numerator/denominator)
}

#'@return derivative of l_nm(z(theta);g,s^2(theta)) w.r.t prior (a,mu)
l_nm_d1_g = function(z,theta,s,a,mu,grid){
  w=softmax(a)
  l_nm_d1_z(z,s,w,mu,grid)*z_d1_g(z,theta,s,a,mu,grid) + cbind(l_nm_d1_a(z,s,a,mu,grid),l_nm_d1_mu(z,s,w,mu,grid))
}

z_d1_g = function(z,theta,s,a,mu,grid){
  w=softmax(a)
  n_a = -s^2*(l_nm_d2_za(z,s,a,mu,grid))
  d_a = 1+s^2*l_nm_d2_z(z,s,w,mu,grid)
  n_mu = -s^2*l_nm_d2_zmu(z,s,w,mu,grid)
  d_mu = 1+d_a
  return(cbind(n_a/d_a,n_mu/d_mu))
}

softmax = function(a){
  exp(a-max(a))/sum(exp(a-max(a)))
}

#' #'posterior mean operator
#' S = function(x,s,w,mu,grid){
#'   K = length(w)
#'   g = normalmix(pi=w,mean=rep(mu,K),sd=grid)
#'   fit.ash = ashr::ash(x,s,g=g,fixg=T)
#'   fit.ash$result$PosteriorMean
#' }
#'
#' #'posterior mean operator
#' S2 = function(x,s,w,mu,grid){
#'   lW = matrix(log(w),nrow=length(x),ncol=length(grid),byrow=T)
#'   pw = lW+dnorm(x,mean=mu,sd=outer(s^2,grid^2,FUN='+'),log=TRUE)
#'   pw = pw - apply(pw,1,max)
#'   pw = exp(pw)/rowSums(exp(pw))
#'   temp  = outer(s^2,grid^2,FUN='/')
#'   pm = x/(1+temp) + mu/(1+1/temp)
#'   return(rowSums(pw*pm))
#' }

#'posterior mean operator
S = function(x,s,w,mu,grid){
  return(x+s^2*l_nm_d1_z(x,s,w,mu,grid))
}

#'posterior mean of exp(mu) operator
S_exp = function(x,s,w,mu,grid){
  n = length(x)
  K= length(grid)
  if(length(s)==1){
    s = rep(s,n)
  }
  w = pmax(w,1e-8)
  lW = matrix(log(w),nrow=n,ncol=K,byrow=T)
  sdmat = sqrt(outer(s^2,grid^2,FUN="+"))
  pw = lW+dnorm(outer(x,rep(1,K),FUN="*"),mu,sdmat,log=TRUE)
  pw = pw - apply(pw,1,max)
  pw = exp(pw)/rowSums(exp(pw))
  pv = 1/outer(1/s^2,1/grid^2,FUN='+')
  temp  = outer(s^2,grid^2,FUN='/')
  pm = x/(1+temp) + mu/(1+1/temp)
  #browser()
  return(rowSums(pw*exp(pm+pv/2)))
}

#'posterior var operator
PV = function(x,s,w,mu,grid){
  return(1+s^2*l_nm_d2_z(x,s,w,mu,grid))
}

S_inv_obj = function(z,theta,s,w,mu,grid){
  return(z+s^2*l_nm_d1_z(z,s,w,mu,grid)-theta)
}

# S_inv_obj = function(t=0,y,theta,s,w,mu,grid,parms=NULL){
#   return(y+s^2*l_nm_d1_z(y,s,w,mu,grid)-theta)
# }


# x = rnorm(1000)
# s = rep(1,1000)
# sigma2k = c(0.1,0.2,1,2,3)
# K = length(sigma2k)
# w = rep(1/K,K)
# mu0=0
# theta = S(x,s,w,mu0,sigma2k)
# S_inv_obj_jac(x,theta,s,w,mu0,sigma2k)
# S_inv2(theta,s,w,mu0,sigma2k)

S_inv_obj_jac = function(z,theta,s,w,mu,grid){
  return(diag(c(1+s^2*l_nm_d2_z(z,s,w,mu,grid))))
}

# $root
# [1]      NaN 1.479702      NaN 1.481018 0.326163      NaN 0.326163      NaN      NaN 1.479702
#
# $f.root
# [1]           NaN  4.596237e-04           NaN -5.040602e-04 -4.440892e-16           NaN -4.440892e-16           NaN
# [9]           NaN  4.596237e-04
#
# $iter
# [1] 1000
#
# $estim.precis
# [1] NaN

#'@title bisection for root finding. Vectorized
#'@description from https://stat.ethz.ch/pipermail/r-help/2012-November/340295.html
bisection= function(f, lower, upper, ...,
                    maxiter =100,
                    tol = 1e-8) {

  fl=f(lower, ...)
  fu=f(upper, ...)
  if(any(fl*fu>0)){
    print(which(fl*fu>0))
    stop('f at lower and upper must have opposite sign')
  }
  for(n in 1:maxiter) {
    mid=(lower+upper)/2
    fmid=f(mid, ...)
    if(all(abs(fmid) < tol)){
      break
    }
    samesign = ((fmid<0)&(fl<0))|((fmid>=0)&(fl>=0))
    lower = ifelse(samesign, mid, lower)
    fl = ifelse(samesign, fmid, fl)
    upper = ifelse(!samesign, mid, upper)
    fu = ifelse(!samesign, fmid, fu)
  }
  return(mid)
}

#'@title inverse PM operator
#'@import rootSolve
#'@import nleqslv
S_inv = function(theta,s,w,mu,grid){
  #print(grid[1]==0)
  #print(all.equal(w[1],1))
  if(grid[1]==0&isTRUE(all.equal(w[1],1,tol=1e-5))){
    return(rep(mu,length(theta)))
  }else{
    n = length(theta)
    lower = ifelse(theta>=mu,theta-1,theta-s^2*10)
    upper = ifelse(theta<=mu,theta+1,theta+s^2*10)
    sol = try(bisection(S_inv_obj,lower = lower,upper = upper,s=s,w=w,mu=mu,grid=grid,theta=theta))
    if(class(sol)=='try-error'){
        sol = multiroot(S_inv_obj,start = theta,
                            #jacfunc=S_inv_obj_jac,
                            jactype = 'bandint',
                            bandup=0,banddown=0,maxiter =100,
                            theta=theta,s=s,w=w,mu=mu,grid=grid)
        #print(sol)
        if(is.nan(sol$estim.precis)|sol$iter==1){
          sol = nleqslv(x=theta,fn=S_inv_obj,jac=S_inv_obj_jac,
                        theta=theta,s=s,w=w,mu=mu,grid=grid,control = list(allowSingular=TRUE))
          return(sol$x)
        }else{
          return(sol$root)
        }
    }else{
      return(sol)
    }
  }
}


#' #' I tried supply the Jacobin matrix but there's always an error says singular matrix. But it works sometimes with numerically calculated Jacobian. Eventhough I checked thy are the same
#' #' rootSolve
#' #' nleqslv
#' S_inv = function(theta,s,w,mu,grid){
#'   if(grid[1]==0&isTRUE(all.equal(w[1],1))){
#'     return(rep(mu,length(theta)))
#'   }else{
#'     sol = multiroot(S_inv_obj,start = theta,
#'                         #jacfunc=S_inv_obj_jac,
#'                         jactype = 'bandint',
#'                         bandup=0,banddown=0,maxiter =100,
#'                         theta=theta,s=s,w=w,mu=mu,grid=grid)
#'     #print(sol)
#'     if(is.nan(sol$estim.precis)|sol$iter==1){
#'       sol = nleqslv(x=theta,fn=S_inv_obj,jac=S_inv_obj_jac,
#'                     theta=theta,s=s,w=w,mu=mu,grid=grid,control = list(allowSingular=TRUE))
#'       if(sol$iter==1){
#'         return(S_inv_loop(theta,s,w,mu,grid))
#'       }else{
#'         return(sol$x)
#'       }
#'
#'     }else{
#'       return(sol$root)
#'     }
#'   }
#' }


# S_inv3 = function(theta,s,w,mu,grid){
#   temp_obj = function(z,theta,s,w,mu,grid){
#     return(sum(z^2/2+s^2*l_nm(z,s,w,mu,grid)-z*theta))
#   }
#   opt = optim(theta,temp_obj,gr = S_inv_obj,theta=theta,s=s,w=w,mu=mu,grid=grid,method = 'L-BFGS-B',control = list(fnscale=-1))
#   return(opt)
# }

#' #########################S_inv using uniroot in a for loop. Slow###############'@title inverse operator of S
#' #'  S^{-1}(theta) returns the z such that S(z) = theta
#' S_inv_loop = function(theta,s,w,mu,grid,z_range=NULL){
#'
#'   # obj = function(z,theta,s,w,mu,grid){
#'   #   return(z+s^2*l_nm_d1_z(z,s,w,mu,grid)-theta)
#'   # }
#'
#'   if(grid[1]==0&isTRUE(all.equal(w[1],1))){
#'     return(rep(mu,length(theta)))
#'   }else{
#'     if(is.null(z_range)){
#'       z_range = c(-10,10)
#'     }
#'     n = length(theta)
#'     z_out = double(n)
#'     for(j in 1:n){
#'       #print(j)
#'       if(theta[j]>=0){
#'         # z_out[j] = uniroot(obj,c(theta[j],theta[j]/median(1/(1+s[j]^2/grid^2))),
#'         #                    theta=theta[j],s=s[j],w=w,grid=grid,extendInt = 'upX')$root
#'         z_out[j] = uniroot(S_inv_obj,c(theta[j],z_range[2]),
#'                            theta=theta[j],s=s[j],w=w,grid=grid,mu=mu,extendInt = 'upX')$root
#'       }else{
#'         # z_out[j] = uniroot(obj,c(theta[j]/median(1/(1+s[j]^2/grid^2)),theta[j]),
#'         #                    theta=theta[j],s=s[j],w=w,grid=grid,extendInt = 'upX')$root
#'         z_out[j] = uniroot(S_inv_obj,c(z_range[1],theta[j]),
#'                            theta=theta[j],s=s[j],w=w,grid=grid,mu=mu,extendInt = 'upX')$root
#'       }
#'
#'     }
#'     return(z_out)
#'   }


#'}
#'
#' #' inverse operator of S, parallel version using mclapply
#' #'  S^{-1}(theta) returns the z such that S(z) = theta
#'
#' library(parallel)
#' S_inv_parallel = function(theta,s,w,mu,grid,z_range=NULL){
#'
#'   obj = function(z,theta,s,w,mu,grid){
#'     return(z+s^2*l_nm_d1_z(z,s,w,mu,grid)-theta)
#'   }
#'
#'   f = function(i,theta,s,w,mu,grid,z_range){
#'     if(theta[i]>=0){
#'       return(uniroot(obj,c(theta[i],z_range[2]),
#'                          theta=theta[i],s=s[i],w=w,grid=grid,mu=mu,extendInt = 'upX')$root)
#'     }else{
#'      return(uniroot(obj,c(z_range[1],theta[i]),
#'                          theta=theta[i],s=s[i],w=w,grid=grid,mu=mu,extendInt = 'upX')$root)
#'     }
#'   }
#'
#'   if(is.null(z_range)){
#'     z_range = range(theta) + c(-5,5)
#'   }
#'
#'   n = length(theta)
#'   z_out = simplify2array(mclapply(1:n,f,theta=theta,s=s,w=w,mu=mu,grid=grid,z_range=z_range,mc.cores = 10))
#'   z_out
#' }

#'@title compound penalty function of ebnm
nm_penalty_compound = function(theta,s,w,mu,grid){
  return(-l_nm(theta,s,w,mu,grid) - (theta-S(theta,s,w,mu,grid))^2/2/s^2)
}
#'@title gradient of compound penalty function
#'@description -l'(theta) - s^2l'(theta)l''(theta)
nm_penalty_compound_grad = function(theta,s,w,mu,grid){
  return(-l_nm_d1_z(theta,s,w,mu,grid) - s^2*l_nm_d1_z(theta,s,w,mu,grid)*l_nm_d2_z(theta,s,w,mu,grid))
}
#'@title hessian of compound penalty function
#'@description -l''(theta) - s^2(l''(theta)^2+l'''(theta)l'(theta))
nm_penalty_compound_hess = function(theta,s,w,mu,grid){
  return(-l_nm_d2_z(theta,s,w,mu,grid) - s^2*(l_nm_d2_z(theta,s,w,mu,grid)^2+l_nm_d3_z(theta,s,w,mu,grid)*l_nm_d1_z(theta,s,w,mu,grid)))
}

#'@title penalty function of ebnm
nm_penalty = function(theta,s,w,mu,grid){
  # if(is.null(z_range)){
  #   z_range = range(theta) + c(-1,1)
  # }
  z = S_inv(theta,s,w,mu,grid)
  original_penalty = nm_penalty_compound(z,s,w,mu,grid)
  return(original_penalty)
}

#'@title gradient of penalty function of ebnm
nm_penalty_grad = function(theta,s,w,mu,grid){
  # if(is.null(z_range)){
  #   z_range = range(theta) + c(-1,1)
  # }
  z = S_inv(theta,s,w,mu,grid)
  return((z-theta)/s^2)
}

#'@title hessian of penalty function of ebnm
nm_penalty_hess = function(theta,s,w,mu,grid){
  # if(is.null(z_range)){
  #   z_range = range(theta) + c(-1,1)
  # }
  z = S_inv(theta,s,w,mu,grid)
  temp = l_nm_d2_z(z,s,w,mu,grid)
  return(-temp/(1+s^2*temp))
}
