#'@title Optimize vga poisson problem
#'@description This function tries Newton's method. If not working, then use bisection.
#'@param init_val initial value for posterior mean
#'@param x,s data and scale factor
#'@param beta,sigma2 prior mean and variance. Their length should be equal to n=length(x)
#'@export
vga_pois_solver = function(init_val,x,s,beta,sigma2,maxiter=1000,tol=1e-5,method = 'newton'){

  n = length(x)
  if(length(sigma2)==1){
    sigma2 = rep(sigma2,n)
  }
  if(length(beta)==1){
    beta = rep(beta,n)
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  if(method=='newton'){
    # use Newton's method fist
    res = try(vga_pois_solver_Newton(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent = TRUE)
    if(class(res)=='try-error'){
      # If Newton failed, use bisection
      res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
      if(class(res)=='try-error'){
        # If bisection also failed, return initial  values with a warning.
        warnings('Both Newton and Bisection methods failed. Returning initial values.')
        return(list(m = init_val,v = sigma2/(sigma2*x+beta+1-init_val)))
      }else{
        return(res)
      }
    }else{
      return(res)
    }
  }else if(method=='bisection'){
    res = try(vga_pois_solver_bisection(x,s,beta,sigma2,maxiter=maxiter,tol=tol),silent=TRUE)
    if(class(res)=='try-error'){
      return(list(m = init_val,v = sigma2/(sigma2*x+beta+1-init_val)))
    }else{
      return(res)
    }
  }else{
    stop('Only Newton and bisection are supported.')
  }


}

h_v = function(v,x,s,beta,sigma2){
  val = (-s*exp(sigma2*x+beta+1-sigma2/v + v/2) - 1/sigma2 + 1/v)
  return(-val)
}

# h_m = function(v,x,s,beta,sigma2){
#   m = sigma2*x + beta + 1 - sigma2/v
#   val = x - s*exp(m+v/2) - (m-beta)/sigma2
#   return(-val)
# }


# vga_pois_solver_m = function(init_val,x,s,beta,sigma2,maxiter=1000,tol=1e-8){
#
#   n = length(x)
#   if(length(sigma2)==1){
#     upper = rep(sigma2,n)
#   }else if(length(sigma2)==n){
#     upper = sigma2
#   }else{
#     stop('check length of sigma2')
#   }
#   return(vga_optimize_m(init_val,x,s,beta,sigma2))
#
# }

vga_pois_solver_bisection = function(x,s,beta,sigma2,maxiter=1000,tol=1e-5){
  n = length(x)
  if(length(sigma2)==1){
    upper = rep(sigma2,n)
  }else if(length(sigma2)==n){
    upper = sigma2
  }else{
    stop('check length of sigma2')
  }
  v = bisection(h_v,
                lower = rep(0,n),upper = upper,
                x=x,s=s,beta=beta,sigma2=sigma2,
                auto_adjust_interval = FALSE,
                maxiter=maxiter,tol=tol)
  m = sigma2*x + beta + 1 - sigma2/v
  return(list(m=m,v=v))
}

vga_pois_solver_Newton = function(m,x,s,beta,sigma2,maxiter=1000,tol=1e-5){

  const0 = sigma2*x+beta + 1
  const1 = 1/sigma2
  const2 = sigma2/2

  # make sure m < sigma2*x+beta
  idx = (m>(const0-1))
  if(sum(x)>0){
    m[idx] =suppressWarnings(vga_pois_solver_bisection(x[idx],s[idx],beta[idx],sigma2[idx],maxiter = 10)$m)
  }


  for(i in 1:maxiter){

    temp = (const0-m)
    sexp = s*exp(m+const2/temp)
    f = x - sexp - (m-beta)/sigma2
    if(max(abs(f))<tol){
      break
    }
    f_grad = -sexp*(1+const2/temp^2)-const1
    m = m - f/f_grad
  }
  if(i>=maxiter){
    warnings('Newton method not converged yet.')
  }
  return(list(m=m,v=sigma2/(const0-m)))

}

#
# vga_pos_Newton_f = function(m,x,s,beta,sigma2){
#   return(x - s*exp(m+v_m(m,x,beta,sigma2)/2) - (m-beta)/sigma2)
# }
#
# vga_pos_Newton_f_grad = function(m,x,s,beta,sigma2){
#   return(-s*exp(m+v_m(m,x,beta,sigma2)/2)*(1+v_dm(m,x,beta,sigma2)/2)-1/sigma2)
# }

####################################
####################################
#################3 notes ###########
# I tried multiroot in rootSOlve and it's not faster than mine implementation.
