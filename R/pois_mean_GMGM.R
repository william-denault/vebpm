#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian Mixture prior, Gaussian Mixture posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param w prior weights
#'@param prior_mean prior mean
#'@param mixsd prior sd
#'@param point_mass whether put a point-mass in the prior
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posterior:}{posteriorMean_latent is the posterior mean of mu, posteriorMean_mean is the posterior of exp(mu)}
#'  \item{fitted_g:}{estimated prior}
#'  \item{obj_value:}{objective function values}
#'@examples
#'\dontrun{
#'n = 1000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'pois_mean_GMGM(x)
#'}
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim \sum_k w_k N(\beta,\sigma_k^2).}
#'@export

pois_mean_GMGM = function(x,
                          s=NULL,
                          w = NULL,
                          prior_mean = NULL,
                          mixsd=NULL,
                          point_mass = TRUE,
                          optim_method = 'L-BFGS-B',
                          maxiter = 1000,
                          tol = 1e-5){

  n = length(x)

  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  # init prior mean,
  if(is.null(prior_mean)){
    beta = log(sum(x)/sum(s))
    est_beta = TRUE
  }else{
    est_beta = FALSE
    beta = prior_mean
  }

  if(is.null(mixsd)){

    ## how to choose grid in this case?
    ## use ebnm method
    sigma2k = ashr:::autoselect.mixsd(data=list(x = log(0.1/s+x/s),s = sqrt(1/(0.1/s+x/s)),lik=list(name='normal')),sqrt(2),mode=0,grange=c(-Inf,Inf),mixcompdist = 'normal')^2
    #sigma2k = (ebnm:::default_smn_scale(log(x/s+1),sqrt(1/(x/s+1)),mode=beta)[-1])^2
    #sigma2k = c(1e-10,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
    #sigma2k = c(1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
  }else{
    sigma2k = mixsd^2
  }
  K = length(sigma2k)


  M = matrix(log(1+x),nrow=n,ncol=K,byrow = F)
  V = matrix(1/n,nrow=n,ncol=K)
  Sigma2k = matrix(sigma2k,nrow=n,ncol=K,byrow=T)
  X = matrix(x,nrow=n,ncol=K,byrow=F)

  # init prior weights
  if(is.null(w)){
    if(point_mass){
      w0 = 1/(K+1)
      w= rep(1/(K+1),K)
    }else{
      w= rep(1/K,K)
      w0 = 0
    }
  }else{
    if(point_mass){
      w0=w[1]
      w= w[-1]
    }else{
      w0=0
    }
  }
  # const in objective function
  const = sum((x-1)*log(s)) - sum(lfactorial(x))

  qz = matrix(0,nrow=n,ncol=K)
  qz0 = rep(0,n)

  obj = rep(0,maxiter+1)
  obj[1] = -Inf

  for(iter in 1:maxiter){

    # #update posterior mean, variances
    # # this can be paralleled?
    # # this is too slow, need a vectorized version
    # for(i in 1:n){
    #   for (k in 1:K) {
    #     temp = pois_mean_GG1(x[i],s[i],beta,sigma2k[k],optim_method,M[i,k],V[i,k])
    #     M[i,k] = temp$m
    #     V[i,k] = temp$v
    #   }
    # }

    # for each K, solve a vectorized version
    for(k in 1:K){
      opt = vga_optimize(c(M[,k],V[,k]),x,s,beta,sigma2k[k])
      M[,k] = opt$m
      V[,k] = opt$v
    }
    # update posterior weights
    qz = X*M-s*exp(M+V/2)+matrix(log(w),nrow=n,ncol=K,byrow=T)-log(Sigma2k)/2-(M^2+V-2*M*beta+beta^2)/Sigma2k/2 + log(V)/2
    if(point_mass){
      qz0 = x*beta - s*exp(beta) + log(w0)
      qz = cbind(qz0,qz)
    }
    #browser()
    qz = qz - apply(qz,1,max)
    qz = exp(qz)
    qz = qz/rowSums(qz)
    qz = pmax(qz,1e-15)

    if(point_mass){
      qz0 = qz[,1]
      qz = qz[,-1]
    }

    # update prior

    if(est_beta){
      if(point_mass){
        beta = optimize_prior_mean_point_mass(beta,x,s,M,qz,qz0,Sigma2k,optim_method)
      }else{
        beta = sum(M*qz/Sigma2k)/sum(qz/Sigma2k)
      }
    }

    w = colMeans(qz)
    w = pmax(w, 1e-15)
    w0 = mean(qz0)
    w0 = pmax(w0,1e-15)


    obj[iter+1] = pois_mean_GMGM_obj(X,x,s,M,V,w,beta,Sigma2k,qz,point_mass,w0,qz0,const)
    if((obj[iter+1] - obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }

  }

  return(list(posterior = list(mean_log = rowSums(qz*M) + qz0*beta,
                               #2nd_moment_log = rowSums(qz*(M^2+V)) + qz0*beta^2,
                               mean = rowSums(qz*exp(M + V/2))+ qz0*exp(beta)),
              fitted_g = list(weight = c(w0,w),mean=beta,var=c(0,sigma2k)),
              obj_value=obj,
              fit = list(M=M,V=V,qz=qz,qz0=qz0)))

  # if(point_mass){
  #   return(list(posteriorMean=rowSums(qz*M) + qz0*beta,posterior2nd_moment= rowSums(qz*(M^2+V)) + qz0*beta^2,priorMean=beta,w=w,w0=w0,M=M,V=V,obj_value=obj,qz=qz,qz0=qz0))
  # }else{
  #   return(list(posteriorMean=rowSums(qz*M),posterior2nd_moment= rowSums(qz*(M^2+V)),M=M,V=V,obj_value=obj,w=w,qz=qz,priorMean=beta))
  # }



}

pois_mean_GMGM_obj = function(X,x,s,M,V,w,beta,Sigma2k,qz,point_mass,w0,qz0,const){
  n = dim(X)[1]
  K = length(w)
  lW = matrix(log(w),nrow=n,ncol=K,byrow=T)

  if(point_mass){
    return(sum(qz*(X*M-s*exp(M+V/2)+lW-log(Sigma2k)/2-(M^2+V-2*M*beta+beta^2)/2/Sigma2k-log(qz)+log(V)/2)) + sum(qz0*(x*beta-s*exp(beta)))+sum(qz0*log(w0))-sum(qz0*log(qz0))+ const)
  }else{
    return(sum(qz*(X*M-s*exp(M+V/2)+lW-log(Sigma2k)/2-(M^2+V-2*M*beta+beta^2)/2/Sigma2k-log(qz)+log(V)/2))+const)
  }

}

optimize_prior_mean_point_mass = function(beta_init,x,s,M,qz,qz0,Sigma2k,optim_method){
  opt = optim(beta_init,
              fn = optimize_prior_mean_point_mass_obj,
              gr = optimize_prior_mean_point_mass_obj_gradient,
              x=x,
              s=s,
              M=M,
              Sigma2k=Sigma2k,
              qz=qz,
              qz0=qz0,
              method = optim_method)

  return(opt$par)
}

optimize_prior_mean_point_mass_obj = function(beta,x,s,M,qz,qz0,Sigma2k){
  val = sum(qz0*(x*beta-s*exp(beta))) + beta*sum(qz/Sigma2k*M) - sum(qz/Sigma2k)/2*beta^2
  return(-val)
}

optimize_prior_mean_point_mass_obj_gradient = function(beta,x,s,M,qz,qz0,Sigma2k){
  val = sum(qz0*x) -sum(qz0*s)*exp(beta) + sum(qz/Sigma2k*M) - sum(qz/Sigma2k)*beta
  return(-val)
}



