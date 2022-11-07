#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian Mixture prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param w prior weights
#'@param prior_mean prior mean
#'@param mixsd prior sd grids
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posteriorMean:}{posterior mean}
#'  \item{posteriorVar:}{posterior var}
#'  \item{obj_value:}{objective function values}
#'  \item{priorMean:}{prior mean}
#'  \item{sigma2:}{prior variance}
#'@examples
#'n = 10000
#'mu = rnorm(n)
#'x = rpois(n,exp(mu))
#'pois_mean_GMG(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim \sum_k N(\beta,\sigma_k^2).}
#'@export

pois_mean_GMG = function(x,
                         s = NULL,
                         w = NULL,
                         prior_mean = NULL,
                         mixsd=NULL,
                         optim_method = 'L-BFGS-B',
                         maxiter = 1000,
                         tol = 1e-5,
                         m_init = NULL,
                         v_init = NULL){
  n = length(x)

  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  if(is.null(prior_mean)){
    beta = log(sum(x)/sum(s))
    est_beta=TRUE
  }else{
    est_beta=FALSE
    beta = prior_mean
  }


  if(is.null(mixsd)){
    sigma2k = ashr:::autoselect.mixsd(data=list(x = log(0.1/s+x/s),s = sqrt(1/(0.1/s+x/s)),lik=list(name='normal')),sqrt(2),mode=0,grange=c(-Inf,Inf),mixcompdist = 'normal')^2
    #sigma2k = (ebnm:::default_smn_scale(log(x/s+1),sqrt(1/(x/s+1)),mode=beta)[-1])^2
    if(min(sigma2k)>1e-4){
      sigma2k =c(1e-4,sigma2k)
    }
    #sigma2k = c(1e-10,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
    #sigma2k = c(1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
  }else{
    sigma2k = mixsd^2
  }
  K = length(sigma2k)

  if(is.null(m_init)){
    m = log(x+1)
  }else{
    m = m_init
  }

  if(is.null(v_init)){
    v = rep(1/n,n)
  }else{
    v = v_init
  }


  Sigma2k = matrix(sigma2k,nrow=n,ncol=K,byrow=T)
  if(is.null(w)){
    w = rep(1/K,K)
  }

  qz = matrix(0,nrow=n,ncol=K)

  obj = rep(0,maxiter+1)
  obj[1] = -Inf

  for(iter in 1:maxiter){

    qz = matrix(log(w),nrow=n,ncol=K,byrow=T)-log(Sigma2k)/2-matrix(m^2+v-2*m*beta+beta^2,nrow=n,ncol=K,byrow=F)/Sigma2k/2
    qz = qz - apply(qz,1,max)
    qz = exp(qz)
    qz = qz/rowSums(qz)
    qz = pmax(qz,1e-15)

    # for(i in 1:n){
    #   temp = pois_mean_GMG1(x[i],s[i],qz[i,],beta,sigma2k,optim_method,m[i],v[i])
    #   m[i] = temp$m
    #   v[i] = temp$s2
    # }

    opt = optim(c(m,log(v)),
                fn = pois_mean_GMG_opt_obj,
                gr = pois_mean_GMG_opt_obj_gradient,
                x=x,
                s=s,
                qz=qz,
                beta=beta,
                Sigma2k=Sigma2k,
                n=n,
                method = optim_method)
    m = opt$par[1:n]
    v = exp(opt$par[(n+1):(2*n)])

    if(est_beta){
      beta = sum(qz/Sigma2k*matrix(m,nrow=n,ncol=K,byrow=F)/sum(qz/Sigma2k))
    }

    w = colMeans(qz)
    w = pmax(w, 1e-15)

    #print(w)
    obj[iter+1] = pois_mean_GMG_obj(x,s,m,v,w,beta,Sigma2k,qz)
    if((obj[iter+1] - obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }


  }

  return(list(posterior = list(mean_log = m,
                               var_log = v,
                               mean = exp(m+v/2)),
              fitted_g = list(weight = w,mean=beta,var=sigma2k),
              obj_value=obj,
              qz=qz))
  #return(list(posteriorMean=m,posteriorVar=v,obj_value=obj,w=w,qz=qz,priorMean=beta))

}
#'calculate objective function
pois_mean_GMG_obj = function(x,s,m,s2,w,beta,Sigma2k,qz){
  n = length(x)
  K = length(w)
  W = matrix(w,nrow=n,ncol=K,byrow=T)
  M = matrix(m,nrow=n,ncol=K,byrow=F)
  S2 = matrix(s2,nrow=n,ncol=K,byrow=F)
  return(sum(x*m-s*exp(m+s2/2))+sum(qz*(log(W)-log(Sigma2k)/2-(M^2+S2-2*M*beta+beta^2)/2/Sigma2k))-sum(qz*log(qz))+sum(log(s2))/2)
}

#' calc obj for optim
pois_mean_GMG_opt_obj = function(theta,x,s,qz,beta,Sigma2k,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-s*exp(m+exp(v)/2)-sum((m^2+exp(v)-2*m*beta)*c(rowSums(qz/Sigma2k)))/2+v/2))
}

#'calculate gradient vector for optim
pois_mean_GMG_opt_obj_gradient = function(theta,x,s,qz,beta,Sigma2k,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+exp(v)/2)-(m-beta)*c(rowSums(qz/Sigma2k)))
  g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2*c(rowSums(qz/Sigma2k)) + 1/2)
  return(c(g1,g2))
}


#' pois_mean_GMG1 = function(x,
#'                           s,
#'                           qz,
#'                           beta,
#'                           sigma2k,
#'                           optim_method,
#'                           m_init  = NULL,
#'                           s2_init = NULL){
#'   # init m, s2
#'   if(is.null(m_init)){
#'     m = 0
#'   }else{
#'     m = m_init
#'   }
#'
#'   if(is.null(s2_init)){
#'     s2 = 1
#'   }else{
#'     s2 = s2_init
#'   }
#'
#'   opt = optim(c(m,log(s2)),
#'               fn = pois_mean_GMG1_obj,
#'               gr = pois_mean_GMG1_obj_gradient,
#'               x=x,
#'               s=s,
#'               qz=qz,
#'               beta=beta,
#'               sigma2k=sigma2k,
#'               method = optim_method)
#'
#'   return(list(m=opt$par[1],s2=exp(opt$par[2]),obj=-opt$value))
#' }
#'
#'
#' #'calculate objective function
#' pois_mean_GMG1_obj = function(theta,x,s,qz,beta,sigma2k){
#'   return(-(x*theta[1]-s*exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2*sum(qz/sigma2k)+log(exp(theta[2]))/2))
#' }
#'
#' #'calculate gradient vector
#' pois_mean_GMG1_obj_gradient = function(theta,x,s,qz,beta,sigma2k){
#'   g1 = -(x-s*exp(theta[1]+exp(theta[2])/2)-theta[1]*sum(qz/sigma2k)+beta*sum(qz/sigma2k))
#'   g2 = -(-exp(theta[2])/2*s*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2*sum(qz/sigma2k) + 1/2)
#'   return(c(g1,g2))
#' }
#'
#'


