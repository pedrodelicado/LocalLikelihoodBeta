# functions for local linear Beta regression
# 

# You also need to read this file:
# source("generic_functions.R")

## Evaluate the local likelihood with 4 parameters: 
# par=c(a0, a1, b0, b1)
# t0: point around which we calculate local likelihood
# t, y vectors of length m with the observed data
# wt vector of length m with the weights of (t[i],y[i])
#    when evaluating local log likelihood at t0
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
local_log_lik <- function(par, t0, t, y, wt, myBeta=FALSE, mv=FALSE){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0 <- t[Iwt]-t0 # vector of length equal to length(t)
  # delta <- a0+a1*t_t0 # vector of length equal to length(t)
  # eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  # alpha <- exp(a0+a1*t_t0) # vector of length equal to length(t)
  # beta  <- exp(b0+b1*t_t0) # vector of length equal to length(t)
  if (myBeta){
    log_fy <- my_dbeta(y[Iwt], par1=a0+a1*t_t0, par2=b0+b1*t_t0, mv=mv, log=TRUE)
  } else {
    log_fy <- dbeta(y[Iwt], shape1=exp(a0+a1*t_t0), shape2=exp(b0+b1*t_t0), log=TRUE)
  } 
  sum(wt[Iwt]*log_fy) 
}

## Evaluate MINUS the local likelihood with 4 parameters: 
#     par=c(a0, a1, b0, b1) 
#  It also computes the gradient and Hessian, which
#  are attached to the result as attributes, to allow
#  the use of function "nlm" to do minimization
#
# t0: point around which we calculate local likelihood
# t, y vectors of length m with the observed data
# wt vector of length m with the weights of (t[i],y[i])
#    when evaluating local log likelihood at t0
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
minus_local_log_lik_grad_hess <- function(par, t0, t, y, wt, 
                                          myBeta=FALSE, mv=FALSE){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0 <- t[Iwt]-t0 # vector of length equal to length(t)
  # delta <- a0+a1*t_t0 # vector of length equal to length(t)
  # eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  alpha <- exp(a0+a1*t_t0) # vector of length equal to length(t)
  beta  <- exp(b0+b1*t_t0) # vector of length equal to length(t)
  if (myBeta){
    log_fy <- my_dbeta(y[Iwt],shape1 = alpha, shape2 = beta, log=TRUE)
  } else {
    log_fy <- dbeta(y[Iwt],shape1 = alpha, shape2 = beta, log=TRUE)
  } 

    res <- -sum(wt[Iwt]*log_fy) 
  
  # gradient
  d_log_fy_delta <- (digamma(alpha+beta)-digamma(alpha)+log(y[Iwt]))*alpha
  d_log_fy_eta <- (digamma(alpha+beta)-digamma(beta)+log(1-y[Iwt]))*beta
  
  d_local_log_lik_a0 <- sum(wt[Iwt]*d_log_fy_delta)
  d_local_log_lik_a1 <- sum(wt[Iwt]*d_log_fy_delta*t_t0)
  
  d_local_log_lik_b0 <- sum(wt[Iwt]*d_log_fy_eta)
  d_local_log_lik_b1 <- sum(wt[Iwt]*d_log_fy_eta*t_t0)
  
  grad <- -c(d_local_log_lik_a0,d_local_log_lik_a1,
           d_local_log_lik_b0,d_local_log_lik_b1)
  
  attr(res, "gradient") <- grad
  
  # Hessian
  d2_log_fy_delta_eta <- 
    trigamma(alpha+beta)*alpha*beta
  d2_log_fy_delta2 <- 
    (trigamma(alpha+beta)-trigamma(alpha))*alpha^2 + d_log_fy_delta
  d2_log_fy_eta2 <- 
    (trigamma(alpha+beta)-trigamma(beta))*beta^2 + d_log_fy_eta
  
  d2_local_log_lik_a0a0 <- sum(wt[Iwt]*d2_log_fy_delta2)
  d2_local_log_lik_a0a1 <- sum(wt[Iwt]*d2_log_fy_delta2*t_t0)
  d2_local_log_lik_a0b0 <- sum(wt[Iwt]*d2_log_fy_delta_eta)
  d2_local_log_lik_a0b1 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0)
  
  d2_local_log_lik_a1a1 <- sum(wt[Iwt]*d2_log_fy_delta2*t_t0^2)
  d2_local_log_lik_a1b0 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0)
  d2_local_log_lik_a1b1 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0^2)
  
  d2_local_log_lik_b0b0 <- sum(wt[Iwt]*d2_log_fy_eta2)
  d2_local_log_lik_b0b1 <- sum(wt[Iwt]*d2_log_fy_eta2*t_t0)
  
  d2_local_log_lik_b1b1 <- sum(wt[Iwt]*d2_log_fy_eta2*t_t0^2)
  
  Hess <- - matrix(c(
    d2_local_log_lik_a0a0,d2_local_log_lik_a0a1,d2_local_log_lik_a0b0, d2_local_log_lik_a0b1,
    d2_local_log_lik_a0a1,d2_local_log_lik_a1a1,d2_local_log_lik_a1b0, d2_local_log_lik_a1b1,
    d2_local_log_lik_a0b0,d2_local_log_lik_a1b0,d2_local_log_lik_b0b0, d2_local_log_lik_b0b1,
    d2_local_log_lik_a0b1,d2_local_log_lik_a1b1,d2_local_log_lik_b0b1, d2_local_log_lik_b1b1
  ), ncol=4)
  
  attr(res, "hessian") <- Hess
  
  return(res)
}

# Gradient of local_log_lik with respect to its 4 parameters
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
grad_local_log_lik <- function(par, t0, t, y, wt, myBeta=FALSE, mv=FALSE){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0 <- t[Iwt]-t0 # vector of length equal to length(t)
  delta <- a0+a1*t_t0 # vector of length equal to length(t)
  eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  alpha <- exp(delta) # vector of length equal to length(t)
  beta  <- exp(eta) # vector of length equal to length(t)
  
  #log_fy <- dbeta(y,shape1 = alpha, shape2 = beta, log=TRUE)
  #sum(wt*log_fy) 
  
  d_log_fy_delta <- (digamma(alpha+beta)-digamma(alpha)+log(y[Iwt]))*alpha
  d_log_fy_eta <- (digamma(alpha+beta)-digamma(beta)+log(1-y[Iwt]))*beta
  
  d_local_log_lik_a0 <- sum(wt[Iwt]*d_log_fy_delta)
  d_local_log_lik_a1 <- sum(wt[Iwt]*d_log_fy_delta*t_t0)
  
  d_local_log_lik_b0 <- sum(wt[Iwt]*d_log_fy_eta)
  d_local_log_lik_b1 <- sum(wt[Iwt]*d_log_fy_eta*t_t0)
  
  return(c(d_local_log_lik_a0,d_local_log_lik_a1,
           d_local_log_lik_b0,d_local_log_lik_b1))
}

# Matrix form of Gradient of local_log_lik with respect to its 4 parameters
matrix_grad_local_log_lik <- function(par, t0, t, y, wt){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0 <- t[Iwt]-t0 # vector of length equal to length(t)
  m <- length(t_t0)
  
  X <- cbind(rep(1,m),t_t0)
  W <- diag(wt[Iwt])
  
  delta <- a0+a1*t_t0 # vector of length equal to length(t)
  eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  alpha <- exp(delta) # vector of length equal to length(t)
  beta  <- exp(eta) # vector of length equal to length(t)
  
  #log_fy <- dbeta(y,shape1 = alpha, shape2 = beta, log=TRUE)
  #sum(wt*log_fy)
  
  d_log_fy_delta <- (digamma(alpha+beta)-digamma(alpha)+log(y[Iwt]))*alpha
  d_log_fy_eta <- (digamma(alpha+beta)-digamma(beta)+log(1-y[Iwt]))*beta
  
  return(c( t(X) %*% W %*% d_log_fy_delta, 
            t(X) %*% W %*% d_log_fy_eta ))
}

# Hessian matrix of local_log_lik with respect to its 4 parameters
Hessian_local_log_lik <- function(par, t0, t, y, wt){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0  <- t[Iwt]-t0  # vector of length equal to length(t)
  delta <- a0+a1*t_t0 # vector of length equal to length(t)
  eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  alpha <- exp(delta) # vector of length equal to length(t)
  beta  <- exp(eta)   # vector of length equal to length(t)
  
  #log_fy <- dbeta(y,shape1 = alpha, shape2 = beta, log=TRUE)
  #sum(wt*log_fy)
  
  d_log_fy_delta <- (digamma(alpha+beta)-digamma(alpha)+log(y[Iwt]))*alpha
  d_log_fy_eta <- (digamma(alpha+beta)-digamma(beta)+log(1-y[Iwt]))*beta
  
  d2_log_fy_delta_eta <- 
    trigamma(alpha+beta)*alpha*beta
  d2_log_fy_delta2 <- 
    (trigamma(alpha+beta)-trigamma(alpha))*alpha^2 + d_log_fy_delta
  d2_log_fy_eta2 <- 
    (trigamma(alpha+beta)-trigamma(beta))*beta^2 + d_log_fy_eta
  
  d2_local_log_lik_a0a0 <- sum(wt[Iwt]*d2_log_fy_delta2)
  d2_local_log_lik_a0a1 <- sum(wt[Iwt]*d2_log_fy_delta2*t_t0)
  d2_local_log_lik_a0b0 <- sum(wt[Iwt]*d2_log_fy_delta_eta)
  d2_local_log_lik_a0b1 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0)

  d2_local_log_lik_a1a1 <- sum(wt[Iwt]*d2_log_fy_delta2*t_t0^2)
  d2_local_log_lik_a1b0 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0)
  d2_local_log_lik_a1b1 <- sum(wt[Iwt]*d2_log_fy_delta_eta*t_t0^2)
  
  d2_local_log_lik_b0b0 <- sum(wt[Iwt]*d2_log_fy_eta2)
  d2_local_log_lik_b0b1 <- sum(wt[Iwt]*d2_log_fy_eta2*t_t0)
  
  d2_local_log_lik_b1b1 <- sum(wt[Iwt]*d2_log_fy_eta2*t_t0^2)
  
  J <- matrix(c(
      d2_local_log_lik_a0a0,d2_local_log_lik_a0a1,d2_local_log_lik_a0b0, d2_local_log_lik_a0b1,
      d2_local_log_lik_a0a1,d2_local_log_lik_a1a1,d2_local_log_lik_a1b0, d2_local_log_lik_a1b1,
      d2_local_log_lik_a0b0,d2_local_log_lik_a1b0,d2_local_log_lik_b0b0, d2_local_log_lik_b0b1,
      d2_local_log_lik_a0b1,d2_local_log_lik_a1b1,d2_local_log_lik_b0b1, d2_local_log_lik_b1b1
                ), ncol=4)
  return(J)
}

# Matrix form of Hessian matrix of local_log_lik with respect to its 4 parameters
matrix_Hessian_local_log_lik <- function(par, t0, t, y, wt){
  a0 <- par[1]
  a1 <- par[2]
  b0 <- par[3]
  b1 <- par[4]
  Iwt <- which(wt>0)
  t_t0 <- t[Iwt]-t0 # vector of length equal to length(t)
  m <- length(t_t0)
  
  X <- cbind(rep(1,m),t_t0)
  W <- diag(wt[Iwt])
  
  delta <- a0+a1*t_t0 # vector of length equal to length(t)
  eta   <- b0+b1*t_t0 # vector of length equal to length(t)
  alpha <- exp(delta) # vector of length equal to length(t)
  beta  <- exp(eta) # vector of length equal to length(t)
  
  #log_fy <- dbeta(y,shape1 = alpha, shape2 = beta, log=TRUE)
  #sum(wt*log_fy)
  
  d_log_fy_delta <- (digamma(alpha+beta)-digamma(alpha)+log(y[Iwt]))*alpha
  d_log_fy_eta <- (digamma(alpha+beta)-digamma(beta)+log(1-y[Iwt]))*beta
  
  d2_log_fy_delta_eta <- 
    trigamma(alpha+beta)*alpha*beta
  d2_log_fy_delta2 <- 
    (trigamma(alpha+beta)-trigamma(alpha))*alpha^2 + d_log_fy_delta
  d2_log_fy_eta2 <- 
    (trigamma(alpha+beta)-trigamma(beta))*beta^2 + d_log_fy_eta
  
  V_delta_2 <- diag(d2_log_fy_delta2)
  V_delta_eta <- diag(d2_log_fy_delta_eta)
  V_eta_2 <- diag(d2_log_fy_eta2)
  
  J <-   rbind(
    cbind( t(X) %*% W %*% V_delta_2 %*% X, t(X) %*% W %*% V_delta_eta %*% X ),   
    cbind( t(X) %*% W %*% V_delta_eta %*% X, t(X) %*% W %*% V_eta_2 %*% X )   
  )
  
  return(J)
}


### Estimation of Beta log(parameters) by maximum local log likelihood
#   We use local linear approximations to 
#   delta(t)=log(alpha(t)) and eta(t)=log(beta(t)).
#   So we use local likelihood with 4 parameters: (a0,a1, b0,b1)
# t, y vectors of length m with the observed data
# newt: vector with points at which we estimate the Beta log(parameters)
# h: bandwidth to compute weights using a Gaussian kernel
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
loc_lik_Beta <- function(t,y,h, 
                         newt=seq(min(t),max(t),length=101),
                         myBeta=FALSE, mv = FALSE, kern="normal", 
                         method="Nelder-Mead", gr=NULL){
  lt <- length(t)
  lnt <- length(newt)
  log_lik <- deltat <- etat <- numeric(lnt)
  ab <- matrix(nrow = lnt, ncol=4)

  for (j in (1:lnt)){
    wt <- Wkern((t-newt[j])/h, kern=kern) # We follow Loader 1999, eq (2.2)
    #    wt <- wt/sum(wt)
    Iwt <- which(wt>0)

    # starting values
    # if (j==1){
    wmy <- sum(y[Iwt]*wt[Iwt])/sum(wt[Iwt])
    wmy2 <- sum((y[Iwt]^2)*wt[Iwt])/sum(wt[Iwt])
    wvy <- wmy2-wmy^2
    if (mv){
      par0 <- c(log(wmy/(1-wmy)),0,log(wvy),0)
    }else{
      aux_den <- wmy*(1-wmy)/wvy - 1
      alpha0 <- wmy*aux_den
      beta0 <- aux_den - alpha0
      # for numeric stability, we avoid too large parameter values 
      par0 <- c(min(log(alpha0),8), 0, min(log(beta0),8), 0)
    }
    # }else{
    #   par0 <- res_opt$par
    # }
    
    res_opt <- optim(par0, method=method, gr=gr, 
                     fn=local_log_lik, myBeta=myBeta, mv=mv,
                     t0=newt[j], t=t[Iwt], y=y[Iwt], wt=wt[Iwt],
                     control=list(trace=0, fnscale=-1))
    
    deltat[j] <- res_opt$par[1]
    etat[j] <- res_opt$par[3]
    ab[j,] <- res_opt$par
  }
  return(list(newt=newt, deltat=deltat, 
              etat=etat, ab=ab))
}

### Estimation of Beta log(parameters) by 
#   maximum local log likelihood.
#
#   The computations for leave-one-out CV are also done
#
#   We use local linear approximations to 
#   delta(t)=log(alpha(t)) and eta(t)=log(beta(t)).
#   So we use local likelihood with 4 parameters: (a0,a1, b0,b1)
# t, y vectors of length m with the observed data
# newt: vector with points at which we estimate the Beta log(parameters)
# h: bandwidth to compute weights using a Gaussian kernel
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
loc_lik_Beta_with_loo <- function(t,y,h, newt=t,
                         method="Nelder-Mead", gr=NULL, 
                         myBeta=FALSE, mv = FALSE, kern="normal"){
  lt <- length(t)
  lnt <- length(newt)
  loo_decrement <- loologlik <- loglik <- deltat <- etat <- numeric(lnt)
  infl <- ab <- matrix(nrow = lnt, ncol=4)

  for (j in (1:lnt)){
    wt <- Wkern((t-newt[j])/h, kern=kern) # We follow Loader 1999, eq (2.2)
    #    wt <- wt/sum(wt)
    Iwt <- which(wt>0)

    # starting values
    # if (j==1){
      wmy <- sum(y[Iwt]*wt[Iwt])/sum(wt[Iwt])
      wmy2 <- sum((y[Iwt]^2)*wt[Iwt])/sum(wt[Iwt])
      wvy <- wmy2-wmy^2
      aux_den <- wmy*(1-wmy)/wvy - 1
      alpha0 <- wmy*aux_den
      beta0 <- aux_den - alpha0
      # for numeric stability, we avoid too large parameter values 
      par0 <- c(min(log(alpha0),8), 0, min(log(beta0),8), 0)
      # }else{
    #   par0 <- res_opt$par
    # }
    
    if (method!="nlm"){
      res_opt <- optim(par0, method=method, gr=gr, 
                       fn=local_log_lik, myBeta=myBeta, mv=mv,
                       t0=newt[j], t=t[Iwt], y=y[Iwt], wt=wt[Iwt],
                       control=list(trace=0, fnscale=-1))
    }else{
      res_opt <- nlm(f=minus_local_log_lik_grad_hess,
                     p=par0, myBeta=myBeta, mv=mv,
                     t0=newt[j], t=t[Iwt], y=y[Iwt], wt=wt[Iwt])
      res_opt$par <- res_opt$estimate
    }
    deltat[j] <- res_opt$par[1]
    etat[j] <- res_opt$par[3]
    ab[j,] <- res_opt$par
    
    # abt <- loclikBeta_t$ab
    # deltat <- loclikBeta_t$deltat
    # etat <- loclikBeta_t$etat
    
    alphat <- exp(deltat[j]) 
    betat <- exp(etat[j])
    d_log_fy_deltat <- (digamma(alphat+betat)-digamma(alphat)+log(y[j]))*alphat
    d_log_fy_etat <- (digamma(alphat+betat)-digamma(betat)+log(1-y[j]))*betat
    grad_log_fy <- cbind(d_log_fy_deltat,d_log_fy_etat)
    
    Jj <- -matrix_Hessian_local_log_lik(par=ab[j,], t[j], t, y, wt)
    infl_j <- solve(Jj)[c(1,3),c(1,3)]*Wkern(0, kern=kern) # We follow Loader 1999, eq (2.28)
    infl[j,] <- infl_j
      
    loo_decrement[j] <- grad_log_fy %*% infl_j %*% t(grad_log_fy)
    if (myBeta){
      loglik[j] <- my_dbeta(y[j],alphat,betat,log=TRUE)  
    } else {
      loglik[j] <- dbeta(y[j],alphat,betat,log=TRUE)  
    } 
    loologlik[j] <- loglik[j] - loo_decrement[j]
  }
  sumloglik <- sum(loglik)
  sumloologlik <- sum(loologlik)
  return(list(newt=newt, 
              deltat=deltat, etat=etat, ab=ab, 
              loglik=loglik, loologlik=loologlik, 
              sumloglik=sumloglik, sumloologlik=sumloologlik,
              loo_decrement=loo_decrement, infl=infl))
}


### Estimation of Beta log(parameters) by maximum local log likelihood
#   using the minimization function "nlm", 
#   which uses a Newton-type algorithm, using gradient and Hessian. 
#
#   We use local linear approximations to 
#   delta(t)=log(alpha(t)) and eta(t)=log(beta(t)).
#   So we use local likelihood with 4 parameters: (a0,a1, b0,b1)
#
# t, y vectors of length m with the observed data
# newt: vector with points at which we estimate the Beta log(parameters)
# h: bandwidth to compute weights using a Gaussian kernel
loc_lik_Beta_nlm <- function(t,y,h,newt=seq(min(t),max(t),length=101), 
                             myBeta=FALSE, mv=FALSE, kern="normal"){
  lt <- length(t)
  lnt <- length(newt)
  log_lik <- deltat <- etat <- numeric(lnt)
  ab <- matrix(nrow = lnt, ncol=4)

  for (j in (1:lnt)){
    wt <- Wkern((t-newt[j])/h, kern=kern) # We follow Loader 1999, eq (2.2)
    #    wt <- wt/sum(wt)
    Iwt <- which(wt>0)
    
    # starting values
    # if (j==1){
      wmy <- sum(y[Iwt]*wt[Iwt])/sum(wt[Iwt])
      wmy2 <- sum((y[Iwt]^2)*wt[Iwt])/sum(wt[Iwt])
      wvy <- wmy2-wmy^2
      aux_den <- wmy*(1-wmy)/wvy - 1
      alpha0 <- wmy*aux_den
      beta0 <- aux_den - alpha0
      # for numeric stability, we avoid too large parameter values 
      par0 <- c(min(log(alpha0),8), 0, min(log(beta0),8), 0)
    # }else{
    #   par0 <- res_nlm$estimate
    # }
    
    res_nlm <- nlm(f=minus_local_log_lik_grad_hess,
                     p=par0, myBeta=myBeta, mv=mv, 
                     t0=newt[j], t=t[Iwt], y=y[Iwt], wt=wt[Iwt])
    
    deltat[j] <- res_nlm$estimate[1]
    etat[j] <- res_nlm$estimate[3]
    ab[j,] <- res_nlm$estimate
  }
  return(list(newt=newt, deltat=deltat, 
              etat=etat, ab=ab))
}



# Choosing the bandwidth h by leave-one-out cross-validation.
#    The leave-one-out process is explicitly done.
#    So this function is too slow.
#
# t, y vectors of length m with the observed data
# vh: vector of possible bandwidth values h. 
#     By default: vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh)
# lh: length of the default value for vh
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
looCV_loclikBeta <- function(t, y, lh=10,
                             vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh),
                             method="Nelder-Mead", gr=NULL, 
                             myBeta=FALSE, mv=FALSE, kern="normal"){
  lvh <- length(vh)
  lt <- length(t)
  loologlik <- numeric(lvh)
  for (k in (1:lvh)){
    # print(paste(k,"of",lvh))
    h <- vh[k]
    for (j in (1:lt)){
      if (method!="nlm"){
        loclikBeta_t0 <- loc_lik_Beta(t[-j],y[-j], h, newt = t[j],
                                      method=method, gr=gr,
                                      myBeta=myBeta, mv=mv, kern=kern)
      }else{
        loclikBeta_t0 <- loc_lik_Beta_nlm(t[-j],y[-j], h, newt = t[j],
                                          myBeta=myBeta, mv=mv, kern=kern)
      }
      deltat0 <- loclikBeta_t0$deltat
      etat0 <- loclikBeta_t0$etat
      alphat0 <- exp(deltat0) 
      betat0 <- exp(etat0)
      if (myBeta){
        loologlik[k] <- loologlik[k] + my_dbeta(y[j],alphat0,betat0,log=TRUE)
      }else{
        loologlik[k] <- loologlik[k] + dbeta(y[j],alphat0,betat0,log=TRUE)
      }
    }
  }
  return(list(vh=vh,loglik=loologlik))
}


# Choosing the bandwidth h by approximate leave-one-out cross-validation.
#    The leave-one-out process is approximated
#    by the influence function.
#    (the derivation follows Problem 4.17 in Loader 1999)
#
# t, y vectors of length m with the observed data
# vh: vector of possible bandwidth values h. 
#     By default: vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh)
# lh: length of the default value for vh
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
infl_looCV_loclikBeta <- function(t, y, lh=10,
                             vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh),
                             method="Nelder-Mead", gr=NULL, 
                             myBeta=FALSE, mv=FALSE, kern="normal"){
  lvh <- length(vh)
  m <- length(t)
  loologlik <- numeric(lvh)
  loo_decrement <- matrix(nrow=lvh, ncol=m)
  infl <- array(dim=c(lvh,m,4))
  
  for (k in (1:lvh)){
    # print(paste(k,"of",lvh))
    h <- vh[k]
    if (method!="nlm"){
      loclikBeta_t <- loc_lik_Beta(t,y, h, newt = t,
                                   method=method, gr=gr,
                                   myBeta=myBeta, mv=mv, kern=kern)
    }else{
      loclikBeta_t <- loc_lik_Beta_nlm(t,y, h, newt = t,
                                       myBeta=myBeta, mv=mv, kern=kern)
    }
    abt <- loclikBeta_t$ab
    deltat <- loclikBeta_t$deltat
    etat <- loclikBeta_t$etat
    alphat <- exp(deltat) 
    betat <- exp(etat)
    d_log_fy_deltat <- (digamma(alphat+betat)-digamma(alphat)+log(y))*alphat
    d_log_fy_etat <- (digamma(alphat+betat)-digamma(betat)+log(1-y))*betat
    grad_log_fy <- cbind(d_log_fy_deltat,d_log_fy_etat)
    
    for (j in (1:m)){
      wtj <- Wkern((t-t[j])/h, kern=kern)
      Jj <- -matrix_Hessian_local_log_lik(par=abt[j,], t[j], t, y, wtj)
      infl_j <- solve(Jj)[c(1,3),c(1,3)]*Wkern(0, kern=kern) # We follow Loader 1999, eq (2.28)
      infl[k,j,] <- infl_j
      
      loo_decrement[k,j] <- t(grad_log_fy[j,]) %*% infl_j %*% grad_log_fy[j,]
      if (myBeta){
        log_dens_j <- my_dbeta(y[j],alphat[j],betat[j],log=TRUE)
      }else{
        log_dens_j <- dbeta(y[j],alphat[j],betat[j],log=TRUE)
      }
      
      loologlik[k] <- loologlik[k] + log_dens_j - loo_decrement[k,j]
    }
  }
  return(list(vh=vh,loologlik=loologlik,
              infl=infl,loo_decrement=loo_decrement))
}


# 2nd version of "infl_looCV_loclikBeta".
# Instead of calling "loc_lik_Beta", 
# now the function "loc_lik_Beta_with_loo" is called 
# to avoid a second for loop over points.
#
# Choosing the bandwidth h by approximate leave-one-out cross-validation.
#    The leave-one-out process is approximated
#    by the influence function.
#    (the derivation follows Problem 4.17 in Loader 1999)
#
# t, y vectors of length m with the observed data
# vh: vector of possible bandwidth values h. 
#     By default: vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh)
# lh: length of the default value for vh
infl_looCV_loclikBeta_2nd <- function(t, y, lh=10,
                                  vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh),
                                  method="Nelder-Mead", gr=NULL, 
                                  myBeta=FALSE, mv=FALSE, kern="normal"){
  lvh <- length(vh)
  m <- length(t)
  loologlik <- numeric(lvh)
  loo_decrement <- matrix(nrow=lvh, ncol=m)
  infl <- array(dim=c(lvh,m,4))
  
  for (k in (1:lvh)){
    # print(paste(k,"of",lvh))
    h <- vh[k]
    loclikBeta_t <- loc_lik_Beta_with_loo(t,y, h, newt = t,
                                   method=method, gr=gr,
                                   myBeta=myBeta, mv=mv, kern=kern)
    loologlik[k] <- loclikBeta_t$sumloologlik
    loo_decrement[k,] <- loclikBeta_t$loo_decrement
    infl[k,,] <- loclikBeta_t$infl
  }
  return(list(vh=vh,loologlik=loologlik,
              infl=infl,loo_decrement=loo_decrement))
}

####

# Choosing the bandwidth h by k-fold cross-validation.
#
# t, y vectors of length m with the observed data
# kCV: Number of folds, in the k-fold CV process (default: 5)
# lh: length of the default value for vh (default: 10)
# vh: vector of possible bandwidth values h. 
#     By default: vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh)
# myBeta: if TRUE use my_dbeta, otherwise
# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
k_fold_CV_loclikBeta <- function(t, y, kCV=5, 
                            lh=10,
                            vh=seq((max(t)-min(t))/20,(max(t)-min(t))/2,length=lh),
                            method="Nelder-Mead", gr=NULL, 
                            myBeta=FALSE, mv=FALSE, kern="normal"){
  lvh <- length(vh)
  lt <- length(t)
  kfCVloglik <- numeric(lvh)
  
  size_fold <- ceiling(lt/kCV)
  Ifold <- sample(rep(1:kCV,size_fold)[1:lt])

  for (k in (1:lvh)){
    # print(paste(k,"of",lvh))
    h <- vh[k]
    for (j in (1:kCV)){
      Ij <- which(Ifold==j)
      if (method!="nlm"){
        loclikBeta_t0 <- loc_lik_Beta(t[-Ij],y[-Ij], h, newt = t[Ij],
                                      method=method, gr=gr,
                                      myBeta=myBeta, mv=mv, kern=kern)
      }else{
        loclikBeta_t0 <- loc_lik_Beta_nlm(t[-Ij],y[-Ij], h, newt = t[Ij],
                                          myBeta=myBeta, mv=mv, kern=kern)
      }
      deltat0 <- loclikBeta_t0$deltat
      etat0 <- loclikBeta_t0$etat
      alphat0 <- exp(deltat0) 
      betat0 <- exp(etat0)
      if (myBeta){
        kfCVloglik[k] <- kfCVloglik[k] + sum(my_dbeta(y[Ij],alphat0,betat0,log=TRUE))
      }else{
        kfCVloglik[k] <- kfCVloglik[k] + sum(dbeta(y[Ij],alphat0,betat0,log=TRUE))
      }
    }
  }
  return(list(vh=vh,loglik=kfCVloglik))
}
