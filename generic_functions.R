# generic functions 
# valid for both, local constant and local linear

# kernel function
Wkern <- function(u,kern=c("normal","triweight"),threshold_dnorm=4){
  if (kern=="normal"){
    # Gaussian kernel
    res <- numeric(length(u))
    Iu <- (abs(u)<threshold_dnorm)
    res[Iu] <- dnorm(u[Iu])
  }else{
    ## rescaled triweight kernel
    res <- numeric(length(u))
    Iu <- (abs(u)<3)
    res[Iu] <- (1-u[Iu]^2/9)^3
  }
  ## Rectangular kernel
  #as.numeric(abs(u)<1)
  res
}

# mv: if TRUE  par1: log(mu/(1-mu)),    par2: log(var)
#     if FLASE par1: log(alpha), par2: log(beta)
my_dbeta<-function(x, par1, par2, mv=FALSE, log=TRUE){
  if (mv){# mu: par1,  var: par2
    mu  <- exp(par1)/(1+exp(par1))
    var <- exp(par2)
    aux_denom <- mu*(1-mu)/var - 1
    shape1 <- mu*aux_denom   # alpha
    shape2 <- aux_denom - shape1 # beta
  }else{
    shape1 <- exp(par1)
    shape2 <- exp(par2)
  }
  res <- lgamma(shape1+shape2)-
    lgamma(shape1)-lgamma(shape2)+
    (shape1-1)*log(x)+(shape2-1)*log(1-x)
  if (!log) res <- exp(res)
  res
}

# log_lik_Beta_val_set: Computes the mean value of the
#       log-likelihood function in the validation set
#           (tval,yval)
#       of a Functional Beta model for which
#       the functions 
#          delta(t)=log(alpha(t))
#          eta(t)=log(beta(t))
#       have been evaluated only at values t in the array newt
# 
# Output: meanloglik_val, the mean value of the
#       log-likelihood function in the validation set
log_lik_Beta_val_set <- function(tval,yval,
                                 newt, deltant, etant){
  # we interpolate the results (deltant, etant)
  # to obtain an approximation of (delta_tval, eta_tval)
  delta_tval <- approx(newt, deltant, xout=tval, rule = 2)$y
  eta_tval   <- approx(newt, etant,   xout=tval, rule = 2)$y
  alpha_tval <- exp(delta_tval) 
  beta_tval  <- exp(eta_tval)
  loglik_tval <- dbeta(yval,alpha_tval,beta_tval,log=TRUE)  
  meanloglik_val <- mean(loglik_tval)
  return(meanloglik_val)
}
