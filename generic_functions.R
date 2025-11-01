# generic functions 
# valid for both, local constant and local linear

# kernel function
Wkern <- function(u,kern=c("normal","triweight")){
  if (kern=="normal"){
    # Gaussian kernel
    res <- dnorm(u)
  }else{
    ## rescaled triweight kernel
    res <- numeric(length(u))
    res[abs(u)<3] <- (1-u[abs(u)<3]^2/9)^3
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