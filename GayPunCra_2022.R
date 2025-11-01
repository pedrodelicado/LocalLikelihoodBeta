## Implementing the proposal of Gaynanova, Punjabi and Crainiceanu (2022)

GayPunCra_2022 <- function(
    lom # list_of_matrices obtained from function list_2_matrix()
    , pve = 0.99 # proportion of variance explained in FPCA (using refund::fpca.face)
    , knots = 35 # number of knots to use in FPCA (using refund::fpca.face) 
    , do_plots =FALSE
){
  n <- length(lom)
  time <- lom[[1]]$time
  # Step 1: moment estimation of \mu_i(t) and \sigma^2_i(t)
  mu_tilde <- t( sapply(lom, function(x){apply(x$CGM,2,mean)}) )
  sigma_tilde_0 <- t( sapply(lom, function(x){apply(x$CGM,2,sd)}) )
  
  # Stem 2: FPCA via FACE for mu_tilde and sigma_tilde
  require(refund)
  fpca_mu <- fpca.face(mu_tilde, argvals = time, pve=pve, knots=knots)
  mu_hat <- fpca_mu$Yhat
  if (do_plots){
    matplot(time, t(mu_tilde),type="l",col=8)
    matlines(time, t(mu_hat), col=1)
    lines(time,fpca_mu$mu,col="red",lwd=4)
    title(main=paste0("FPCA for mu_tilde_i(t). Number of PCs: ",fpca_mu$npc))
  }
  sigma_tilde <- t(sapply(1:n, 
                          function(i,lom,mu_hat){
                            ni <- dim(lom[[i]]$CGM)[1]
                            sqrt( 
                              apply((lom[[i]]$CGM-rep(1,ni) %*% t(mu_hat[i,]))^2,2,sum)/(ni-1)
                            )
                          }, lom=lom, mu_hat = fpca_mu$Yhat
  ) 
  )
  fpca_sigma <- fpca.face(sigma_tilde, argvals = time, pve=pve, knots=knots)
  sigma_hat <- fpca_sigma$Yhat
  if (do_plots){
    matplot(time, t(sigma_tilde),type="l",col=8)
    matlines(time, t(sigma_hat), col=1)
    lines(time,fpca_sigma$mu,col="red",lwd=4)
    title(main=paste0("FPCA for sigma_tilde_i(t). Number of PCs: ",fpca_sigma$npc))
  }
  if (do_plots){
    matplot(time, fpca_mu$efunctions[,1:min(4,fpca_mu$npc)], type="l", 
            main=paste("fpca_mu: First", min(4,fpca_mu$npc)," principal functions"))
    abline(h=0,col=8)
    matplot(time, fpca_sigma$efunctions[,1:min(4,fpca_sigma$npc)], type="l", 
            main=paste("fpca_sigma: First", min(4,fpca_mu$npc)," principal functions"))
    abline(h=0,col=8)
    mat_aux <- cbind(fpca_mu$scores[,1:min(4,fpca_mu$npc)], 
                     fpca_sigma$scores[,1:min(4,fpca_sigma$npc)])
    colnames(mat_aux) <- c( paste("fpca_mu",1:min(4,fpca_mu$npc),sep="_"),
                            paste("fpca_sigma",1:min(4,fpca_sigma$npc),sep="_"))
    pairs(mat_aux)
    print(round(cor(mat_aux)[1:min(4,fpca_mu$npc),
                             (min(4,fpca_mu$npc)+1):(min(4,fpca_mu$npc)+min(4,fpca_sigma$npc))],2))
    print(round(cor(mat_aux)[1:min(4,fpca_mu$npc),
                             (min(4,fpca_mu$npc)+1):(min(4,fpca_mu$npc)+min(4,fpca_sigma$npc))]^2,2))
    image(1:dim(mat_aux)[2],1:dim(mat_aux)[2],cor(mat_aux)^2)
  }
  
  # Step 3: Beta estimation
  # Subject specific support: mi, Mi
  cm <- .99
  cM <- 1.01
  mi <- cm * as.numeric(lapply(lom,function(x){min(x$CGM)}))
  Mi <- cM * as.numeric(lapply(lom,function(x){max(x$CGM)}))
  if (do_plots){
    plot(Mi,ylim=c(0,410), ylab="mi and Mi",main="Min and Max of CGM values for each patient")
    points(mi, col=2)
  }
  J <- length(time)
  alpha_plus_beta <- ( mu_hat - (mi%*%t(rep(1,J))) )*( (Mi%*%t(rep(1,J))) - mu_hat )/
    sigma_hat^2 - 1
  alpha <- alpha_plus_beta * ( mu_hat - (mi%*%t(rep(1,J))) ) / ((Mi-mi)%*%t(rep(1,J))) 
  beta <- alpha_plus_beta - alpha 
  
  # quantiles
  q_025 <- (mi%*%t(rep(1,J))) + 
    ((Mi-mi)%*%t(rep(1,J))) * qbeta(.025, shape1 = alpha, shape2=beta)
  q_5 <- (mi%*%t(rep(1,J))) + 
    ((Mi-mi)%*%t(rep(1,J))) * qbeta(.5, shape1 = alpha, shape2=beta)
  q_975 <- (mi%*%t(rep(1,J))) + 
    ((Mi-mi)%*%t(rep(1,J))) * qbeta(.975, shape1 = alpha, shape2=beta)
  if (do_plots){
    matplot( time, t(q_025), type="l",col=2,lty=1,ylab="quantiles",ylim=c(min(mi),max(Mi)))
    matlines(time, t(q_975), type="l",col=4,lty=1)    
    matlines(time, t(q_5)  , type="l",col=1,lty=1)
    title(main="Quantiles 0.025, 0.5 and 0.975 for each patient and each time t")
  }
  
  if (do_plots){
    set.seed(123456)
    ng <- 6
    sample_Patients <- sort(sample(n,ng))
    mm <- min(mi) 
    MM <- max(Mi)
    op <- par(mfrow=c(2,3))
    for (i in 1:ng){
      si <- sample_Patients[i]
      matplot(time, t(lom[[si]]$CGM), type="l", col=8, ylim=c(mm,MM), 
              main=paste("si=",si,", Patient",names(lom)[si]))
      lines(time, q_025[si,], col=2, lwd=2)
      lines(time, q_975[si,], col=4, lwd=2)
      lines(time, q_5[si,]  , col=1, lwd=2)
      lines(time, mu_hat[si,],col="darkgreen", lwd=3, lty=2)
    }
    par(op)
  }
  
  # return:
  return(list(
    mu_tilde=mu_tilde,
    fpca_mu=fpca_mu,
    mu_hat=mu_hat,
    sigma_tilde_0=sigma_tilde_0,
    sigma_tilde=sigma_tilde,
    fpca_sigma=fpca_sigma,
    sigma_hat=sigma_hat,
    mi=mi,
    Mi=Mi,
    alpha=alpha,
    beta=beta,
    q_025=q_025,
    q_5=q_5,
    q_975=q_975
  )
  )
}
