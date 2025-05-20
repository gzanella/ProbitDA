#auxiliary
beta_sampler <- function(mu, chol_Sigma, p){
  #mu = mean of beta
  #chol_Sigma = Cholesky decomposition variance
  #p = dimension of mu
  
  beta <- mu + chol_Sigma%*%matrix(rnorm(p), p, 1)
  beta <- matrix(beta, p, 1)
  return(beta)
}
z_sampler <- function(mu, lower_lim, upper_lim, n){
  #mu = mean of z
  #lower_lim = lower limit of truncations
  #upper_lim = upper limit of truncations
  #n = length of mu
  
  z <- rtruncnorm(n, lower_lim, upper_lim, as.vector(mu))
  z <- matrix(z, n, 1)
  return(z)
}
rwm_step_beta <- function(old_int, old_pred, y, var_int, sd_update, kappa, t){
  #updating intercept
  #old_pred = old value predictors
  #y = observed data
  #var_int = prior variance intercept. Prior variance assumed to be diagonal
  #sd_update = standard deviation proposal
  #kappa = learning rate
  #t = iteration
  
  #proposal
  shift <- sd_update*rnorm(1)
  prop_pred <- old_pred+shift
  
  #acceptance probability
  old_loglik <- sum(y*pnorm(old_pred, log = T)) + sum((1-y)*pnorm(old_pred, log = T, lower.tail = F))
  prop_loglik <- sum(y*pnorm(prop_pred, log = T)) + sum((1-y)*pnorm(prop_pred, log = T, lower.tail = F))
  log_prop <- prop_loglik - old_loglik-shift*old_int/var_int-shift^2/(2*var_int)
  alpha <- min(c(1, exp(log_prop)))
  
  #accepting
  new_loc <- 0
  accepted <- 0
  if(log(runif(1)) < log_prop){
    new_loc <- shift
    accepted <- 1
  }
  
  #update sd
  log_sd <- log(sd_update) + t^(-kappa)*(alpha-0.23)
  
  return(list(new_loc, exp(log_sd), accepted))
}
rwm_step_pred <- function(old_pred, y, Prec_one, one_Prec_one, sd_update, kappa, t){
  #updating intercept
  #old_pred = old value predictors
  #y = observed data
  #Prec_one = prior precision matrix pred times a vector of ones
  #one_Prec_one = vector ones times prior precision matrix pred times a vector of ones
  #sd_update = standard deviation proposal
  #kappa = learning rate
  #t = iteration
  
  #proposal
  shift <- sd_update*rnorm(1)
  prop_pred <- old_pred+shift
  
  #acceptance probability
  old_loglik <- sum(y*pnorm(old_pred, log = T)) + sum((1-y)*pnorm(old_pred, log = T, lower.tail = F))
  prop_loglik <- sum(y*pnorm(prop_pred, log = T)) + sum((1-y)*pnorm(prop_pred, log = T, lower.tail = F))
  log_prop <- prop_loglik - old_loglik-shift*t(old_pred)%*%Prec_one-shift^2/2*one_Prec_one
  alpha <- min(c(1, exp(log_prop)))
  
  #accepting
  new_pred <- old_pred
  accepted <- 0
  if(log(runif(1)) < log_prop){
    new_pred <- prop_pred
    accepted <- 1
  }
  
  #update sd
  log_sd <- log(sd_update) + t^(-kappa)*(alpha-0.23)
  
  return(list(new_pred, exp(log_sd), accepted))
}
rwm_step_z <- function(z_vec, Q_one, one_Q_one, upper_lim, lower_lim, sd_update, kappa, t){
  #updating intercept
  #z_vec = values of z (column)
  #Q_one <- Q*ones
  #one_Q_one <- t(ones)*Q*ones
  #sd_update = standard deviation proposal
  #kappa = learning rate
  #t = iteration
  
  #proposal
  prop <- sd_update*rnorm(1)
  
  #acceptance probability
  log_prop <- -prop*t(z_vec)%*%Q_one-prop^2/2*one_Q_one
  #checking the limits
  log_prop <- log_prop + ifelse(all(as.vector(z_vec)+prop < upper_lim), 0, -Inf)
  log_prop <- log_prop + ifelse(all(as.vector(z_vec)+prop > lower_lim), 0, -Inf)
  alpha <- min(c(1, exp(log_prop)))
  
  #accepting
  new_loc <- 0
  accepted <- 0
  if(log(runif(1)) < log_prop){
    new_loc <- prop
    accepted <- 1
  }
  
  #update sd
  log_sd <- log(sd_update) + t^(-kappa)*(alpha-0.23)
  
  return(c(new_loc, exp(log_sd), accepted))
}


#DA
coupling_DA_n_bigger_p <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix
  V <- chol2inv(chol(t(X)%*%X+solve(Sigma0)))
  chol_V <- t(chol(V)) #to have it lower triangular
  inv_chol_V <- solve(chol_V)
  VX_t <- V%*%t(X)
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #starting points
  beta1 <- beta_sampler(0, t(chol(Sigma0)), p)
  z1 <- z_sampler(X%*%beta1, lower_lim, upper_lim, n)
  for(i in 1:L){
    z1 <- z_sampler(X%*%beta1, lower_lim, upper_lim, n)
    beta1 <- beta_sampler(VX_t%*%z1, chol_V, p)
  }
  
  beta2 <- beta_sampler(0, t(chol(Sigma0)), p)
  z2 <- z_sampler(X%*%beta2, lower_lim, upper_lim, n)
  
  
  
  #distance
  dist2 <- sum((beta1-beta2)^2) + sum((z1-z2)^2)
  #dist2 <- sum((beta1-beta2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #z
      U <- runif(n)
      z1 <- qtruncnorm(U, lower_lim, upper_lim, X%*%beta1)
      z1 <- matrix(z1, n, 1)
      z2 <- qtruncnorm(U, lower_lim, upper_lim, X%*%beta2)
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- matrix(rnorm(p), p, 1)
      beta1 <- VX_t%*%z1 + chol_V%*%Z
      beta2 <- VX_t%*%z2 + chol_V%*%Z
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #z
      pred1 <- as.vector(X%*%beta1)
      pred2 <- as.vector(X%*%beta2)
      #sample z1
      z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
      log1 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred1)))
      log2 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred2)))
      if(log(runif(1)) <= log2-log1){
        z2 <- z1
      }
      else{
        reject <- T
        while(reject){
          #propose z2
          prop <- z_sampler(pred2, lower_lim, upper_lim, n)
          log1 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred1)))
          log2 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred2)))
          if(log(runif(1)) <= log1 - log2){
            z2 <- prop
            reject <- F
          }
        }
      }
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- inv_chol_V%*%VX_t%*%(z1-z2)
      e <- Z/sqrt(sum(Z^2))
      X_p <-  matrix(rnorm(p), p, 1)
      if(log(runif(1)) <= -0.5*t(Z)%*%(2*X_p+Z)){
        Y_p <- X_p + Z
      }
      else{
        Y_p <- X_p -2*as.vector(t(e)%*%X_p)*e
      }
      beta1 <- chol_V%*%X_p + VX_t%*%z1
      beta2 <- chol_V%*%Y_p + VX_t%*%z2
    }
    #compute distance and check
    dist2 <- sum((beta1-beta2)^2) + sum((z1-z2)^2)
    #dist2 <- sum((beta1-beta2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 10), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, beta1))
  }
  return(n_iter)
}

coupling_DA_p_bigger_n <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix
  X_Sigma_X_T <- X%*%Sigma0%*%t(X)
  Q <- chol2inv(chol(diag(n)+X_Sigma_X_T))
  W <- X_Sigma_X_T%*%(diag(n)-Q%*%X_Sigma_X_T)
  chol_W <- t(chol(W)) #to have it lower triangular
  inv_chol_W <- solve(chol_W)
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #starting points
  pred1 <- X%*%beta_sampler(0, chol(Sigma0), p)
  z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
  for(i in 1:L){
    z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
    pred1 <- beta_sampler(W%*%z1, chol_W, n)
  }
  
  pred2 <- X%*%beta_sampler(0, chol(Sigma0), p)
  z2 <- z_sampler(pred2, lower_lim, upper_lim, n)
  
  
  
  #distance
  dist2 <- sum((pred1-pred2)^2) + sum((z1-z2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #z
      U <- runif(n)
      z1 <- qtruncnorm(U, lower_lim, upper_lim, pred1)
      z1 <- matrix(z1, n, 1)
      z2 <- qtruncnorm(U, lower_lim, upper_lim, pred2)
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- matrix(rnorm(n), n, 1)
      pred1 <- W%*%z1 + chol_W%*%Z
      pred2 <- W%*%z2 + chol_W%*%Z
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #z
      #sample z1
      z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
      log1 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred1)))
      log2 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred2)))
      if(log(runif(1)) <= log2-log1){
        z2 <- z1
      }
      else{
        reject <- T
        while(reject){
          #propose z2
          prop <- z_sampler(pred2, lower_lim, upper_lim, n)
          log1 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred1)))
          log2 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred2)))
          if(log(runif(1)) <= log1 - log2){
            z2 <- prop
            reject <- F
          }
        }
      }
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- inv_chol_W%*%W%*%(z1-z2)
      e <- Z/sqrt(sum(Z^2))
      X_p <-  matrix(rnorm(n), n, 1)
      if(log(runif(1)) <= -0.5*t(Z)%*%(2*X_p+Z)){
        Y_p <- X_p + Z
      }
      else{
        Y_p <- X_p -2*as.vector(t(e)%*%X_p)*e
      }
      pred1 <- chol_W%*%X_p + W%*%z1
      pred2 <- chol_W%*%Y_p + W%*%z2
    }
    #compute distance and check
    dist2 <- sum((pred1-pred2)^2) + sum((z1-z2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 3), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, pred1))
  }
  return(n_iter)
}

#DA noncentered
coupling_DA_n_bigger_p_noncIntercept <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix
  V <- chol2inv(chol(t(X)%*%X+solve(Sigma0)))
  chol_V <- t(chol(V)) #to have it lower triangular
  inv_chol_V <- solve(chol_V)
  VX_t <- V%*%t(X)
  var_int <- Sigma0[1,1] #prior variance of intercept
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #intercept
  SD_update <- rep(0, L)
  Accepted <- rep(0, L)
  sd_update <- 2.4
  kappa <- 0.6 #power learning rate
  
  #starting points
  beta1 <- beta_sampler(0, t(chol(Sigma0)), p)
  z1 <- z_sampler(X%*%beta1, lower_lim, upper_lim, n)
  for(i in 1:L){
    old_int <- beta1[1,1]
    pred <- X%*%beta1
      
    #update intercept
    res <- rwm_step_beta(old_int, pred, y, var_int, sd_update, kappa, i)
    shift <- res[[1]]
    #diagnostic on the adaptive MCMC
    sd_update <- res[[2]]
    SD_update[i] <- sd_update
    Accepted[i] <- res[[3]]
    
    pred <- pred + shift
    z1 <- z_sampler(pred, lower_lim, upper_lim, n)
    beta1 <- beta_sampler(VX_t%*%z1, chol_V, p)
  }
  if(verbose){
    print("End of lag")
    plot(SD_update, type = "l", lwd = 2, main = "SD of proposal", xlab = "", ylab = "")
    print(paste("Proportion of acceptance of RwM:", mean(Accepted)))
  }
  
  beta2 <- beta_sampler(0, t(chol(Sigma0)), p)
  z2 <- z_sampler(X%*%beta2, lower_lim, upper_lim, n)
  
  
  
  #distance
  dist2 <- sum((beta1-beta2)^2) + sum((z1-z2)^2)
  #dist2 <- sum((beta1-beta2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    #update intercept: only maximal (reflection) coupling
    old_int1 <- beta1[1,1]
    old_int2 <- beta2[1,1]
    z <- 1/sd_update*(old_int1-old_int2)
    x_p <- rnorm(1)
    if(log(runif(1)) <= -0.5*z*(2*x_p+z)){
      y_p <- x_p + z
    }
    else{
      y_p <- -x_p
    }
    new1 <- sd_update*x_p + old_int1
    new2 <- sd_update*y_p + old_int2
    
    u <- runif(1)
    
    pred1 <- X%*%beta1
    shift1 <- new1 - old_int1
    prop_pred1 <- pred1 + shift1
    
    pred2 <- X%*%beta2
    shift2 <- new2 - old_int2
    prop_pred2 <- pred2 + shift2
    
    #acceptance probability
    old_loglik1 <- sum(y*pnorm(pred1, log = T)) + sum((1-y)*pnorm(pred1, log = T, lower.tail = F))
    prop_loglik1 <- sum(y*pnorm(prop_pred1, log = T)) + sum((1-y)*pnorm(prop_pred1, log = T, lower.tail = F))
    log_prop1 <- prop_loglik1 - old_loglik1 - shift1*old_int1/var_int-shift1^2/(2*var_int)
    
    old_loglik2 <- sum(y*pnorm(pred2, log = T)) + sum((1-y)*pnorm(pred2, log = T, lower.tail = F))
    prop_loglik2 <- sum(y*pnorm(prop_pred2, log = T)) + sum((1-y)*pnorm(prop_pred2, log = T, lower.tail = F))
    log_prop2 <- prop_loglik2 - old_loglik2 - shift2*old_int2/var_int-shift2^2/(2*var_int)
    
    #accepting
    if(log(u) < log_prop1){
      beta1[1,1] <- new1
      pred1 <- prop_pred1
    }
    if(log(u) < log_prop2){
      beta2[1,1] <- new2
      pred2 <- prop_pred2
    }
    
    dist2 <- sum((beta1-beta2)^2) + sum((z1-z2)^2)
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #z
      U <- runif(n)
      z1 <- qtruncnorm(U, lower_lim, upper_lim, pred1)
      z1 <- matrix(z1, n, 1)
      z2 <- qtruncnorm(U, lower_lim, upper_lim, pred2)
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- matrix(rnorm(p), p, 1)
      beta1 <- VX_t%*%z1 + chol_V%*%Z
      beta2 <- VX_t%*%z2 + chol_V%*%Z
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #z
      pred1 <- as.vector(X%*%beta1)
      pred2 <- as.vector(X%*%beta2)
      #sample z1
      z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
      log1 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred1)))
      log2 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred2)))
      if(log(runif(1)) <= log2-log1){
        z2 <- z1
      }
      else{
        reject <- T
        while(reject){
          #propose z2
          prop <- z_sampler(pred2, lower_lim, upper_lim, n)
          log1 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred1)))
          log2 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred2)))
          if(log(runif(1)) <= log1 - log2){
            z2 <- prop
            reject <- F
          }
        }
      }
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- inv_chol_V%*%VX_t%*%(z1-z2)
      e <- Z/sqrt(sum(Z^2))
      X_p <-  matrix(rnorm(p), p, 1)
      if(log(runif(1)) <= -0.5*t(Z)%*%(2*X_p+Z)){
        Y_p <- X_p + Z
      }
      else{
        Y_p <- X_p -2*as.vector(t(e)%*%X_p)*e
      }
      beta1 <- chol_V%*%X_p + VX_t%*%z1
      beta2 <- chol_V%*%Y_p + VX_t%*%z2
    }
    #compute distance and check
    dist2 <- sum((beta1-beta2)^2) + sum((z1-z2)^2)
    #dist2 <- sum((beta1-beta2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 10), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, beta1))
  }
  return(n_iter)
}

coupling_DA_p_bigger_n_noncIntercept <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta, assumed so that the intercept is independent of everything else
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix
  X_Sigma_X_T <- X%*%Sigma0%*%t(X)
  Q <- chol2inv(chol(diag(n)+X_Sigma_X_T))
  W <- X_Sigma_X_T%*%(diag(n)-Q%*%X_Sigma_X_T)
  chol_W <- chol(W) 
  Prec <- chol2inv(chol_W) #precision matrix
  Prec_one <- Prec%*%matrix(1, n, 1)
  one_Prec_one <- matrix(1, 1, n)%*%Prec_one
  chol_W <- t(chol_W) #to have it lower triangular
  inv_chol_W <- solve(chol_W)
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #intercept
  SD_update <- rep(0, L)
  Accepted <- rep(0, L)
  sd_update <- 2.4
  kappa <- 0.6 #power learning rate
  
  #starting points
  pred1 <- X%*%beta_sampler(0, t(chol(Sigma0)), p)
  z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
  for(i in 1:L){
    
    #update intercept
    res <- rwm_step_pred(pred1, y, Prec_one, one_Prec_one, sd_update, kappa, i)
    pred1 <- res[[1]]
    #diagnostic on the adaptive MCMC
    sd_update <- res[[2]]
    SD_update[i] <- sd_update
    Accepted[i] <- res[[3]]
    
    #update z1
    z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
    #update pred1
    pred1 <- beta_sampler(W%*%z1, chol_W, n)
  }
  
  pred2 <- X%*%beta_sampler(0, chol(Sigma0), p)
  z2 <- z_sampler(pred2, lower_lim, upper_lim, n)
  
  if(verbose){
    print("End of lag")
    plot(SD_update, type = "l", lwd = 2, main = "SD of proposal", xlab = "", ylab = "")
    print(paste("Proportion of acceptance of RwM:", mean(Accepted)))
  }
  
  
  
  #distance
  dist2 <- sum((pred1-pred2)^2) + sum((z1-z2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #update of the intercepts: always W2 coupling
      
      #proposal
      common_shift <- sd_update*rnorm(1)
      prop_pred1 <- pred1+common_shift
      prop_pred2 <- pred2+common_shift
      
      #acceptance probabilities
      old_loglik1 <- sum(y*pnorm(pred1, log = T)) + sum((1-y)*pnorm(pred1, log = T, lower.tail = F))
      prop_loglik1 <- sum(y*pnorm(prop_pred1, log = T)) + sum((1-y)*pnorm(prop_pred1, log = T, lower.tail = F))
      log_prop1 <- prop_loglik1 - old_loglik1 -common_shift*t(pred1)%*%Prec_one-common_shift^2/2*one_Prec_one
      
      old_loglik2 <- sum(y*pnorm(pred2, log = T)) + sum((1-y)*pnorm(pred2, log = T, lower.tail = F))
      prop_loglik2 <- sum(y*pnorm(prop_pred2, log = T)) + sum((1-y)*pnorm(prop_pred2, log = T, lower.tail = F))
      log_prop2 <- prop_loglik2 - old_loglik2-common_shift*t(pred2)%*%Prec_one-common_shift^2/2*one_Prec_one
      
      
      #accepting
      u <- runif(1)
      if(log(u) < log_prop1){
        pred1 <- prop_pred1
      }
      if(log(u) < log_prop2){
        pred2 <- prop_pred2
      }
      
      #z
      U <- runif(n)
      z1 <- qtruncnorm(U, lower_lim, upper_lim, pred1)
      z1 <- matrix(z1, n, 1)
      z2 <- qtruncnorm(U, lower_lim, upper_lim, pred2)
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- matrix(rnorm(n), n, 1)
      pred1 <- W%*%z1 + chol_W%*%Z
      pred2 <- W%*%z2 + chol_W%*%Z
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #z
      #sample z1
      z1 <- z_sampler(pred1, lower_lim, upper_lim, n)
      log1 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred1)))
      log2 <- sum(log(dtruncnorm(z1, lower_lim, upper_lim, pred2)))
      if(log(runif(1)) <= log2-log1){
        z2 <- z1
      }
      else{
        reject <- T
        while(reject){
          #propose z2
          prop <- z_sampler(pred2, lower_lim, upper_lim, n)
          log1 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred1)))
          log2 <- sum(log(dtruncnorm(prop, lower_lim, upper_lim, pred2)))
          if(log(runif(1)) <= log1 - log2){
            z2 <- prop
            reject <- F
          }
        }
      }
      z2 <- matrix(z2, n, 1)
      
      #beta
      Z <- inv_chol_W%*%W%*%(z1-z2)
      e <- Z/sqrt(sum(Z^2))
      X_p <-  matrix(rnorm(n), n, 1)
      if(log(runif(1)) <= -0.5*t(Z)%*%(2*X_p+Z)){
        Y_p <- X_p + Z
      }
      else{
        Y_p <- X_p -2*as.vector(t(e)%*%X_p)*e
      }
      pred1 <- chol_W%*%X_p + W%*%z1
      pred2 <- chol_W%*%Y_p + W%*%z2
    }
    #compute distance and check
    dist2 <- sum((pred1-pred2)^2) + sum((z1-z2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 3), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, pred1))
  }
  return(n_iter)
}

#GS
coupling_GS_n_bigger_p <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix of the beta
  X_t <- t(X)
  V <- chol2inv(chol(X_t%*%X+solve(Sigma0)))
  VX_t <- V%*%X_t
  #conditional_variances
  h <- rep(0,n)
  for(i in 1:n){
    h[i] <- sum(VX_t[,i]*X[i,])
  }
  cond_var <- 1/(1-h)
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #starting points
  beta1 <- beta_sampler(0, t(chol(Sigma0)), p)
  z1 <- z_sampler(X%*%beta1, lower_lim, upper_lim, n)
  #computing the vector B = VX_Tz
  B1 <- VX_t%*%matrix(z1, n, 1)
  for(t in 1:L){
    #each iteration is given by n steps
    choices <- sample(1:n, n, replace = T)
    for(l in 1:n){
      i <- choices[l]
      #compute conditional mean
      mean_i <- cond_var[i]*sum(X[i,]*B1)-h[i]*cond_var[i]*z1[i]
      #update old
      i_old1 <- i
      z_old1 <- z1[i]
      #sample z_i
      z1[i] <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
      #update B
      B1 <- B1 + VX_t[,i]*(z1[i]-z_old1)
    }
  }
  
  beta2 <- beta_sampler(0, t(chol(Sigma0)), p)
  z2 <- z_sampler(X%*%beta2, lower_lim, upper_lim, n)
  #computing the vector B = VX_Tz
  B2 <- VX_t%*%matrix(z2, n, 1)
  
  
  
  #distance
  dist2 <- sum((z1-z2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #each iteration is given by n steps
      choices <- sample(1:n, n, replace = T)
      U <- runif(n)
      for(l in 1:n){
        i <- choices[l]
        
        #z1
        #compute conditional mean
        mean_i <- cond_var[i]*sum(X[i,]*B1)-h[i]*cond_var[i]*z1[i]
        #update old
        i_old1 <- i
        z_old1 <- z1[i]
        #sample z_i
        z1[i] <- qtruncnorm(U[l], a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
        #update B
        B1 <- B1 + VX_t[,i]*(z1[i]-z_old1)
        
        #z2
        #compute conditional mean
        mean_i <- cond_var[i]*sum(X[i,]*B2)-h[i]*cond_var[i]*z2[i]
        #update old
        i_old2 <- i
        z_old2 <- z2[i]
        #sample z_i
        z2[i] <- qtruncnorm(U[l], a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
        #update B
        B2 <- B2 + VX_t[,i]*(z2[i]-z_old2)
      }
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #each iteration is given by n steps
      choices <- sample(1:n, n, replace = T)
      U <- runif(n)
      for(l in 1:n){
        i <- choices[l]
        
        #z1
        #compute conditional mean
        mean_i1 <- cond_var[i]*sum(X[i,]*B1)-h[i]*cond_var[i]*z1[i]
        #update old
        i_old1 <- i
        z_old1 <- z1[i]
        
        
        #z2
        #compute conditional mean
        mean_i2 <- cond_var[i]*sum(X[i,]*B2)-h[i]*cond_var[i]*z2[i]
        #update old
        i_old2 <- i
        z_old2 <- z2[i]
        
        #sample z1
        z1[i] <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i]))
        log1 <- log(dtruncnorm(z1[i], a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i])))
        log2 <- log(dtruncnorm(z1[i], a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i])))
        if(log(runif(1)) <= log2-log1){
          z2[i] <- z1[i]
        }
        else{
          reject <- T
          while(reject){
            #propose z2
            prop <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i]))
            log1 <- log(dtruncnorm(prop, a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i])))
            log2 <- log(dtruncnorm(prop, a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i])))
            if(log(runif(1)) <= log1 - log2){
              z2[i] <- prop
              reject <- F
            }
          }
        }
        
        #update B
        B1 <- B1 + VX_t[,i]*(z1[i]-z_old1)
        B2 <- B2 + VX_t[,i]*(z2[i]-z_old2)
        
      }
    }
    #compute distance and check
    dist2 <- sum((z1-z2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 3), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, z1))
  }
  return(n_iter)
}

coupling_GS_p_bigger_n <- function(X, y, Sigma0, threshold, L, max_iter = 500, verbose = F, return_sample = F){
  
  #Sigma0 = prior covariance of beta
  # the prior mean of beta is assumed to be the zero vector
  # threshold = value to pass from W2 to maximal coupling
  # L = lag
  # max_iter + number of iterations to run
  
  n <- nrow(X)
  p <- ncol(X)
  
  #pre-computation
  
  #covariance matrix of the beta
  #precision matrix
  Q <- chol2inv(chol(X%*%Sigma0%*%t(X)+diag(n)))
  #conditional variances
  cond_var <- as.vector(1/diag(Q))
  #lower and upper limits
  lower_lim <- ifelse(y == 1, 0, -Inf)
  upper_lim <- ifelse(y == 0, 0, Inf)
  
  #starting points
  beta1 <- beta_sampler(0, t(chol(Sigma0)), p)
  z1 <- z_sampler(X%*%beta1, lower_lim, upper_lim, n)
  z1 <- as.vector(z1)
  for(t in 1:L){
    #each iteration is given by n steps
    choices <- sample(1:n, n, replace = T)
    for(l in 1:n){
      i <- choices[l]
      #compute conditional mean
      mean_i <- -cond_var[i]*sum(Q[i,-i]*z1[-i])
      #sample z_i
      z1[i] <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
    }
  }
  
  beta2 <- beta_sampler(0, t(chol(Sigma0)), p)
  z2 <- z_sampler(X%*%beta2, lower_lim, upper_lim, n)
  z2 <- as.vector(z2)
  
  
  
  #distance
  dist2 <- sum((z1-z2)^2)
  
  #flag
  flag <- "W2"
  
  n_iter <- L
  different <- T
  while(different & n_iter < max_iter){
    n_iter <- n_iter + 1
    
    if(dist2 >= threshold^2){
      #W2 coupling
      flag <- "W2"
      
      #each iteration is given by n steps
      choices <- sample(1:n, n, replace = T)
      U <- runif(n)
      for(l in 1:n){
        i <- choices[l]
        
        #z1
        #compute conditional mean
        mean_i <- -cond_var[i]*sum(Q[i,-i]*z1[-i])
        #sample z_i
        z1[i] <- qtruncnorm(U[l], a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
        
        #z2
        #compute conditional mean
        mean_i <- -cond_var[i]*sum(Q[i,-i]*z2[-i])
        #sample z_i
        z2[i] <- qtruncnorm(U[l], a = lower_lim[i], b = upper_lim[i], mean = mean_i, sd = sqrt(cond_var[i]))
      }
    }
    else if(dist2 < threshold^2){
      #maximal coupling
      flag <- "Maximal"
      
      #each iteration is given by n steps
      choices <- sample(1:n, n, replace = T)
      for(l in 1:n){
        i <- choices[l]
        
        #z1
        #compute conditional mean
        mean_i1 <- -cond_var[i]*sum(Q[i,-i]*z1[-i])
        
        
        #z2
        #compute conditional mean
        mean_i2 <- -cond_var[i]*sum(Q[i,-i]*z2[-i])
        
        #sample z1
        z1[i] <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i]))
        log1 <- log(dtruncnorm(z1[i], a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i])))
        log2 <- log(dtruncnorm(z1[i], a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i])))
        if(log(runif(1)) <= log2-log1){
          z2[i] <- z1[i]
        }
        else{
          reject <- T
          while(reject){
            #propose z2
            prop <- rtruncnorm(1, a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i]))
            log1 <- log(dtruncnorm(prop, a = lower_lim[i], b = upper_lim[i], mean = mean_i1, sd = sqrt(cond_var[i])))
            log2 <- log(dtruncnorm(prop, a = lower_lim[i], b = upper_lim[i], mean = mean_i2, sd = sqrt(cond_var[i])))
            if(log(runif(1)) <= log1 - log2){
              z2[i] <- prop
              reject <- F
            }
          }
        }
        
      }
    }
    #compute distance and check
    dist2 <- sum((z1-z2)^2)
    if(verbose){
      print(paste("Dist. iter.", n_iter-L, ":", round(sqrt(dist2), 3), "with", flag))
    }
    #check if they are equal
    different <- (dist2 > 1e-15)
  }
  
  if(return_sample){
    return(c(n_iter, z1))
  }
  return(n_iter)
}