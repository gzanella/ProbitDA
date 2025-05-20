library(truncnorm)
source("Functions_coupling.R")


r <- 3 #n/p
P <- c(50, 100, 200) #number of covariates p
N <- ceiling(r*P) #number of observations n
C <- c(1) #prior parameter c

L_val <- c(200, 200, 200, 200) #lags
max_iter <- 1000 #maximum number of iterations for coupling
threshold <- 1/10 #threshold for coupling
print_message <- T
len_to_print <- 100


replicas <- 500 #number of replications
meeting_DA <- list()
meeting_GS <- list()
bound_DA <- list()
bound_GS <- list()
for(j in 1:length(C)){
  meeting_DA[[j]] <- matrix(0, replicas, length(P))
  meeting_GS[[j]] <- matrix(0, replicas, length(P))
  bound_DA[[j]] <- matrix(0, replicas, length(P))
  bound_GS[[j]] <- matrix(0, replicas, length(P))
}

for(j in 1:length(C)){
  c <- C[j]
  for(k in 1:length(P)){
    n <- N[k]
    p <- P[k]
    L <- L_val[k]
    
    #generate X
    X <- matrix(rnorm(n*(p-1)), n, p-1)
    X <- X/sqrt(p)
    X <- cbind(1, X)
    
    #prior variance
    Sigma0 <- c*diag(p)
    
    #generate y: COMMENT THE ONE YOU DO NOT USE
    #all equal y
    y <- rep(1,n)
    #average y
    beta <- t(chol(Sigma0))%*%rnorm(p, 0, 1)
    y <- ifelse(X%*%beta + rnorm(n, 0, 1) > 0, 1, 0)
    
    for(r in 1:replicas){
      
      #compute theoretical bound: DA
      H <- X%*%Sigma0%*%t(X)
      H <- 0.5*H+0.5*t(H)
      gap <- 1/(1+max(eigen(H)$values))
      bound_DA[[j]][r, k] <- 1/gap
      
      
      #compute theoretical bound: GS
      Q <- chol2inv(chol(diag(n)+H))
      D <- diag(1/sqrt(diag(Q)))
      R <- D%*%Q%*%D
      R <- 0.5*R+0.5*t(R)
      gap <- min(eigen(R)$values)
      bound_GS[[j]][r, k] <- 1/gap
      
      if(n >= p){
        #run DA
        res_DA <- coupling_DA_n_bigger_p(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][r, k] <- res_DA
        
        
        #run GS
        res_GS <- coupling_GS_n_bigger_p(X, y, Sigma0, threshold/100, L, max_iter, verbose = F)
        meeting_GS[[j]][r, k] <- res_GS
      }
      else{
        #run DA
        res_DA <- coupling_DA_p_bigger_n(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][r, k] <- res_DA
        
        
        #run GS
        res_GS <- coupling_GS_p_bigger_n(X, y, Sigma0, threshold/100, L, max_iter, verbose = F)
        meeting_GS[[j]][r, k] <- res_GS
      }
      
      
      
      
      #print
      if(print_message & r %% len_to_print == 0){
        print(paste("Replica", r, "with n =", n, "p =", p, "and", "c =", c))
        to_print <- c(bound_DA[[j]][r, k], meeting_DA[[j]][r, k]-L,  bound_GS[[j]][r, k], meeting_GS[[j]][r, k]-L)
        names(to_print) <- c("Bound DA", "DA", "Bound GS",  "GS")
        print(round(to_print, 1))
        cat("\n")
      }
    }
  }
}
#saveRDS(meeting_DA, "meeting_DA_noint_imb.rds")
#saveRDS(meeting_GS, "meeting_GS_noint_imb.rds")
#saveRDS(meeting_DA, "meeting_DA_noint_average.rds")
#saveRDS(meeting_GS, "meeting_GS_noint_average.rds")

#tables
epsilon <- 0.1

#c = 1
ind_c <- 1
perf_DA1 <- meeting_DA[[ind_c]][,1]
perf_GS1 <- meeting_GS[[ind_c]][,1]
perf_DA2 <- meeting_DA[[ind_c]][,2]
perf_GS2 <- meeting_GS[[ind_c]][,2]
perf_DA3 <- meeting_DA[[ind_c]][,3]
perf_GS3 <- meeting_GS[[ind_c]][,3]

Times <- seq(0, 1000, 1)
TV_boundsDA1 <- rep(0, length(Times))
TV_boundsGS1 <- rep(0, length(Times))
TV_boundsDA2 <- rep(0, length(Times))
TV_boundsGS2 <- rep(0, length(Times))
TV_boundsDA3 <- rep(0, length(Times))
TV_boundsGS3 <- rep(0, length(Times))
for(i in 1:length(Times)){
  t <- Times[i]
  
  #1
  L <- L_val[1]
  maxs <- pmax(0, ceiling((perf_DA1-L-t)/L))
  TV_boundsDA1[i] <- mean(maxs)
  
  maxs <- pmax(0, ceiling((perf_GS1-L-t)/L))
  TV_boundsGS1[i] <- mean(maxs)
  
  #2
  L <- L_val[2]
  maxs <- pmax(0, ceiling((perf_DA2-L-t)/L))
  TV_boundsDA2[i] <- mean(maxs)
  
  maxs <- pmax(0, ceiling((perf_GS2-L-t)/L))
  TV_boundsGS2[i] <- mean(maxs)
  
  #3
  L <- L_val[3]
  maxs <- pmax(0, ceiling((perf_DA3-L-t)/L))
  TV_boundsDA3[i] <- mean(maxs)
  
  maxs <- pmax(0, ceiling((perf_GS3-L-t)/L))
  TV_boundsGS3[i] <- mean(maxs)
}

#DA
which(TV_boundsDA1 <= epsilon)[1] - 1
which(TV_boundsDA2 <= epsilon)[1] - 1
which(TV_boundsDA3 <= epsilon)[1] - 1

#GS
which(TV_boundsGS1 <= epsilon)[1] - 1
which(TV_boundsGS2 <= epsilon)[1] - 1
which(TV_boundsGS3 <= epsilon)[1] - 1
