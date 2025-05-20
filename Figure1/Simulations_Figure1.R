library(truncnorm)
source("Functions_coupling.R")
P <- c(200, 400, 600) #number of covariates p
N <- c(20, 40, 60) #number of observations n
C <- c(1) #prior parameter c

L_val <- c(200, 200, 200) #lags
max_iter <- 1000 #maximum number of iterations for coupling
threshold <- 1/10 #threshold for coupling
print_message <- T
len_to_print <- 100


replicas <- 200 #number of replications
meeting_DA <- list()
meeting_DA_nonC <- list()
bound_DA <- list()
for(j in 1:length(C)){
  meeting_DA[[j]] <- matrix(0, replicas, length(P))
  meeting_DA_nonC[[j]] <- matrix(0, replicas, length(P))
  bound_DA[[j]] <- matrix(0, replicas, length(P))
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
    
    #generate the y: COMMENT THE ONE YOU DO NOT USE
    #all equal y
    y <- rep(1,n)
    #average y
    #beta <- t(chol(Sigma0))%*%rnorm(p, 0, 1)
    #y <- ifelse(X%*%beta + rnorm(n, 0, 1) > 0, 1, 0)
    for(r in 1:replicas){
      
      #compute theoretical bounds: DA
      H <- X%*%Sigma0%*%t(X)
      H <- 0.5*H+0.5*t(H)
      gap <- 1/(1+max(eigen(H)$values))
      bound_DA[[j]][r, k] <- 1/gap
      
      
      
      if(n >= p){
        #run DA
        res_DA <- coupling_DA_n_bigger_p(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][r, k] <- res_DA
        
        #run DA nonC
        res_DA_nonC <- coupling_DA_n_bigger_p_noncIntercept(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA_nonC[[j]][r, k] <- res_DA_nonC
        
      }
      else{
        #run DA
        res_DA <- coupling_DA_p_bigger_n(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][r, k] <- res_DA
        
        #run DA nonC
        res_DA_nonC <- coupling_DA_p_bigger_n_noncIntercept(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA_nonC[[j]][r, k] <- res_DA_nonC
        
      }
      
      
      
      
      #print
      if(print_message & r %% len_to_print == 0){
        print(paste("Replica", r, "with n =", n, "p =", p, "and", "c =", c))
        to_print <- c(bound_DA[[j]][r, k], meeting_DA[[j]][r, k]-L,meeting_DA_nonC[[j]][r, k]-L)
        names(to_print) <- c("Bound DA", "DA", "DA_nonC")
        print(round(to_print, 1))
        cat("\n")
      }
    }
  }
}
#saveRDS(meeting_DA, "meeting_DA_withint_average.rds")
#saveRDS(meeting_DA_nonC, "meeting_DA_nonC_withint_average.rds")

#plots


#c <- 1
ind_c <- 1
perf_DA1 <- meeting_DA[[ind_c]][,1]
perf_DA_nonC1 <- meeting_DA_nonC[[ind_c]][,1]
perf_DA2 <- meeting_DA[[ind_c]][,2]
perf_DA_nonC2 <- meeting_DA_nonC[[ind_c]][,2]
perf_DA3 <- meeting_DA[[ind_c]][,3]
perf_DA_nonC3 <- meeting_DA_nonC[[ind_c]][,3]

Times <- seq(0, 100, 1)
TV_boundsDA1 <- rep(0, length(Times))
TV_boundsDA_nonC1 <- rep(0, length(Times))
TV_boundsDA2 <- rep(0, length(Times))
TV_boundsDA_nonC2 <- rep(0, length(Times))
TV_boundsDA3 <- rep(0, length(Times))
TV_boundsDA_nonC3 <- rep(0, length(Times))
for(i in 1:length(Times)){
  t <- Times[i]
  
  #1
  L <- L_val[1]
  maxs <- pmax(0, ceiling((perf_DA1-L-t)/L))
  TV_boundsDA1[i] <- mean(maxs)
  
  maxs <- pmax(0, ceiling((perf_DA_nonC1-L-t)/L))
  TV_boundsDA_nonC1[i] <- mean(maxs)
  
  #2
  L <- L_val[2]
  maxs <- pmax(0, ceiling((perf_DA2-L-t)/L))
  TV_boundsDA2[i] <- mean(maxs)
  
  maxs <- pmax(0, ceiling((perf_DA_nonC2-L-t)/L))
  TV_boundsDA_nonC2[i] <- mean(maxs)
  
  
  #3
  L <- L_val[3]
  maxs <- pmax(0, ceiling((perf_DA3-L-t)/L))
  TV_boundsDA3[i] <- mean(maxs)
  
  
  maxs <- pmax(0, ceiling((perf_DA_nonC3-L-t)/L))
  TV_boundsDA_nonC3[i] <- mean(maxs)
  
}


#only DA
Max <- max(c(TV_boundsDA1, TV_boundsDA2, TV_boundsDA3))
Min <- min(c(TV_boundsDA1, TV_boundsDA2, TV_boundsDA3))
plot(Times, TV_boundsDA1, type = "b", lwd = 2, pch = 19, col = gray(0), ylim = c(Min, Max), xlab = "Iteration", ylab = "TV distance from stat.", cex = 1.5, cex.lab = 1.4, cex.axis = 2)
points(Times, TV_boundsDA2, type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.6), cex = 1.5)
points(Times, TV_boundsDA3, type = "b", lwd = 2, pch = 17, ylim = c(Min, Max), col = gray(0.3), cex = 1.5)
legend(40, 0.8, legend = c("n = 20, p = 200", "n = 40, p = 400", "n = 60, p = 600"), 
       col = c(gray(0), gray(0.6), gray(0.3)), lwd = 2, pch = c(19, 18, 17), bty = "n", cex = 1.5, y.intersp = 2)



#only DA nonC
Max <- max(c(TV_boundsDA_nonC1, TV_boundsDA_nonC2, TV_boundsDA_nonC3))
Min <- min(c(TV_boundsDA_nonC1, TV_boundsDA_nonC2, TV_boundsDA_nonC3))
plot(Times, TV_boundsDA_nonC1, type = "b", lwd = 2, pch = 19, col = gray(0), ylim = c(Min, Max), xlab = "Iteration", ylab = "TV distance from stat.", cex = 1.5, cex.lab = 1.4, cex.axis = 2)
points(Times, TV_boundsDA_nonC2, type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.6), cex = 1.5)
points(Times, TV_boundsDA_nonC3, type = "b", lwd = 2, pch = 17, ylim = c(Min, Max), col = gray(0.3), cex = 1.5)
legend(40, 0.8, legend = c("n = 20, p = 200", "n = 40, p = 400", "n = 60, p = 600"), 
       col = c(gray(0), gray(0.6), gray(0.3)), lwd = 2, pch = c(19, 18, 17), bty = "n", cex = 1.5, y.intersp = 2)
