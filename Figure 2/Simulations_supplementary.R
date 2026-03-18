library(truncnorm)
source("Functions_coupling.R")


r <- c(0.3, 0.6, 0.9, 1, 1.5, 2, 3) #n/p
p <- c(200) #number of covariates p
N <- ceiling(r*p) #number of observations n
C <- c(10) #prior parameter c

L_val <- rep(200, length(r))
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
  meeting_DA[[j]] <- matrix(0, replicas, length(r))
  meeting_GS[[j]] <- matrix(0, replicas, length(r))
  bound_DA[[j]] <- matrix(0, replicas, length(r))
  bound_GS[[j]] <- matrix(0, replicas, length(r))
}

for(j in 1:length(C)){
  c <- C[j]
  for(k in 1:length(r)){
    n <- N[k]
    L <- L_val[k]
    
    #generate X
    X <- matrix(rnorm(n*p), n, p)
    X <- X/sqrt(p)
    
    #prior variance
    Sigma0 <- c*diag(p)
    
    #generate y: COMMENT THE ONE YOU DO NOT USE
    #all equal y
    #y <- rep(1,n)
    #average y
    beta <- t(chol(Sigma0))%*%rnorm(p, 0, 1)
    y <- ifelse(X%*%beta + rnorm(n, 0, 1) > 0, 1, 0)
    
    for(t in 1:replicas){
      
      #compute bounds: DA
      H <- X%*%Sigma0%*%t(X)
      H <- 0.5*H+0.5*t(H)
      gap <- 1/(1+max(eigen(H)$values))
      bound_DA[[j]][t, k] <- 1/gap
      
      
      #compute bounds: GS
      Q <- chol2inv(chol(diag(n)+H))
      D <- diag(1/sqrt(diag(Q)))
      R <- D%*%Q%*%D
      R <- 0.5*R+0.5*t(R)
      gap <- min(eigen(R)$values)
      bound_GS[[j]][t, k] <- 1/gap
      #bound_GS2 <- (1+max(eigen(H)$values))/(1+min(eigen(H)$values))
      
      if(n >= p){
        #run DA
        res_DA <- coupling_DA_n_bigger_p(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][t, k] <- res_DA
        
        
        #run GS
        res_GS <- coupling_GS_n_bigger_p(X, y, Sigma0, threshold/100, L, max_iter, verbose = F)
        meeting_GS[[j]][t, k] <- res_GS
      }
      else{
        #run DA
        res_DA <- coupling_DA_p_bigger_n(X, y, Sigma0, threshold, L, max_iter, verbose = F)
        meeting_DA[[j]][t, k] <- res_DA
        
        
        #run GS
        res_GS <- coupling_GS_p_bigger_n(X, y, Sigma0, threshold/100, L, max_iter, verbose = F)
        meeting_GS[[j]][t, k] <- res_GS
      }
      
      
      
      
      #print
      if(print_message & t %% len_to_print == 0){
        print(paste("Replica", t, "with n =", n, "p =", p, "and", "c =", c))
        to_print <- c(bound_DA[[j]][t, k], meeting_DA[[j]][t, k]-L,  bound_GS[[j]][t, k], meeting_GS[[j]][t, k]-L)
        names(to_print) <- c("Bound DA", "DA", "Bound GS",  "GS")
        print(round(to_print, 1))
        cat("\n")
      }
    }
  }
}
#saveRDS(meeting_DA, "meeting_DA_average_supp.rds")
#saveRDS(meeting_GS, "meeting_GS_average_supp.rds")

#plot
epsilon <- 0.1

#c = 10
ind_c <- 1
perf_DA <- meeting_DA[[ind_c]]
perf_GS <- meeting_GS[[ind_c]]

Times <- seq(0, 150, 1)
TV_boundsDA <- matrix(0, length(Times), length(r))
TV_boundsGS <- matrix(0, length(Times), length(r))
for(i in 1:length(Times)){
  t <- Times[i]
  
  for(j in 1:length(r)){
    L <- L_val[j]
    maxs <- pmax(0, ceiling((perf_DA[,j]-L-t)/L))
    TV_boundsDA[i, j] <- mean(maxs)
    
    maxs <- pmax(0, ceiling((perf_GS[,j]-L-t)/L))
    TV_boundsGS[i, j] <- mean(maxs)
  }
}

#mixing times
mixing_DA <- apply(TV_boundsDA, 2, function(x) which(x <= epsilon)[1] - 1)
mixing_GS <- apply(TV_boundsGS, 2, function(x) which(x <= epsilon)[1] - 1)

#plot
Max <- max(c(mixing_DA, mixing_GS))
Min <- min(c(mixing_DA, mixing_GS))
plot(r, mixing_DA, type = "b", lwd = 2, pch = 19, col = gray(0), ylim = c(Min, Max), xlab = "n/p", ylab = "Upper bound mixing times", cex = 1.5, cex.lab = 1.4, cex.axis = 2)
points(r, mixing_GS, type = "b", lwd = 2, pch = 18, ylim = c(Min, Max), col = gray(0.6), cex = 1.5)
legend(1.7, 80, legend = c("DA", "Collapsed"), 
       col = c(gray(0), gray(0.6)), lwd = 2, pch = c(19, 18), bty = "n", cex = 1.5, y.intersp = 2)

