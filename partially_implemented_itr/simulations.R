source("algo1.R")

expit <- function(x) 1/(1+exp(-x))

library(pracma)
library(mvtnorm)
library(boot)

set.seed(1934)

# Create a large super population
largeN <- 2*10^6

# True parameters for Ya0|X=x
alpha <- c(0, -3, -.5, 5, -1.5, -2, 0) / 10

# True parameters for Ya1|X=x
beta <- c(0, -2, .5, 3, -1, -1, 0) / 10

# True parameters for S|X=x
gamma <- c(0, 0, 0, 0, 0, 0, 1)

# True parameters for r(X)
delta <- c(.05, -.5, .5, -.5, .5, .0, 0)

# True parameters for Psi(X)  
zeta <- delta 

#
nvar <- 6
O <- randortho(nvar, type = "orthonormal")
D_mat <- matrix(0, nrow = nvar, ncol = nvar)
diag(D_mat) <- seq(1, by=.2, length=nvar)
Sigma <- O %*% D_mat %*% t(O)


dgp <- function(n, alpha, beta, gamma, delta, zeta, Sigma) {
  
  # Function that takes sample size (n) and true value of the parameters and outputs
  # a simulated data in dataframe format
  
  # Generate Xs
  nvar <- 6
  X_temp <- rmvnorm(n, mean = rep(0, nrow(Sigma)), sigma = Sigma)
  X_temp[,1:2] <- as.numeric(X_temp[,1:2] < 0)
  X_temp[,3:5] <- exp(X_temp[,3:5])
  X <- cbind(1, X_temp)
  
  # Generate Ss
  rho <- expit(X %*% gamma)
  S <- rbinom(n, 1, rho)
  
  # Generate r(X)
  r <-  as.numeric(X %*% delta < 0)
  
  # Generate potential outcomes A_s0 and A_s1
  psi <-  expit(X %*% zeta)
  A_s0 <- rbinom(n, 1, psi)
  
  A_s1 <- r
  
  # Generate potential outcomes Y_a0 and Y_a0
  Pr_Ya0 <- expit(X %*% alpha)
  Pr_Ya1 <- expit(X %*% beta)
  
  Y_a0 <- rbinom(n, 1, Pr_Ya0)
  Y_a1 <- rbinom(n, 1, Pr_Ya1)
  
  # Generate potential outcomes Y_s0 and Y_s0
  Y_s0 <- A_s0*Y_a1 + (1-A_s0)*Y_a0
  Y_s1 <- A_s1*Y_a1 + (1-A_s1)*Y_a0
  
  # Generate observed outcome Y (under partial implementation of an ITR)
  Y <- S*Y_s1 + (1-S)*Y_s0
  
  # Generate observed outcome A (under partial implementation of an ITR)
  A <- S*A_s1 + (1-S)*A_s0
  
  # Create C
  C <- as.numeric(r==A)
  
  return(data.frame(X=X, rho=rho, S=S, Pr_Ya1=Pr_Ya1, Pr_Ya0=Pr_Ya0, Y_a1=Y_a1, Y_a0=Y_a0, Y_s1=Y_s1, Y_s0=Y_s0, Y=Y, r=r, A_s0=A_s0, A=A, psi=psi, C=C))
}

seed <- 1934
set.seed(seed)

### Simulation setup ###

n_cores <- parallel::detectCores(logical = FALSE)

n_sim <- 10
boot_resample <- 9
prec <- 6
sample_sizes <- c(200, 800, 2000)

start_time <- NULL
end_time <- NULL

a_var <- "A"
r_var <- "r"
y_var <- "Y"
mu_vars <- c("X.2", "X.3", "X.4", "X.5", "X.6")
psi_vars <- c("X.1", "X.2", "X.3", "X.4", "X.5") 
rho_vars <- c("X.7")

### Psi and r very different i.e., zeta = delta ### 

superpop_psi_different <- dgp(largeN, alpha, beta, gamma, delta, zeta=delta, Sigma)
TrueVals_psi_different <- apply(superpop_psi_different[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_psi_different
TrueEstimands_psi_different <- data.frame(True_MIG=TrueVals_psi_different["Y_s1"]-TrueVals_psi_different["Y"],
                                          True_ARE=TrueVals_psi_different["Y_s1"]-TrueVals_psi_different["Y_s0"],
                                          True_AIE=TrueVals_psi_different["Y"]-TrueVals_psi_different["Y_s0"])

TrueEstimands_psi_different

estimates_psi_different <- list()
lb_ci_psi_different <- list()
ub_ci_psi_different <- list()
#n_viol_different <- list()

ss_it <- 0
for (sample_size in sample_sizes){
  ss_it <- ss_it + 1
  
  estimates_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  lb_ci_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  ub_ci_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  
  colnames(estimates_psi_different[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(lb_ci_psi_different[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(ub_ci_psi_different[[ss_it]]) <- c("MIG", "ARE", "AIE")
  
  x <- 1:nrow(superpop_psi_different);
  y <- seq(from = 1, to = nrow(superpop_psi_different), by = sample_size);
  indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)
  
  for (i in 1:n_sim) {
    
    cat("last iteration took", as.numeric(difftime(end_time, start_time, unit="min")), "minutes\n")
    start_time <- Sys.time()
    cat("Currently running iteration", i, "of", n_sim, 
        "\nSample Size: ", sample_size,
        "\nr and psi different (situation 1)\n")
    
    res <- boot(superpop_psi_different[indices[, i], ], parallel="snow", ncpus=n_cores,
                estimators_boot,
                R=boot_resample,
                mu_vars=mu_vars, psi_vars=psi_vars, rho_vars=rho_vars, a_var=a_var, r_var=r_var, y_var=y_var, maxit=10, epsi = 10^-3, verbose=TRUE)
    
    estimates_psi_different[[ss_it]][i,] <- round(res$t0[,1:3], prec)
    ci <- round(emp_boot_ci(t0 = res$t0[,1:3], t = res$t[,1:3]), prec)
    lb_ci_psi_different[[ss_it]][i,] <- ci[1,]
    ub_ci_psi_different[[ss_it]][i,] <- ci[2,]
      
    end_time <- Sys.time()
    cat('\014')
  }
}

#save.image(file="end_first_scenario.RData")

gc() # Free memory

### Psi is random i.e., zeta = rep(0,7) ### 

superpop_psi_random <- dgp(largeN, alpha, beta, gamma, delta, zeta=rep(0,7), Sigma)
TrueVals_psi_random <- apply(superpop_psi_random[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_psi_random
TrueEstimands_psi_random <- data.frame(True_MIG=TrueVals_psi_random["Y_s1"]-TrueVals_psi_random["Y"],
                                       True_ARE=TrueVals_psi_random["Y_s1"]-TrueVals_psi_random["Y_s0"],
                                       True_AIE=TrueVals_psi_random["Y"]-TrueVals_psi_random["Y_s0"])
TrueEstimands_psi_random

estimates_psi_random <- list()
lb_ci_psi_random <- list()
ub_ci_psi_random <- list()

ss_it <- 0
for (sample_size in sample_sizes){
  ss_it <- ss_it + 1
  
  estimates_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  lb_ci_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  ub_ci_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  
  colnames(estimates_psi_random[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(lb_ci_psi_random[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(ub_ci_psi_random[[ss_it]]) <- c("MIG", "ARE", "AIE")
  
  x <- 1:nrow(superpop_psi_random);
  y <- seq(from = 1, to = nrow(superpop_psi_random), by = sample_size);
  indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)
  
  for (i in 1:n_sim) {
    
    cat("last iteration took", as.numeric(difftime(end_time, start_time, unit="min")), "minutes\n")
    start_time <- Sys.time()
    cat("Currently running iteration", i, "of", n_sim, 
        "\nSample Size: ", sample_size,
        "\npsi random (situation 2)\n")
    
    res <- boot(superpop_psi_random[indices[, i], ], parallel="snow", ncpus=n_cores,
                estimators_boot,
                R=boot_resample,
                mu_vars=mu_vars, psi_vars=psi_vars, rho_vars=rho_vars, a_var=a_var, r_var=r_var, y_var=y_var, maxit=10, epsi = 10^-3, verbose=TRUE)
    
    estimates_psi_random[[ss_it]][i,] <- round(res$t0[,1:3], prec)
    ci <- round(emp_boot_ci(t0 = res$t0[,1:3], t = res$t[,1:3]), prec)
    lb_ci_psi_random[[ss_it]][i,] <- ci[1,]
    ub_ci_psi_random[[ss_it]][i,] <- ci[2,]
    
    end_time <- Sys.time()
    cat('\014')
  }
}

#save.image(file="end_second_scenario.RData")

### Psi is similar to r i.e., zeta = - delta ### 

superpop_psi_similar <- dgp(largeN, alpha, beta, gamma, delta, zeta=- delta, Sigma)
TrueVals_psi_similar <- apply(superpop_psi_similar[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_psi_similar
TrueEstimands_psi_similar <- data.frame(True_MIG=TrueVals_psi_similar["Y_s1"]-TrueVals_psi_similar["Y"],
                                        True_ARE=TrueVals_psi_similar["Y_s1"]-TrueVals_psi_similar["Y_s0"],
                                        True_AIE=TrueVals_psi_similar["Y"]-TrueVals_psi_similar["Y_s0"])
TrueEstimands_psi_similar

estimates_psi_similar <- list()
lb_ci_psi_similar <- list()
ub_ci_psi_similar <- list()

ss_it <- 0
for (sample_size in sample_sizes){
  ss_it <- ss_it + 1
  
  estimates_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  lb_ci_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  ub_ci_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=3)
  
  colnames(estimates_psi_similar[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(lb_ci_psi_similar[[ss_it]]) <- c("MIG", "ARE", "AIE")
  colnames(ub_ci_psi_similar[[ss_it]]) <- c("MIG", "ARE", "AIE")
  
  x <- 1:nrow(superpop_psi_similar);
  y <- seq(from = 1, to = nrow(superpop_psi_similar), by = sample_size);
  indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)
  
  for (i in 1:n_sim) {
    
    cat("last iteration took", as.numeric(difftime(end_time, start_time, unit="min")), "minutes\n")
    start_time <- Sys.time()
    cat("Currently running iteration", i, "of", n_sim, 
        "\nSample Size: ", sample_size,
        "\nr and psi similar (situation 3)\n")
    
    res <- boot(superpop_psi_similar[indices[, i], ], parallel="snow", ncpus=n_cores,
                estimators_boot,
                R=boot_resample,
                mu_vars=mu_vars, psi_vars=psi_vars, rho_vars=rho_vars, a_var=a_var, r_var=r_var, y_var=y_var, maxit=10, epsi = 10^-2, verbose=TRUE)
    
    estimates_psi_similar[[ss_it]][i,] <- round(res$t0[,1:3], prec)
    ci <- round(emp_boot_ci(t0 = res$t0[,1:3], t = res$t[,1:3]), prec)
    lb_ci_psi_similar[[ss_it]][i,] <- ci[1,]
    ub_ci_psi_similar[[ss_it]][i,] <- ci[2,]
    
    end_time <- Sys.time()
    cat('\014')
  }
}

save.image(file="end_third_scenario.RData")

