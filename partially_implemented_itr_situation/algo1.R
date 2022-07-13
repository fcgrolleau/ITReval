estimators_boot <- function(d, i=1:nrow(d),
                            mu_vars, psi_vars, rho_vars, a_var, r_var, y_var, maxit, epsi, verbose) {
  
  # Bootstrap function that explicitly contains the ests_fun function 
  # so as to make Parallel Processing compatible with both Windows and Max / Linux OS.
  
  z<-d[i,]
  
  expit <- function(x) 1/(1+exp(-x))
  
  ests_fun <- function(dat, mu_vars, psi_vars, rho_vars, a_var, r_var, y_var, maxit, epsi, verbose=TRUE)
  {
  
  # Function that takes as input data originating from a partially implemented ITR situation
  # and using an EM algorithm outputs MIG, ARE, and AIE estimates.
  
  ### INPUTS ###  
  # dat: (dataframe)
  # a_var: (character) treatment variable name in dat
  # r_var: (character) ITR variable name in dat
  # y_var: (character) outcome variable name in dat
  # mu_vars: (concatenation of characters) prognostic variable names in dat
  # psi_vars: (concatenation of characters) propensity score in the absence of ITR implementation variable names in dat
  # rho_vars: (concatenation of characters) stochastic implementation function variable name in dat
  # maxit: (num) maximum number of EM iterations
  # epsi: (num) convergence threshold for the parameters ζ
  # verbose: (logi) print each iteration
    
  ### OUTPUTS ###
  # (num [1, 1:3]) estimates of MIG, ARE and AIE
    
    A <- dat[,a_var]
    r <- dat[,r_var]
    Y <- dat[, y_var]
    
    # Initialize the prior probabilities associated with the nodes of the tree 
    g_0 <- .5; g_1 <- .5; 
    
    # Initialize the parameters ζ of the expert πs=0(·) at random
    psi_param <- rnorm(length(psi_vars), 0, 10000)
    
    # Compute individual contributions to r’s likelihood 
    P_1 <- r^A * (1-r)^(1-A)
    
    # Compute individual predictions from the initiated expert network πs=0(·)
    p_0 <- expit(as.matrix(dat[psi_vars]) %*% psi_param)
    
    # Iterate until convergence on the parameters ζ
    
    for (it in 1:maxit)
    {
      
      if(it>1){
        params_prev <- params_now
      }
      
      # Compute individual contributions to psi’s likelihood 
      P_0 <- p_0^A * (1-p_0)^(1-A)
      
      # Compute the posterior probabilities associated with the nodes of the tree
      h_0 <- (g_0 * P_0) / (g_0 * P_0 + g_1 * P_1)
      h_1 <- (g_1 * P_1) / (g_0 * P_0 + g_1 * P_1)
      
      h_0[P_0==0] <- 1 # for numerical stability as P_0==0 not possible from expit
      h_1[P_0==0] <- 0 # idem as above
      
      # For the gating network ρ(·) estimate parameters γ by solving the IRLS problem
      f_rho_1 <- formula(paste0("h_1", " ~ -1 +", paste(rho_vars, collapse=" + ")))
      mod_rho_1 <- glm(f_rho_1, data = dat, family = "quasibinomial", method = "glm.fit") # "glm.fit" for IRLS solver
      
      # For the expert network πs=0(·) estimate parameters ζ by solving the IRLS problem
      f_psi <- formula(paste0("A", " ~ -1 +", paste(psi_vars, collapse=" + ")))
      mod_psi <- glm(f_psi, weights = h_0, data = dat, family = "quasibinomial", method = "glm.fit") # "glm.fit" for IRLS solver
      
      # Update the prior probabilities associated with the nodes of the tree 
      g_1 <- predict(mod_rho_1, type="response")
      g_0 <- 1 - g_1
      
      # Update the predictions from the expert network πs=0(·)
      p_0 <- predict(mod_psi, newdata = dat, type="response")
      
      # Check for convergence on parameters ζ
      params_now <- c(coef(mod_rho_1), coef(mod_psi))
      if(it>1){
        params_dif <- params_now - params_prev
        if(sqrt(sum(params_dif %*% params_dif)) < epsi) {
          break }
      }
      last_it <- it
    }
    
    
    f_mu <- formula(paste0(y_var, " ~ -1 +", paste(mu_vars, collapse=" + ")))
    mod_mu0 <- glm(f_mu, data=dat[dat[,a_var]==0,], family = "binomial")
    mod_mu1 <- glm(f_mu, data=dat[dat[,a_var]==1,], family = "binomial")
    
    mu0_hat <- predict(mod_mu0, newdata = dat, type="response")
    mu1_hat <- predict(mod_mu1, newdata = dat, type="response")  
    
    Q_1 <- r * mu1_hat + (1-r) * mu0_hat
    Q_0 <- p_0 * mu1_hat + (1-p_0) * mu0_hat
    
    E_Ys_1 <- mean(Q_1)
    E_Y <- mean(Y)
    E_Ys_0 <- mean(Q_0)
    
    MIG <- E_Ys_1 - E_Y
    ARE <- E_Ys_1 - E_Ys_0
    AIE <- E_Y  - E_Ys_0 
    
    if(verbose){
      if(last_it<maxit){
        cat("No. iterations =", last_it, "\n") }
      else{
        cat("Max iterations reached:", "maxit =", maxit, "\n")
      }
    }
    
    return(cbind(MIG=MIG, ARE=ARE, AIE=AIE) )
  }
  
  return(ests_fun(z, mu_vars, psi_vars, rho_vars, a_var, r_var, y_var, maxit, epsi, verbose))
}

emp_boot_ci <-  function(t0, t, alpha=.05){
  # Empirical Bootstrap CI
  # as described in John Rice, Mathematical Statistics and Data Analysis, 3rd edition, p. 285.
  temp <- rbind(
    2*t0 - apply(t , 2, function(x) quantile(x, probs = 1-alpha/2) ),
    2*t0 - apply(t , 2, function(x) quantile(x, probs = alpha/2) )
  )
  
  return( temp ) 
}


## Example with Parallel Processing on Windows
## On a Mac or Linux replace parallel="snow" by parallel="multicore"

# library(boot)
# res <- boot(my_data,
#            estimators_boot,
#            a_var="A", r_var="r", y_var="Y",
#            mu_vars=c("X.2", "X.3", "X.4", "X.5", "X.6"),
#            psi_vars=c("X.1", "X.2", "X.3", "X.4", "X.5"),
#            rho_vars=c("X.7"),
#            maxit=10, epsi=10^-3, verbose=TRUE,
#            R=999, parallel="snow", ncpus=parallel::detectCores(logical = FALSE))

# estimates <- res$t0[,1:3]
# ci <- emp_boot_ci(t0 = res$t0[,1:3], t = res$t[,1:3])


