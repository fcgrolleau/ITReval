cb_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it_boot <<- it_boot + 1 
  
  ###
  # Fit propensity score model
  ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
                data=z, family = "binomial")
  
  # Get propensity score predictions
  z$ps_hat <-  ps_mod$fitted.values
  
  
  # Fit prognostic model
  pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
                data=z, family = "binomial")
  
  # Get predictions for the ITE
  newdata <- z[, c("a1", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
  newdata[,"a1"] <- 1
  z$E_hat_Ya1_given_x <- predict(pr_mod, newdata, type="response")
  
  newdata[,"a1"] <- 0
  z$E_hat_Ya0_given_x <- predict(pr_mod, newdata, type="response")
  
  z$ITE_hat <- with(z, E_hat_Ya1_given_x - E_hat_Ya0_given_x)
  
  ###
  
  # cognitive bias scenario

  z$rho <- with(z, (1- abs(r - ps_hat) )^(.5*log((alpha_cb+1)/(1-alpha_cb))) )
  Lambda <-  mean( with(z, rho*(r - ps_hat ) * ITE_hat))
  
  if ( (100*it_boot/resamples) %% 10 == 0) {
    print(paste0('Iteration alpha ', it, ' Iteration bootstrap ', it_boot, ': ', 100*it_boot/resamples, '%')) }
  return (Lambda)
}

cl_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it_boot <<- it_boot + 1 
  
  ###
  # Fit propensity score model
  ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
                data=z, family = "binomial")
  
  # Get propensity score predictions
  z$ps_hat <-  ps_mod$fitted.values
  
  
  # Fit prognostic model
  pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
                data=z, family = "binomial")
  
  # Get predictions for the ITE
  newdata <- z[, c("a1", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
  newdata[,"a1"] <- 1
  z$E_hat_Ya1_given_x <- predict(pr_mod, newdata, type="response")
  
  newdata[,"a1"] <- 0
  z$E_hat_Ya0_given_x <- predict(pr_mod, newdata, type="response")
  
  z$ITE_hat <- with(z, E_hat_Ya1_given_x - E_hat_Ya0_given_x)
  
  ###
  
  # confidence level scenario
  
  z$rho <- with(z, (iard - qnorm(1 - alpha_cl/2))*(iard + qnorm(1 - alpha_cl/2)*ise) > 0 )
  Lambda <-  mean( with(z, rho*(r - ps_hat ) * ITE_hat)) 
  
  if ( (100*it_boot/resamples) %% 10 == 0) {
    print(paste0('Iteration alpha ', it, ' Iteration bootstrap ', it_boot, ': ', 100*it_boot/resamples, '%')) }
  return (Lambda)
} 

emp_boot_ci <-  function(t0, t, alpha=.05){
  # Empirical Bootstrap CI
  # as described in John Rice, Mathematical Statistics and Data Analysis, 3rd edition, p. 285.        
  return( as.numeric( c(2*t0 - quantile(t, probs = 1-alpha/2), 2*t0 - quantile(t, probs = alpha/2)) ) )
}
