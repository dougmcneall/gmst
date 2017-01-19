# gmstARIMAX.R
# predict Global Mean Surface Temperature using a Bayesian
# ARIMAX model (ARIMA with external regressors)
# Doug McNeall & Andrew Parnell 

library(R2jags)
arimax_code = '
model
{
  # Set up residuals
  for(t in 1:q) {
  eps[t] <- y[t] - alpha
  }
  # Likelihood
  for (t in (q+1):T) {
  y[t] ~ dnorm(alpha + ar_mean[t] + ma_mean[t] + reg_mean[t], tau)
  ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
  ar_mean[t] <- inprod(phi, y[(t-p):(t-1)])
  reg_mean[t] <- inprod(beta, x[t,])
  eps[t] <- y[t] - alpha - ar_mean[t] - ma_mean[t] - reg_mean[t]
  }
  
  # Priors
  alpha ~ dnorm(0.0,0.01)
  for (i in 1:q) {
  theta[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:p) {
  phi[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:k) {
  beta[i] ~ dnorm(0.0,0.01)
  }
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
}
'

gmstARIMAX = function(model_code, xreg, hadcrut, ny, q=1, p=1){
  
  # xreg need to be as long as hadcrut$Year+ny
  #stopifnot (nrow(xreg)  ==  (length(hadcrut$Year)+ny), 'time mismatch') 
  year_future = (max(hadcrut$Year)+1):(max(hadcrut$Year)+ny)
  model_parameters =  c("y","alpha","theta","phi","beta","sigma")
  k = ncol(xreg)+1
  
  data_future = with(hadcrut,
                     list(T = nrow(hadcrut) + ny,
                          y = c(Anomaly, rep(NA,ny)),
                          x = matrix(cbind(c(Year,year_future), xreg), ncol = k),
                          q = q,
                          p = p,
                          k = k))
  
  run_future = jags(data = data_future,
                    parameters.to.save = model_parameters,
                    model.file = textConnection(model_code),
                    n.chains = 4,
                    n.iter = 10000,
                    n.burnin = 2000,
                    n.thin = 8)
  
  # Get the future values
  y_all = run_future$BUGSoutput$sims.list$y
  y_all_mean = apply(y_all,2,'mean')
  # Also create the upper/lower 95% CI values
  y_all_low = apply(y_all,2,'quantile',0.025) 
  y_all_high = apply(y_all,2,'quantile',0.975) 
  year_all = c(hadcrut$Year,year_future)
  
  alpha = run_future$BUGSoutput$sims.list$alpha
  beta = run_future$BUGSoutput$sims.list$beta
  theta = run_future$BUGSoutput$sims.list$theta
  phi = run_future$BUGSoutput$sims.list$phi
  sigma = run_future$BUGSoutput$sims.list$sigma
  
  return(list(year=year_all,
              mean=y_all_mean, lower=y_all_low, y = y_all, upper=y_all_high,
              alpha=alpha, beta=beta, theta=theta, phi = phi, sigma = sigma,
              all_out = run_future)
         )
}
