# alternative_orders.R
# Try the ARIMAX with a couple of different orders, to make
# sure we're not missing a better timeseries model.

library(R2jags)

# load ARIMAX prediction function
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/gmstARIMAX.R')

# load data
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/load_gmst_data.R')

rank.unc = function(pred, obs){
  # pred is sample of prediction
  ranks = rank(c(obs,pred))
  out = ranks[1] / length(pred)
  out
}


# ---------------------------------------------------------------------
# One-step-ahead prediction for 1985 - 1999, using ARIMAX
# models of different order
# ---------------------------------------------------------------------
# Here is the standard prediction


osayears = 1985:1999
obs = hadcrut$Anomaly[which(hadcrut$Year%in%osayears)]

# setup containers for the predictions.
osamean101 = rep(NA,15)
osalower101 = rep(NA,15)
osaupper101 = rep(NA,15)
osarank101 = rep(NA,15)


for(i in 1:(length(osayears))){
  ny = 1
  trainyears = 1892:((osayears[1] - 2)+i)
  predyears = 1892:((osayears[1] - 2)+i+ny)
  # offset enso by a year
  enso34years = 1891:((osayears[1] - 3)+i+ny)
  
  # end of last year's ENSO
  nino_pred = nino34ts[which(ninoyear%in%enso34years)]
  
  # last year's natural forcing and amo
  nat_forc  = predictors2[which(predictors2$YEAR%in%predyears), c(4,5) ]
  amo       = predictors2[which(predictors2$YEAR%in%enso34years), 2 ]
  
  # The projected anthropogenic forcing
  anth_forc = predictors2[which(predictors2$YEAR%in%predyears), 6 ]
  
  # Temperature data for the training period
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]
  
  # Bind the exogenous variables together
  xreg=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))
  
  pred = gmstARIMAX(model_code=arimax_code, xreg=xreg, hadcrut=hadcrut_train, ny=1,p=1, q=1)
  
  osamean101[i] = tail(pred$mean,1)
  osaupper101[i] = tail(pred$upper,1)
  osalower101[i] = tail(pred$lower,1)
  osarank101[i] = rank.unc(pred$y[, ncol(pred$y)], obs[i])
}

# Basic prediction error statistics
arimax_mae_101 = mean(abs(obs - osamean101))
arimax_rmse_101 = sqrt(mean((obs - osamean101)^2))


# setup containers for the predictions.
osamean303 = rep(NA,15)
osalower303 = rep(NA,15)
osaupper303 = rep(NA,15)
osarank303 = rep(NA,15)

for(i in 1:(length(osayears))){
  ny = 1
  trainyears = 1892:((osayears[1] - 2)+i)
  predyears = 1892:((osayears[1] - 2)+i+ny)
  # offset enso by a year
  enso34years = 1891:((osayears[1] - 3)+i+ny)
  
  # end of last year's ENSO
  nino_pred = nino34ts[which(ninoyear%in%enso34years)]
  
  # last year's natural forcing and amo
  nat_forc  = predictors2[which(predictors2$YEAR%in%predyears), c(4,5) ]
  amo       = predictors2[which(predictors2$YEAR%in%enso34years), 2 ]
  
  # The projected anthropogenic forcing
  anth_forc = predictors2[which(predictors2$YEAR%in%predyears), 6 ]
  
  # Temperature data for the training period
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]
  
  # Bind the exogenous variables together
  xreg=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))
  
  pred = gmstARIMAX(model_code=arimax_code, xreg=xreg, hadcrut=hadcrut_train, ny=1,p=3, q=3)
  
  osamean303[i] = tail(pred$mean,1)
  osaupper303[i] = tail(pred$upper,1)
  osalower303[i] = tail(pred$lower,1)
  osarank303[i] = rank.unc(pred$y[, ncol(pred$y)], obs[i])
}

# Basic prediction error statistics
arimax_mae_303 = mean(abs(obs - osamean303))
arimax_rmse_303 = sqrt(mean((obs - osamean303)^2))





# Here is when the Moving Average term is zero 
arimax_100_code = '
model
{
  # Likelihood
  for (t in (p+1):T) {
  y[t] ~ dnorm(alpha + ar_mean[t] + reg_mean[t], tau)
  ar_mean[t] <- inprod(phi, y[(t-p):(t-1)])
  reg_mean[t] <- inprod(beta, x[t,])
  }
  
  # Priors
  alpha ~ dnorm(0.0,0.01)
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

  gmstARIMAX = function(model_code, xreg, hadcrut, ny, p=1){
    
    # xreg need to be as long as hadcrut$Year+ny
    #stopifnot (nrow(xreg)  ==  (length(hadcrut$Year)+ny), 'time mismatch') 
    year_future = (max(hadcrut$Year)+1):(max(hadcrut$Year)+ny)
    model_parameters =  c("y","alpha","phi","beta","sigma")
    k = ncol(xreg)+1
    
    data_future = with(hadcrut,
                       list(T = nrow(hadcrut) + ny,
                            y = c(Anomaly, rep(NA,ny)),
                            x = matrix(cbind(c(Year,year_future), xreg), ncol = k),
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
    phi = run_future$BUGSoutput$sims.list$phi
    sigma = run_future$BUGSoutput$sims.list$sigma
    
    return(list(year=year_all,
                mean=y_all_mean, lower=y_all_low, y = y_all, upper=y_all_high,
                alpha=alpha, beta=beta, phi = phi, sigma = sigma,
                all_out = run_future)
    )
  }

# setup containers for the predictions.
osamean100 = rep(NA,15)
osalower100 = rep(NA,15)
osaupper100 = rep(NA,15)
osarank100 = rep(NA,15)

for(i in 1:(length(osayears))){
  ny = 1
  trainyears = 1892:((osayears[1] - 2)+i)
  predyears = 1892:((osayears[1] - 2)+i+ny)
  # offset enso by a year
  enso34years = 1891:((osayears[1] - 3)+i+ny)
  
  # end of last year's ENSO
  nino_pred = nino34ts[which(ninoyear%in%enso34years)]
  
  # last year's natural forcing and amo
  nat_forc  = predictors2[which(predictors2$YEAR%in%predyears), c(4,5) ]
  amo       = predictors2[which(predictors2$YEAR%in%enso34years), 2 ]
  
  # The projected anthropogenic forcing
  anth_forc = predictors2[which(predictors2$YEAR%in%predyears), 6 ]
  
  # Temperature data for the training period
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]
  
  # Bind the exogenous variables together
  xreg=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))
  
  pred = gmstARIMAX(model_code=arimax_100_code, xreg=xreg, hadcrut=hadcrut_train, ny=1)
  
  osamean100[i] = tail(pred$mean,1)
  osaupper100[i] = tail(pred$upper,1)
  osalower100[i] = tail(pred$lower,1)
  osarank100[i] = rank.unc(pred$y[, ncol(pred$y)], obs[i])
}

# Basic prediction error statistics
arimax_mae_100 = mean(abs(obs - osamean100))
arimax_rmse_100 = sqrt(mean((obs - osamean100)^2))

arimax_001_code = '
model
{
  # Set up residuals
  for(t in 1:q) {
  eps[t] <- y[t] - alpha
  }
  # Likelihood
  for (t in (q+1):T) {
  y[t] ~ dnorm(alpha + ma_mean[t] + reg_mean[t], tau)
  ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
  reg_mean[t] <- inprod(beta, x[t,])
  eps[t] <- y[t] - alpha - ma_mean[t] - reg_mean[t]
  }
  
  # Priors
  alpha ~ dnorm(0.0,0.01)
  for (i in 1:q) {
  theta[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:k) {
  beta[i] ~ dnorm(0.0,0.01)
  }
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
}
'

gmstARIMAX = function(model_code, xreg, hadcrut, ny, q=1){
  
  # xreg need to be as long as hadcrut$Year+ny
  #stopifnot (nrow(xreg)  ==  (length(hadcrut$Year)+ny), 'time mismatch') 
  year_future = (max(hadcrut$Year)+1):(max(hadcrut$Year)+ny)
  model_parameters =  c("y","alpha","theta","beta","sigma")
  k = ncol(xreg)+1
  
  data_future = with(hadcrut,
                     list(T = nrow(hadcrut) + ny,
                          y = c(Anomaly, rep(NA,ny)),
                          x = matrix(cbind(c(Year,year_future), xreg), ncol = k),
                          q = q,
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
  sigma = run_future$BUGSoutput$sims.list$sigma
  
  return(list(year=year_all,
              mean=y_all_mean, lower=y_all_low, y = y_all, upper=y_all_high,
              alpha=alpha, beta=beta, theta=theta, sigma = sigma,
              all_out = run_future)
  )
}

# setup containers for the predictions.
osamean001 = rep(NA,15)
osalower001 = rep(NA,15)
osaupper001 = rep(NA,15)
osarank001 = rep(NA,15)

for(i in 1:(length(osayears))){
  ny = 1
  trainyears = 1892:((osayears[1] - 2)+i)
  predyears = 1892:((osayears[1] - 2)+i+ny)
  # offset enso by a year
  enso34years = 1891:((osayears[1] - 3)+i+ny)
  
  # end of last year's ENSO
  nino_pred = nino34ts[which(ninoyear%in%enso34years)]
  
  # last year's natural forcing and amo
  nat_forc  = predictors2[which(predictors2$YEAR%in%predyears), c(4,5) ]
  amo       = predictors2[which(predictors2$YEAR%in%enso34years), 2 ]
  
  # The projected anthropogenic forcing
  anth_forc = predictors2[which(predictors2$YEAR%in%predyears), 6 ]
  
  # Temperature data for the training period
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]
  
  # Bind the exogenous variables together
  xreg=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))
  
  pred = gmstARIMAX(model_code=arimax_001_code, xreg=xreg, hadcrut=hadcrut_train, ny=1)
  
  osamean001[i] = tail(pred$mean,1)
  osaupper001[i] = tail(pred$upper,1)
  osalower001[i] = tail(pred$lower,1)
  osarank001[i] = rank.unc(pred$y[, ncol(pred$y)], obs[i])
}

# Basic prediction error statistics
arimax_mae_001 = mean(abs(obs - osamean001))
arimax_rmse_001 = sqrt(mean((obs - osamean001)^2))

plot(obs, type = 'o')
points(osamean101, type = 'o', col = 'red')
points(osamean100, type = 'o', col = 'blue')
points(osamean001, type = 'o', col = 'orange')

# Are the regressors doing all the work?
# Here is what a straight ARIMA fit to the data looks like

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Set up residuals
  for(t in 1:q) {
  eps[t] <- z[t] - alpha
  }
  # Likelihood
  for (t in (q+1):T) {
  z[t] ~ dnorm(alpha + ar_mean[t] + ma_mean[t], tau)
  ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
  ar_mean[t] <- inprod(phi, z[(t-p):(t-1)])
  eps[t] <- z[t] - alpha - ar_mean[t] - ma_mean[t]
  }
  
  # Priors
  alpha ~ dnorm(0.0,0.01)
  for (i in 1:q) {
  theta[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:p) {
  phi[i] ~ dnorm(0.0,0.01)
  }
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
}
'

gmstARIMA = function(model_code, hadcrut, ny=1, p=1,q=1, d=1){
  
  # Create some predictions off into the future - this time do it within jags
  # A neat trick - just increase T and add on NAs into y!
  T_future = ny # Number of future data points
  real_data_future = with(hadcrut,
                          list(T = nrow(hadcrut) + T_future - d,
                               z = c(diff(hadcrut$Anomaly,
                                          differences = d),
                                     rep(NA,T_future)),
                               q = 1,
                               p = 1))
  
  # Just watch y now
  model_parameters =  c("z")
  
  # Run the model
  real_data_run_future = jags(data = real_data_future,
                              parameters.to.save = model_parameters,
                              model.file=textConnection(model_code),
                              n.chains=4,
                              n.iter=1000,
                              n.burnin=200,
                              n.thin=2)
  
  # Get the future values
  z_all = real_data_run_future$BUGSoutput$sims.list$z
  # If you look at the above object you'll see that the first columns are all identical because they're the data
  z_all_mean = apply(z_all,2,'mean')
  y_all_mean = cumsum(c(hadcrut$Anomaly[1],z_all_mean))
  year_all = c(hadcrut$Year,(max(hadcrut$Year)+1):(max(hadcrut$Year)+T_future))
  
return(list(year_all = year_all, z_all_mean = z_all_mean, y_all_mean = y_all_mean))

}


# setup containers for the predictions.
osamean_arima = rep(NA,15)
osalower_arima = rep(NA,15)
osaupper_arima = rep(NA,15)
osarank_arima = rep(NA,15)

for(i in 1:(length(osayears))){
  trainyears = 1892:((osayears[1] - 2)+i)
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]
  pred = gmstARIMA(model_code, hadcrut=hadcrut_train)
  osamean_arima[i] = tail(pred$y_all_mean,1)
  
}

pdf(file = 'orders.pdf', width = 8, height = 6)
plot(osayears, obs, type = 'o', bty = 'n', col = 'black', pch = 19,
     xlab = 'Year', ylab = 'Anomaly')
points(osayears,osamean_arima, type = 'o', col = 'purple', pch = 19)
points(osayears,osamean303, type = 'o', col = 'pink', pch = 19)
points(osayears,osamean101, type = 'o', col = 'red', pch = 19)
points(osayears,osamean100, type = 'o', col = 'blue', pch = 19)
points(osayears,osamean001, type = 'o', col = 'orange', pch = 19)
legend('topleft', legend = c('obs', 'arima101', 'arimax303', 'arimax101', 'arimax100', 'arimax001'),
       col = c('black', 'purple','pink', 'red', 'blue', 'orange'),
       lty = 'solid',
       pch = 19)
dev.off()

# Need to difference - previous work suggested ARIMA(3,1,3) would be a good fit.

