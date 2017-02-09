# This code repeats Doug's alternative_orders code but uses the forecast package instead to produce predictions

# Clear workspace
rm(list=ls())

# Load in packages
pkgs = c("tidyverse", "forecast")
lapply(pkgs, library, character.only = TRUE)


# Load in data ------------------------------------------------------------

# load data
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/load_gmst_data.R')
# This loads in:
# gmt_predictions - the previous predictions on Global mean temperature including lower and upper estimates and WMO?
# hadcrut - the year and anomaly temperature values from 1850 to 2016
# hadcrut_all - Same as above but monthly? (missing a month)
# hadcrutfile - The string of the web address for the hadcrut data
# nino34 and nino34file - The ENSO monthly values
# nino34ts and ninoyear - the above but (I think) averaged over month with associated year
# predictors1 and predictors2 - two very slightly different versions of forcings for each year. predictors2 seems to be used more often
# rcp and rcpurl - a large data frame continaing regional concentration pathway values - not sure what all of these mean

# Don't know what this function does but looks like some kind of ranking uncertainty?
rank.unc = function(pred, obs){
  # pred is sample of prediction
  ranks = rank(c(obs,pred))
  out = ranks[1] / length(pred)
  out
}

# One step ahead predictions ----------------------------------------------

# One-step-ahead prediction for 1985 - 1999, using auto.arima
# The auto.arima function automatically estimates the best model based on AIC

# Years to include in one step ahead (osa) analysis
osayears = 1985:1999

# Find out which years in the hadcrut data are in the osa years
obs = hadcrut$Anomaly[which(hadcrut$Year%in%osayears)]

# setup containers for the predictions.
osamean202 = rep(NA,15)
osalower202 = rep(NA,15)
osaupper202 = rep(NA,15)
osarank202 = rep(NA,15)

# Start looping through the osa years
for(i in 1:(length(osayears))){
  # Number of years ahead (always 1)
  ny = 1

  # Get the training years
  trainyears = 1892:((osayears[1] - 2)+i)

  # And the year to predict
  predyears = 1892:((osayears[1] - 2)+i+ny)

  # offset enso by a year (why?)
  enso34years = 1891:((osayears[1] - 3)+i+ny)

  # end of last year's ENSO
  nino_pred = nino34ts[which(ninoyear%in%enso34years)]

  # last year's natural forcing and amo
  nat_forc = predictors2[which(predictors2$YEAR%in%predyears), c(4,5) ]
  amo = predictors2[which(predictors2$YEAR%in%enso34years), 2 ]

  # The projected anthropogenic forcing
  anth_forc = predictors2[which(predictors2$YEAR%in%predyears), 6 ]

  # Temperature data for the training period
  hadcrut_train = hadcrut[which(hadcrut$Year%in%trainyears), ]

  # Bind the exogenous variables together
  xreg=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))

  # Now run auto.arima on these data
  y = ts(hadcrut_train$Anomaly, start = min(hadcrut_train$Year), end = max(hadcrut_train$Year))
  # curr_model = auto.arima(y = y,
  #                         xreg = xreg[1:length(y),])
  # Or use standard arima if you want to fix the parameters
  curr_model = Arima(x = y,
                     xreg = xreg[1:length(y),],
                     order = c(2, 0, 2))
  #plot(curr_model) - optional

  # Forecast into the future
  pred = forecast(curr_model, h = 1, xreg = xreg[nrow(xreg),,drop=FALSE])
  plot(pred)

  # Look at the model it chose
  print(osayears[i])
  print(pred$method)

  # Get the means and upper/lower 95% uncertainty interval values
  osamean202[i] = pred$mean
  osaupper202[i] = pred$upper[2]
  osalower202[i] = pred$lower[2]
  # Still not sure what this line is doing
  osarank202[i] = rank.unc(pred$x[length(pred$x)], obs[i])
}
# Looks to me like 1,0,2 or 2,0,1 usually wins so suggest 2,0,2 to be safe

# Basic prediction error statistics
arimax_mae_202 = mean(abs(obs - osamean202))
arimax_rmse_202 = sqrt(mean((obs - osamean202)^2))

# Create a plot of how the forecasts do over time
pdf(file = 'orders_AP.pdf', width = 8, height = 6)
plot(osayears, obs, type = 'o', bty = 'n', col = 'black', pch = 19,
     xlab = 'Year', ylab = 'Anomaly')
points(osayears,osamean202, type = 'o', col = 'red', pch = 19)
legend('topleft', legend = c('obs', 'arima202'),
       col = c('black', 'red'),
       lty = 'solid',
       pch = 19)
dev.off()

# Create a nicer ggplot
df1 = data.frame(Year = rep(osayears,4),
                 Anomaly = c(obs,osamean202, osalower202, osaupper202),
                 Type2 = c(rep('Observations', length(osayears)),
                          rep('ARIMA (2,0,2) mean', length(osayears)),
                          rep('ARIMA (2,0,2) lower', length(osayears)),
                          rep('ARIMA (2,0,2) upper', length(osayears))),
                 Type = c(rep('Observations', length(osayears)),
                           rep('ARIMA (2,0,2)', 3*length(osayears))))
ggplot(df1, aes(x = Year, y = Anomaly, group = Type2, colour = Type)) +
  geom_line(aes(linetype = Type)) +
  theme_bw() +
  theme(legend.position = 'top')
