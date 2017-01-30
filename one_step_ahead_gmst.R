# one_step_ahead_gmst.R
# Do a one-step-ahead forecast of global mean surface 
# temperature using an ARIMAX model, and compare it against
# the Met Office prediction

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
# One-step-ahead prediction for 2000 - 2016
# ---------------------------------------------------------------------

# setup containers for the predictions.
osamean = rep(NA,17)
osalower = rep(NA,17)
osaupper = rep(NA,17)
osarank = rep(NA,17)
osayears = 2000:2016
obs = tail(hadcrut$Anomaly,17)

for(i in 1:17){
  ny = 1
  trainyears = 1892:(1998+i)
  predyears = 1892:(1998+i+ny)
  # offset enso by a year
  enso34years = 1891:(1997+i+ny)
  
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
  
  pred = gmstARIMAX(model_code=arimax_code, xreg=xreg, hadcrut=hadcrut_train, ny=1)
  
  osamean[i] = tail(pred$mean,1)
  osaupper[i] = tail(pred$upper,1)
  osalower[i] = tail(pred$lower,1)
  osarank[i] = rank.unc(pred$y[, ncol(pred$y)], obs[i])
}

# Basic prediction error statistics
mo_mae = mean(abs(gmt_predictions$MEAN - gmt_predictions$WMO))
mo_rmse = sqrt(mean((gmt_predictions$MEAN - gmt_predictions$WMO)^2 ))

arimax_mae = mean(abs(obs - osamean))
arimax_rmse = sqrt(mean((obs - osamean)^2))

# Plot the Met Office and ARIMAX year-ahead prediction, along with the
# observed data.
dev.new()
par(las=1, mar = c(5,4,2,1), mfrow = c(2,1))

wmo.col = 'darkred'
had.col = 'darkblue'
mopred.col ='tomato'
arimax.col = 'skyblue2'

plot(osayears, obs, type = 'o', pch = 19, ylim = c(0.2,1.1),
     col = had.col, bty = 'n',
     xlab = 'year', ylab='global mean surface temperature anomaly')
points(osayears, osamean, pch=19, col = arimax.col)
points(gmt_predictions$YEAR, gmt_predictions$WMO+0.2,
       col=wmo.col, pch=19, type = 'o')
points(gmt_predictions$YEAR, gmt_predictions$MEAN+0.2, 
       col= mopred.col, pch = 19)

legend('topleft', c('WMO (offset)','Met Office', 'HadCRUT4', 'ARIMAX'),
       col=c('darkred', 'tomato', 'darkblue', 'skyblue2'), pch = 19,
       bty = 'n', horiz = TRUE)

# error through time
plot(osayears, osamean - obs, ylim = c(-0.2, 0.2),
     col = arimax.col, pch=19, bty='n',
     xlab = 'year', ylab='forecast error')
points(osayears, gmt_predictions$MEAN - gmt_predictions$WMO, col='tomato', pch=19)
abline(h=0, col = 'darkgrey')
segments(x0 = osayears, y0 = 0, col = arimax.col,
         x1 = osayears, y1 = osamean - obs)
segments(x0 = osayears, y0 = 0,
         x1 = osayears, y1 = gmt_predictions$MEAN - gmt_predictions$WMO,
         col = 'tomato')
legend('topleft', c('Met Office', 'ARIMAX'), col=c(mopred.col, arimax.col), pch=19,
       bty = 'n')

# x-y plot of predicted vs observed
par(las = 1)
plot(obs,osamean, pch=19, col=arimax.col,
     xlab='measured', ylab='predicted', bty = 'n', xlim = c(0.2,0.9), ylim = c(0.2,0.9) )
segments(x0 = obs, y0 = osalower, x1 = obs, y1 = osaupper, col = arimax.col)

points(gmt_predictions$WMO, gmt_predictions$MEAN, pch = 19, col = mopred.col)
segments(x0 = gmt_predictions$WMO, y0 = gmt_predictions$LOWER, 
         x1 = gmt_predictions$WMO, y1 = gmt_predictions$UPPER,
          col = mopred.col)
abline(0,1)
legend('topleft', c('Met Office', 'ARIMAX'), col=c('tomato', 'skyblue2'), pch=19,
       bty = 'n')

