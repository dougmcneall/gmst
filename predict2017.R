# predict2017.R
# Predict 2017, using training data to 2016
# Doug McNeall

library(R2jags)

# load ARIMAX prediction function
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/gmstARIMAX.R')

# load data
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/load_gmst_data.R')


# ---------------------------------------------------------------------
# GMST Prediction for 2017 
# ---------------------------------------------------------------------

p = 2
q = 2

trainyears = 1892:2016
predyears =  1892:2017
# Offset the ENSO years, so that we get the last 6 months of last year's
enso34years = 1891:2016

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

pred = gmstARIMAX(model_code=arimax_code, xreg=xreg, hadcrut=hadcrut_train,
                  ny=1, p=p, q=q)

pdf(file='predict2017.pdf', width = 6, height = 5)
par(las = 1)
plot(head(tail(pred$year,30),29), head(tail(pred$mean, 30),29), 
     axes = FALSE, xlab = 'Year', ylab = 'Anomaly',
     type = 'o', pch = 19, xlim = c(1988, 2018))
points(tail(pred$year,1), tail(pred$mean,1), pch = 19)
segments(x0 = tail(pred$year,1), y0 = tail(pred$lower,1), x1 = tail(pred$year,1),y1 = tail(pred$upper,1))
axis(1)
axis(2)
dev.off()



