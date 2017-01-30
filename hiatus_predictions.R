# hiatus_predictions.R
# Predicting the "hiatus" in global mean surface temperature
# using an ARIMAX model
# Doug McNeall

library(R2jags)

# load ARIMAX prediction function
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/gmstARIMAX.R')

# load data
source('https://raw.githubusercontent.com/dougmcneall/gmst/master/load_gmst_data.R')


# Use temperature and forcing data up to 1999 to train the statistical
# model, and predict the rest of the data, keeping one factor at a time flat.

ny = 17
trainyears = 1892:1999
predyears = 1892:(1999+ny)
# offset enso etc. by a year
enso34years = 1891:(1998+ny)

# This is the 'base' prediction (for 2017), with everything lined up right

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


