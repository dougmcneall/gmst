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
xreg_all=cbind(nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))

pred_all = gmstARIMAX(model_code=arimax_code, xreg=xreg_all, hadcrut=hadcrut_train, ny=ny)

# --------------------------------------------------------
# Replace the elements of xreg_all one-at-a-time
# --------------------------------------------------------

# Replace ENSO during prediction years with neutral (zero)
no_nino_pred = c(nino34ts[which(ninoyear%in%trainyears)], rep(0,ny))
xreg_no_nino = cbind(no_nino_pred, as.matrix(cbind(nat_forc,amo,anth_forc)))
pred_no_nino = gmstARIMAX(model_code=arimax_code, xreg=xreg_no_nino,
                          hadcrut=hadcrut_train, ny=ny)

# Replace Volcanic forcing with last value
volc_train = predictors2[which(predictors2$YEAR%in%trainyears),4]
volc_pred = rep(tail(volc_train,1),ny)
volc = c(volc_train,volc_pred)

xreg_no_volc = xreg_all
xreg_no_volc[, 'VOLC'] = volc

pred_no_volc = gmstARIMAX(model_code=arimax_code, xreg=xreg_no_volc,
                          hadcrut=hadcrut_train, ny=ny)


# Replace Solar forcing with last value
solar_train = predictors2[which(predictors2$YEAR%in%trainyears), 5]
solar_pred = rep(tail(solar_train,1),ny)
solar = c(solar_train, solar_pred)

xreg_no_solar = xreg_all
xreg_no_solar[, 'SOLAR'] = solar

pred_no_solar = gmstARIMAX(model_code=arimax_code, xreg=xreg_no_solar,
                          hadcrut=hadcrut_train, ny=ny)


# Replace AMO with last value
amo_train =  head(amo, -ny)
amo_pred = rep(tail(amo_train,1),ny)
amo_flat = c(amo_train, amo_pred)
xreg_no_amo = xreg_all
xreg_no_amo[, 'amo'] = amo_flat

pred_no_amo = gmstARIMAX(model_code=arimax_code, xreg=xreg_no_amo,
                           hadcrut=hadcrut_train, ny=ny)

# replace anthropogenic forcing
anth_forc_train = predictors2[which(predictors2$YEAR%in%trainyears), 6 ]
anth_forc_pred = rep(tail(anth_forc_train,1), ny)
anth_forc_flat = c(anth_forc_train, anth_forc_pred)

xreg_no_anth_forc = xreg_all
xreg_no_anth_forc[, 'anth_forc'] = anth_forc_flat

pred_no_anth_forc = gmstARIMAX(model_code=arimax_code, xreg=xreg_no_anth_forc,
                         hadcrut=hadcrut_train, ny=ny)

# -------------------------------------------------------
# Plot the results
# -------------------------------------------------------


ix = which(hadcrut$Year%in%predyears)
obs_col = 'black'
all_col = 'red'
enso_col = 'skyblue'
solar_col = 'orange'
volc_col = 'grey'
amo_col = 'purple'
anth_col = 'blue'

pdf(file='hiatus_forecast.pdf', width=5, height=5)
lwd = 2.5
par(mar = c(5,4,2,1), las=1)
plot(hadcrut$Year[ix],hadcrut$Anomaly[ix], ylim=c(0.3,0.8), xlim=c(1995, 2022),
     type='n', bty='n', lwd=lwd, ylab='Temperature Anomaly (C)', xlab='year')
lines(hadcrut$Year[ix], pred_all$mean, col=all_col, lwd=lwd)
lines(hadcrut$Year[ix], pred_no_nino$mean, col=enso_col, lwd=lwd)
lines(hadcrut$Year[ix], pred_no_solar$mean, col=solar_col, lwd=lwd)
lines(hadcrut$Year[ix], pred_no_volc$mean, col=volc_col, lwd=lwd)
lines(hadcrut$Year[ix], pred_no_amo$mean, col=amo_col, lwd=lwd)
lines(hadcrut$Year[ix], pred_no_anth_forc$mean, col=anth_col, lwd=lwd)
lines(hadcrut$Year[ix],hadcrut$Anomaly[ix], col = obs_col, lwd = lwd)

text(2015, 0.78 , 'HADCRUT4', col=obs_col, pos=2)
text(2013,0.67, 'No ENSO', col=enso_col, pos=2)
text(2005, 0.57, 'All', col=all_col, pos=2)
text(2016, tail(pred_no_anth_forc$mean,1), 'No Anthro', col=anth_col, pos=4)
text(2015, 0.63, 'No AMO', col=amo_col, pos=4)
text(2016, 0.8, 'No solar', col=solar_col, pos=4)
text(2016, 0.75, 'No Volc', col=volc_col, pos=4)
dev.off()


# How much difference did each factor make?
obs_target = tail(pred_all$mean,ny) 
enso_effect = obs_target - tail(pred_no_nino$mean,ny)
solar_effect = obs_target - tail(pred_no_solar$mean,ny)
volc_effect = obs_target - tail(pred_no_volc$mean,ny)
amo_effect = obs_target - tail(pred_no_amo$mean,ny)
anth_effect = obs_target - tail(pred_no_anth_forc$mean,ny)
missing_effect = tail(hadcrut$Anomaly, ny) - obs_target

target_years = 2000:2016
pdf(file='effect_size.pdf', width=5, height=5)
par(las = 1)
lwd = 2.5
plot(target_years, enso_effect, type='l', 
     col=enso_col, lwd = lwd, bty = 'n', ylim = c(-0.15, 0.25),
     ylab = 'Effect size (C)', xlab = 'Year')
abline(h = 0)
lines(target_years, solar_effect, col=solar_col, lwd=lwd)
lines(target_years, volc_effect, col=volc_col, lwd = lwd)
lines(target_years, amo_effect, col=amo_col, lwd = lwd)
lines(target_years, anth_effect, col=anth_col, lwd = lwd)
lines(target_years, missing_effect, col=all_col, lwd = lwd)

text(2011,-0.14, 'ENSO', col=enso_col)
text(2004, -0.09, 'Missing', col=all_col)
text(2008, 0.2, 'Anthro', col=anth_col, pos=4)
text(2011, 0.1, 'AMO', col=amo_col, pos=4)
text(2006, 0.04, 'Solar', col=solar_col, pos=4)
text(2011, 0.02, 'Volc', col=volc_col, pos=4)
dev.off()


