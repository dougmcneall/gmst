# load_gmst_data.R
# Load up all the data we might need to build an ARIMAX
# model for global temperature prediction.
# Doug McNeall

gmt_predictions = read.table('https://raw.githubusercontent.com/dougmcneall/gmst/master/gmt_predictions.txt',
                             head=TRUE)

hadcrutfile = 'http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.5.0.0.annual_ns_avg.txt'
hadcrut_all = read.table(hadcrutfile)
hadcrut = hadcrut_all[ , 1:2]
colnames(hadcrut) = c('Year', 'Anomaly')

predictors1 = read.table('/Users/dougmcneall/Documents/work/R/Rcode/useful/predictors1.txt',
                         header=TRUE)
predictors2 = read.table('/Users/dougmcneall/Documents/work/R/Rcode/useful/predictors2.txt', header=TRUE)

# Nino3.4 from 1870 to 2016.
nino34file='https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/nino34.long.anom.data'
nino34 = read.table(file=nino34file, skip=1, nrows=147, na.strings='-99.99')

# just taking the mean of the last 6 months at the moment
nino34ts = apply(nino34[,8:13],1, mean,na.rm = TRUE )
ninoyear = nino34[,1]

rcpurl = 'http://www.pik-potsdam.de/~mmalte/rcps/data/RCP85_MIDYEAR_RADFORCING.DAT'
rcp = read.table(file = rcpurl, skip = 59, header=TRUE)