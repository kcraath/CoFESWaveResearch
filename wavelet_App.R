#=====================================================================================#
# PURPOSE : Implementing various wavelet tools using quantmod data                    #
# AUTHOR  : Kim Raath                                                                 #
# DATE    : 02 February, 2018                                                         #
# VERSION : Ver 0.1.0                                                                 #
#=====================================================================================#
require(vars)
require(biwavelet)
require(timeSeries)
require(glmnet)
require(quantmod)
require(PerformanceAnalytics)
require(WaveletComp)
require(xts)
require(dygraphs)
require(RColorBrewer)
require(TSA)
require(cointReg)
require(e1071)
library(reshape)
library(ggplot2)
library(dplyr)
library(colorRamps)
library(scales)
library(wavScalogram)    ##### This section will not work yet  - package not available ####

#Specify tickers
tickers <- c("WTI","HTG")
getSymbols('DCOILWTICO',src='FRED') #src: from FRED 
WTI <- DCOILWTICO['2017-01-01::2018-01-16'] #Dates need to be specified
getSymbols("HTG",src="google", from = "2017-01-01",to = "2018-01-16") #src: from google
AllPrices <- do.call(merge, lapply(tickers, function(x) get(x)))

#There are missing values in the series and pound to usd conversion
AllPrices$HTG.Close <- AllPrices$HTG.Close/10*0.75 #average usd gbx exchange over the past year
price.pair <- cbind(AllPrices$DCOILWTICO ,AllPrices$HTG.Close) #table
price.pair <- removeNA(price.pair)
colnames(price.pair) <- c("WTI", "HTG")
chart.TimeSeries(price.pair, colorset = rainbow10equal, legend.loc = "bottomright", ylab = "Prices", main = "WTI and HTG Prices", las = 3, 
                 lwd = 0.5,cex.labels = 1, cex.main = 1, cex.legend = 1, pch = " ")

#Wavelet
#Set up the correct data frame
###################### I added returns here
rWTI <- as.data.frame(returns(price.pair$WTI))
rHTG <- as.data.frame(returns(price.pair$HTG))

#Retrieve the dates
date12 <- index(price.pair)[-1]
#save Prices in Matrix
r1 <- cbind(1:(length(price.pair$WTI)), rWTI$WTI[1: length(price.pair$WTI)])
r2 <- cbind(1:(length(price.pair$HTG)), rHTG$HTG[1: length(price.pair$HTG)])

#Wavelet Coherence Plot
wtc.r12=wtc(r1[-1,], r2[-1,], quiet = TRUE, nrands = 100)
#wtc.r12$xaxis <- date12
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 5, 5) + 0.1)
plot(wtc.r12, plot.cb=TRUE, plot.phase=TRUE, xlab = "Days", ylab = "Period (Days)", cex = 1, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2)
n = length(r1[, 1])
abline(v = seq(250, n, 250), h = 1:16, col = "brown", lty = 1, lwd = 1)
title("WSC: WTI and HTG Prices", cex = 1)

# WSD - Package not avail ################

# Parameters
wname <- "MORLET"
dt <- 1
compensated <- 1
delta_t <- 1
commutative <- TRUE
e_density <- TRUE
nt <- length(r1)/2 # number of time points
t <- 1:nt # time vector

nrand <- 100  # MonteCarlo repetitions (for significance contours)
windowrad <- floor(nt/60); # % time radius for windowed scalogram (width)
rdist <- floor(nt/300); # % Scale radius for distance (height)

# Defining the scales (Torrence and Compo's way)
s0 <- 2*dt
Dj <- 12
waverad <- 3 # Morlet wavelet radius
#J <- log2((t[nt]-t[1]-2*windowrad)/(2*waverad*s0*dt))
smax <- (nt-1-2*windowrad)/(2*waverad)
scales <- c(s0,smax, Dj)
#scales <- s0*2^((0:(J*Dj))/Dj)
nesc <- length(scales)

system.time(wsd <-wsd(signal1 = r1[-1,2], signal2 = r2[-1,2], scales = scales, delta_t = 1,
                      windowrad = windowrad, rdist = rdist, mc_nrand = nrand,
                      wname = "MORLET", parallel = TRUE, makefigure = TRUE))

#The lighter the color the higher the similarity of change. Scale is level of similarity. 

#=====================================================================================#
# PURPOSE : Application of Wavelet-ARIMA hybrid model for forecasting time series     #
# AUTHOR  : Ranjit Kumar Paul and Sandipan Samanta                                    #
# DATE    : 24 October, 2017                                                          #
# VERSION : Ver 0.1.0                                                                 #
# Adaption AUTHOR: Kim Raath                                                          #
# MODIFICATION: 06 December 2017                                                      #
#=====================================================================================#

require(fracdiff)
require(forecast)
require(stats)
require(wavelets)

#---------------------------------------------------------------------------------------#
# Computing Wavelet Coefficients using MODWT algorithm using haar filter                #
#---------------------------------------------------------------------------------------#

WaveletFitting <- function(ts,Wvlevels,bndry,FFlag)
{
  mraout <- wavelets::modwt(ts, filter='haar', n.levels=Wvlevels,boundary=bndry, fast=FFlag)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
}


WaveletFittingarma<- function(ts,Waveletlevels,boundary,FastFlag,MaxARParam,MaxMAParam,NForecast)
  
{
  WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL
  AllWaveletForecast_Lower <- NULL
  AllWaveletForecast_Upper <- NULL
  #-----------------------------------------------------------#
  # Fitting of ARIMA model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    #acf(ts)
    WaveletARMAFit <- forecast::auto.arima(x=as.ts(ts), d=NA, D=NA, max.p=MaxARParam, max.q=MaxMAParam,stationary=FALSE,
                                           seasonal=FALSE,ic=c("aicc"), allowdrift=FALSE, allowmean=TRUE,stepwise = TRUE)
    print(WaveletARMAFit)
    WaveletARIMAForecast <- forecast::forecast(WaveletARMAFit,h=NForecast)
    plot(WaveletARIMAForecast)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletARIMAForecast$mean))
    ##### These are only for h=1 ######
    AllWaveletForecast_Lower <- cbind(AllWaveletForecast_Lower, as.matrix(WaveletARIMAForecast$lower[2])) #lower 95% prediction interval
    AllWaveletForecast_Upper <- cbind(AllWaveletForecast_Upper, as.matrix(WaveletARIMAForecast$upper[2])) #Upper 95% prediction interval
  }
  Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
  Finalforecast_Lower <- rowSums(AllWaveletForecast_Lower,na.rm = T)
  Finalforecast_Upper <- rowSums(AllWaveletForecast_Upper,na.rm = T)
  return(list(Finalforecast, Finalforecast_Lower, Finalforecast_Upper))
}

############################################################

#One 

#Wavelet
#Set up the correct data frame
###################### Can add returns here instead of prices
rWTI <- as.data.frame((price.pair$WTI))
rHTG <- as.data.frame((price.pair$HTG))

#Retrieve the dates
date12 <- index(price.pair)
#save Prices in Matrix
r1 <- cbind(1:(length(price.pair$WTI)), rWTI$WTI[1: length(price.pair$WTI)])
r2 <- cbind(1:(length(price.pair$HTG)), rHTG$HTG[1: length(price.pair$HTG)])

forecastlength <- floor(length(r1[,2])/2)
splitTSOne <- r1[1:forecastlength,2]
Waveletlevels <- floor(log(length(splitTSOne))) # to obtain the maximum level for wavelet decomposition
newTSOne <- splitTSOne
newTSOne_Lower <- splitTSOne
newTSOne_Upper <- splitTSOne

#Testing
wfit <- wavelets::modwt(newTSOne, filter='haar', n.levels=Waveletlevels,boundary='periodic', fast = T)
wfit2 <- WaveletFitting(ts(newTSOne), Waveletlevels, 'periodic', TRUE)
plot.ts(wfit2$WaveletSeries, main = "Maximal Decomp Wavelet Coefficients and max Scaling Coefficient")
wfit3 <- WaveletFittingarma(ts(splitTSOne), Waveletlevels=Waveletlevels, boundary='periodic', FastFlag = T, MaxARParam=5,MaxMAParam=5,NForecast=1)
par(mfrow=c(1,1))

#Simulation for whole series 
system.time(for(i in 1:ceiling(length(r1[,2])/2)){  #
  WS <- WaveletFittingarma(ts(splitTSOne), Waveletlevels=Waveletlevels, boundary='periodic', FastFlag = T, MaxARParam=5,MaxMAParam=5,NForecast=1)
  newTSOne <- append(newTSOne,WS[1])
  newTSOne_Lower <- append(newTSOne_Lower,WS[2])
  newTSOne_Upper <- append(newTSOne_Upper,WS[3])
  forecastlength = forecastlength + 1 # Expanding window forecastlength + 1
  splitTSOne <- r1[i:forecastlength,2] # Rolling window is when we replace 1 with i
  i = i + 1
}
)

#TS One      
#Calculations for MSE and SE                                   
errorOne <- r1[forecastlength:length(r1[,2]),2] - as.numeric(newTSOne[forecastlength:length(r1[,2])])
MSE_One <- mean(errorOne^2)
Sd_UP_One <- newTSOne_Upper
Sd_DOWN_One <- newTSOne_Lower

splitTSOne2 <- as.xts(r1[,2], order.by = date12) 
colnames(splitTSOne2) <- "Original Series"
newTSOne2 <- as.xts(as.numeric(newTSOne), order.by = date12)
colnames(newTSOne2) <- "Forecast"
Sd_UP_One <- as.xts(as.numeric(Sd_UP_One), order.by = date12)
colnames(Sd_UP_One) <- "Upper 95% Prediction Interval"
Sd_DOWN_One <- as.xts(as.numeric(Sd_DOWN_One), order.by = date12)
colnames(Sd_DOWN_One) <- "Lower 95% Prediction Interval"
chart.TimeSeries(cbind(splitTSOne2, newTSOne2, Sd_UP_One, Sd_DOWN_One), lty = c(2,1,1,1), colorset = rainbow10equal, legend.loc = "topright", ylab = "Prices", main = "One-Step ahead Wavelet Arima Forecast | One", las = 3,lwd =1)
text(locator(1),"MSE One = 0.33", 4)   #Subject to change given value above

#TS Two 
forecastlength <- floor(length(r2[,2])/2)
splitTSTwo <- r2[1:forecastlength,2]
Waveletlevels <- floor(log(length(splitTSTwo)))
newTSTwo <- splitTSTwo
newTSTwo_Lower <- splitTSTwo
newTSTwo_Upper <- splitTSTwo
system.time(
  for(i in 1:ceiling(length(r2[,2])/2)){
    WS <- WaveletFittingarma(ts(splitTSTwo), Waveletlevels=Waveletlevels, boundary='periodic', FastFlag = T, MaxARParam=5,MaxMAParam=5,NForecast=1)
    newTSTwo <- append(newTSTwo,WS[1])
    newTSTwo_Lower <- append(newTSTwo_Lower,WS[2])
    newTSTwo_Upper <- append(newTSTwo_Upper,WS[3])
    forecastlength = forecastlength + 1
    splitTSTwo <- r2[i:forecastlength,2]
    i = i + 1
  }
)

#Calculations for MSE and SE
errorTwo <- r2[forecastlength:length(r2[,2]),2] - as.numeric(newTSTwo[forecastlength:length(r2[,2])])
MSE_Two <- mean(errorTwo^2)
Sd_UP_Two <- newTSTwo_Upper
Sd_DOWN_Two <- newTSTwo_Lower

#Do a one step-ahead forecast
splitTSTwo2 <- as.xts(r2[,2], order.by = date1) 
colnames(splitTSTwo2) <- "Original Series"
newTSTwo2 <- as.xts(as.numeric(newTSTwo), order.by = date1)
colnames(newTSTwo2) <- "Forecast"
Sd_UP_Two <- as.xts(as.numeric(newTSTwo_Upper), order.by = date1)
colnames(Sd_UP_Two) <- "One Std Deviation Up"
Sd_DOWN_Two <- as.xts(as.numeric(newTSTwo_Lower), order.by = date1)
colnames(Sd_DOWN_Two) <- "One Std Deviation Down"
chart.TimeSeries(cbind(splitTSTwo2, newTSTwo2, Sd_UP_Two, Sd_DOWN_Two), lty = c(2,1,1,1), colorset = rainbow10equal, legend.loc = "bottomright", ylab = "Prices", main = "One-Step ahead Wavelet Arima Forecast | Two", las = 3, lwd =1)
text(locator(1),"MSE Two = 0.063", 4)   #Subject to change given value above
par(mfrow=c(1,1))

png("/Users/kcr2/Dropbox/Kim Folder/Resume and Research/Research/Projects/Paper/Plots/WaveletForecast__Two.png",width = 680, height = 480)
plot(ts(splitTSTwo), xlab = "Date", ylab = "Price", main = "One-Step ahead Wavelet Arima Forecast")
lines(ts(newTSTwo), col = c("red"))
dev.off()
newTSTwo <- newTSTwo
splitTSTwo <- splitTSTwo

save.image(file = "Two_One")
testing_table <- cbind(do.call(cbind,wfit2$WVSeries@W), wfit2$WVSeries@V[[Waveletlevels]])
nTS <- rowSums(testing_table)
SumsTS <- as.xts(nTS, order.by = date1[1:forecastlength])
OriginalTS <- as.xts(newTSOne, order.by = date1[1:forecastlength])
plot.ts(cbind(nTS, newTSOne), main = "Reconstructed Coefficient Sums and Original Time Series", ylab = c("Sums", "Original"))

