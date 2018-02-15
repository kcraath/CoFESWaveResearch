#=====================================================================================#
# PURPOSE : Application of Wavelet-ARIMA hybrid model for forecasting time series     #
# AUTHOR  : Ranjit Kumar Paul and Sandipan Samanta                                    #
# DATE    : 24 October, 2017                                                          #
# VERSION : Ver 0.1.0                                                                 #
# Adaption AUTHOR: Kim Raath                                                      #
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

#To do: (1) save predictive intervals and plot alongside the predictions and then also (2) RMSE, MAE and formal tests such as Diebold & Mariano.

#https://insightr.wordpress.com/2017/11/09/formal-ways-to-compare-forecasting-models-rolling-windows/
#Do a one step-ahead forecast - similar to TS Class for babies forecast (see these plots)

#One 

#Wavelet
#Set up the correct data frame
###################### I added returns here
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

errorOne <- r1[forecastlength:length(r1[,2]),2] - as.numeric(newTSOne[forecastlength:length(r1[,2])])
MSE_One <- mean(errorOne^2)
Sd_UP_One <- newTSOne_Upper
Sd_DOWN_One <- newTSOne_Lower

#png("/Users/kcr2/Dropbox/Kim Folder/Resume and Research/Research/Projects/Paper/Plots/WaveletForecast__One.png",width = 680, height = 480)
plot(ts(r1[,2]), xlab = "Date", ylab = "Price", main = "One-Step ahead Wavelet Arima Forecast")
#axis(1, at = r1[,1], labels = date1, las = 1)
lines(ts(newTSOne), col = c("red"))
lines(as.numeric(Sd_UP_One), col = "lightblue", lty = 2)
lines(as.numeric(Sd_DOWN_One), col = "lightblue", lty = 2)
text(locator(1),"MSE One = 0.95", 4)
#dev.off()

splitTSOne2 <- as.xts(r1[,2], order.by = date12) 
colnames(splitTSOne2) <- "Original Series"
newTSOne2 <- as.xts(as.numeric(newTSOne), order.by = date12)
colnames(newTSOne2) <- "Forecast"
Sd_UP_One <- as.xts(as.numeric(Sd_UP_One), order.by = date12)
colnames(Sd_UP_One) <- "Upper 95% Prediction Interval"
Sd_DOWN_One <- as.xts(as.numeric(Sd_DOWN_One), order.by = date12)
colnames(Sd_DOWN_One) <- "Lower 95% Prediction Interval"
chart.TimeSeries(cbind(splitTSOne2, newTSOne2, Sd_UP_One, Sd_DOWN_One), lty = c(2,1,1,1), colorset = rainbow10equal, legend.loc = "topright", ylab = "Prices", main = "One-Step ahead Wavelet Arima Forecast | One", las = 3,lwd =1)
text(locator(1),"MSE One = 0.33", 4)

#Two 

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

#Calculations

errorTwo <- r2[1305:length(r2[,2]),2] - as.numeric(newTSTwo[1305:2609])
MSE_Two <- mean(errorTwo^2)
Sd_UP_Two <- newTSTwo_Upper
Sd_DOWN_Two <- newTSTwo_Lower


#Do a one step-ahead forecast - similar to TS Class for babies forecast
splitTSTwo2 <- as.xts(r2[,2], order.by = date1) 
colnames(splitTSTwo2) <- "Original Series"
newTSTwo2 <- as.xts(as.numeric(newTSTwo), order.by = date1)
colnames(newTSTwo2) <- "Forecast"
Sd_UP_Two <- as.xts(as.numeric(newTSTwo_Upper), order.by = date1)
colnames(Sd_UP_Two) <- "One Std Deviation Up"
Sd_DOWN_Two <- as.xts(as.numeric(newTSTwo_Lower), order.by = date1)
colnames(Sd_DOWN_Two) <- "One Std Deviation Down"
chart.TimeSeries(cbind(splitTSTwo2[2305:2609], newTSTwo2[2305:2609], Sd_UP_Two[2305:2609], Sd_DOWN_Two[2305:2609]), lty = c(2,1,1,1), colorset = rainbow10equal, legend.loc = "bottomright", ylab = "Prices", main = "One-Step ahead Wavelet Arima Forecast | Two", las = 3, lwd =1)
text(locator(1),"MSE Two = 0.063", 4)
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
SumsTS <- as.xts(nTS, order.by = date1[1:1304])
OriginalTS <- as.xts(newTSOne, order.by = date1[1:1304])
plot.ts(cbind(nTS, newTSOne), main = "Reconstructed Coefficient Sums and Original Time Series", ylab = c("Sums", "Original"))
#SideShow

