#This code can be used to reproduce the forecasts of the M4 Competition STATISTICAL Benchmarks and evaluate their accuracy

library(forecast) #Requires v8.2
library(DTScanF)

#library(foreach)
#library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
#Datasets info
# Daily (4227) | fh=14, frq=1
# Hourly (414) | fh=48, frq=24  DONE!
# Monthly (48k) | fh=18, frq=12
# Yearly (23k) | fh=6, frq=1
# Quarterly (24k) | fh=8, frq=4 
# Weekly (359) | fh=13, frq=1 DONE!
#################################################################################
fh <- 48 #The forecasting hori zon examined
frq <- 24 #The frequency of the data

read_data <- function(path, frequency=F){
  df <- read.csv(path, header=T, row.names='V1')
  df <- t(df)
  df <- lapply(seq_len(ncol(df)), function(i) df[,i])
  return(df)
}

#dataset <- 'Weekly'
#################################################################################

smape_cal <- function(outsample, forecasts){
  #Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts){
  #Used to estimate MASE
  frq <- frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
}

naive_seasonal <- function(input, fh){
  #Used to estimate Seasonal Naive
  frcy <- frequency(input)
  frcst <- naive(input, h=fh)$mean 
  if (frcy>1){ 
    frcst <- head(rep(as.numeric(tail(input,frcy)), fh), fh) + frcst - frcst
  }
  return(frcst)
}

Theta.classic <- function(input, fh){
  #Used to estimate Theta classic
  
  #Set parameters
  wses <- wlrl<-0.5 ; theta <- 2
  #Estimate theta line (0)
  observations <- length(input)
  xt <- c(1:observations)
  xf <- c((observations+1):(observations+fh))
  train <- data.frame(input=input, xt=xt)
  test <- data.frame(xt = xf)
  
  estimate <- lm(input ~ poly(xt, 1, raw=TRUE))
  thetaline0In <- as.numeric(predict(estimate))
  thetaline0Out <- as.numeric(predict(estimate,test))
  
  #Estimate theta line (2)
  thetalineT <- theta*input+(1-theta)*thetaline0In
  sesmodel <- ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  
  #Theta forecasts
  forecastsIn <- (thetaline2In*wses)+(thetaline0In*wlrl)
  forecastsOut <- (thetaline2Out*wses)+(thetaline0Out*wlrl)
  
  #Zero forecasts become positive
  for (i in 1:length(forecastsOut)){
    if (forecastsOut[i]<0){ forecastsOut[i]<-0 }
  }
  
  output=list(fitted = forecastsIn, mean = forecastsOut,
              fitted0 = thetaline0In, mean0 = thetaline0Out,
              fitted2 = thetaline2In, mean2 = thetaline2Out)
  
  return(output)
}

SeasonalityTest <- function(input, ppy){
  #Used to determine whether a time series is seasonal
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )
    
    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }
  
  return(test_seasonal)
}


dtsf <- function(ts, poli, best, window, forecast, reg = FALSE, aggregation = 'median') {
  
  ts_reg <- ts  # If reg = FALSE, ts_reg remains equal to ts
  
  # Regularization
  if(reg) {
    for(i in 2:(length(ts)-1)) {
      ts_reg[i] <- 0.95*ts[i] + 0.025*ts[i-1] + 0.025*ts[i+1]
    }
  }
  
  # Handling polynomial order
  wind <- cbind(1, ts_reg)
  while(dim(wind)[2]-1!=poli) wind <- cbind(wind, ts_reg^dim(wind)[2])
  
  # Run window through data
  ans <- rcpp_seqr2(wind, window)
  colnames(ans$sequence) <- c("start", "stop", "R2")
  sequence <- as.data.frame(ans$sequence)
  y = ans$y
  
  position <- order(sequence$R2, decreasing=T)
  
  #cat(best)
  wind_forecast <- c()
  
  best_index <- 1
  while(best_index <= best && best_index <= nrow(sequence)){
    #for(best_index in 1:best){ # Run it until the highest best match
    pos_train <- sequence$start[position[best_index]] : sequence$stop[position[best_index]]
    model <- lm(y ~ poly(x, poli, raw = T), data=data.frame(x = ts_reg[pos_train]))
    pos_test <- (sequence$stop[position[best_index]]+1):(sequence$stop[position[best_index]]+forecast)
    wind_forecast <- rbind(wind_forecast, predict(model, newdata=data.frame(x = ts[pos_test])))
    best_index <- best_index + 1
  }
  
  if (aggregation == 'median'){
    agg_f = median }
  else {
    agg_f = mean }
  
  return(list("windows" = sequence, 
              "forecasts" = wind_forecast,
              "forecast" = apply(wind_forecast, 2, agg_f)))
}

gsDTScan <- function(input, fh){
  # Define the grid range for each parameter
  poli_grid <- c(1,2,3) #
  best_grid <- c(3,5,7,10,15,25,50)
  window_grid <- fh*c(1, 1.25, 1.5, 2)
  agg_grid <- c('median','mean')
  
  
  # Set first parameter as the best
  poli_sel <- poli_grid[1]
  best_sel <- best_grid[1]
  window_sel <- window_grid[1]
  agg_sel <- agg_grid[1]
  
  # split train/validation sets
  train <- head(input,length(input)-fh)
  test <- tail(input, fh)
  
  err <- Inf
  
  for (poli in poli_grid){
    for (best in best_grid){
      for (window in window_grid){
        for (agg in agg_grid){
          pred <- dtsf(ts=train, poli=poli, best=best, window=window, aggregation=agg, forecast=fh)$forecast
          err_new <- smape_cal(test, pred)
          if (mean(err_new) < mean(err)){
            err <- err_new
            poli_sel <- poli
            best_sel <- best
            window_sel <- window
            agg_sel <- agg
          }
        }
      }
    }
    
  }
  cat('\npoli: ',poli_sel)
  cat('\nbest: ',best_sel)
  cat('\nagg: ',agg_sel)
  cat('\nwindow: ',window_sel)
  return(dtsf(ts=input, poli=poli_sel, best=best_sel, window=window_sel, aggregation=agg_sel, forecast=fh)$forecast)
}



Benchmarks <- function(input, fh){
  #Used to estimate the statistical benchmarks of the M4 competition
  #Estimate seasonaly adjusted time series
  ppy <- frequency(input) ; ST <- F
  #input <- ts(input[!is.na(input)],frequency=ppy)
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==T){
    Dec <- decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    des_input <- input ; SIout <- rep(1, fh)
  }
  
  start_t <- Sys.time()
  f1 <- naive(input, h=fh)$mean #Naive
  t1 <<- t1 + Sys.time() - start_t; start_t <- Sys.time()
  
  f2 <- naive_seasonal(input, fh=fh) #Seasonal Naive
  t2 <<- t2 + Sys.time() - start_t; start_t <- Sys.time()
  
  f3 <- naive(des_input, h=fh)$mean*SIout #Naive2
  t3 <<- t3 + Sys.time() - start_t; start_t <- Sys.time()
  
  f4 <- ses(des_input, h=fh)$mean*SIout #Ses
  t4 <<- t4 + Sys.time() - start_t; start_t <- Sys.time()
  
  f5 <- holt(des_input, h=fh, damped=F)$mean*SIout #Holt
  t5 <<- t5 + Sys.time() - start_t; start_t <- Sys.time()
  
  f6 <- holt(des_input, h=fh, damped=T)$mean*SIout #Dampedwarn
  t6 <<- t6 + Sys.time() - start_t; start_t <- Sys.time()
  
  f7 <- Theta.classic(input=des_input, fh=fh)$mean*SIout #Theta
  t7 <<- t7 + Sys.time() - start_t; start_t <- Sys.time()
  
  f8 <- (f4+f5+f6)/3 #Comb

  f9 <- dtsf(ts=as.numeric(input), poli=1, best=10, window=fh, forecast=fh)$forecast
  t9 <<- t9 + Sys.time() - start_t; start_t <- Sys.time()
  
#  f10 <- gsDTScan(input, fh)
#  t10 = t10 + Sys.time() - start_t; start_t <- Sys.time()
  
  
  return(list(f1,f2,f3,f4,f5,f6,f7,f8,f9)) #,f10
}

#Methods, Horizon, time-series

# 
# ######################
# #setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# ######################
# # MUltithreading forecast
# forecasts <- foreach(i=1:length(data_train), .combine=cbind, 
#                      .packages=c('forecast','DTScanF')) %dopar% {
#   insample <- data_train[[i]]
#   forecast = Benchmarks(input=insample, fh=fh)
#   
#   forecast
# }
# stopCluster(cl)
# 
# 
# 
# # Calculate metrics
# for (i in 1:length(data_train)){
#   insample <- data_train[[i]]
#   outsample <- data_test[[i]]
# 
#   #sMAPE
#   for (j in 1:length(Names_benchmarks)){
#     Total_smape[j,,i] <- smape_cal(outsample, forecasts[[j]]) #j the # of the benchmark
#   }
#   #MASE
#  for (j in 1:length(Names_benchmarks)){
#     Total_mase[j,,i] <- mase_cal(insample, outsample, forecasts[[j]]) #j the # of the benchmark
#   }
# }

datasets <- c('Hourly','Daily','Weekly','Monthly','Quarterly','Yearly')
horizons <- c(48,14,13,18,8,6)
frequencies <- c(24,1,1,12,4,1)

t1 <- t2 <- t3 <- t4 <- t5 <- t6 <- t7 <- t9 <- t10 <- 0

for (dn in seq(1,6)){

  # Get dataset, horizon and frequency
  dataset <- datasets[dn]
  fh <- horizons[dn]
  frq <- frequencies[dn]
  
  print(dataset)

  # Read data
  data_train <- read_data(paste('data/m4/Train/',dataset,'-train.csv',sep=''))
  data_train <- lapply(data_train, function(x) ts(x[!is.na(x)], frequency=frq)) # Transform into TimeSeries object with frequency, drops na
  data_test <- read_data(paste('data/m4/Test/',dataset,'-test.csv',sep=''))
  
  
  # Prepare dataset
  Names_benchmarks <- c("Naive", "sNaive", "Naive2", "SES", "Holt", "Damped", "Theta", "Comb", "DTScanF") #, 'DTScanF2'
  Total_smape=Total_mase <- array(NA,dim = c(length(Names_benchmarks), fh, length(data_train)))
  
  pb <- txtProgressBar(min=0, max=length(data_train), style=3)
  for (i in 1:length(data_train)){
   setTxtProgressBar(pb, i)
  
   insample <- data_train[[i]]
   outsample <- data_test[[i]]
   forecasts <- Benchmarks(input=insample, fh=fh)
  
   #sMAPE
   for (j in 1:length(Names_benchmarks)){
     Total_smape[j,,i] <- smape_cal(outsample, forecasts[[j]]) #j the # of the benchmark
   }
   #MASE
  for (j in 1:length(Names_benchmarks)){
     Total_mase[j,,i] <- mase_cal(insample, outsample, forecasts[[j]]) #j the # of the benchmark
   }
  }
  close(pb)
  
  print(dataset)
  print("########### sMAPE ###############")
  for (i in 1:length(Names_benchmarks)){
    print(paste(Names_benchmarks[i], round(mean(Total_smape[i,,]), 3)))
  }
  print("########### MASE ################")
  for (i in 1:length(Names_benchmarks)){
    print(paste(Names_benchmarks[i], round(mean(Total_mase[i,,]), 3)))
  }
  print("########### OWA ################")
  for (i in 1:length(Names_benchmarks)){
    print(paste(Names_benchmarks[i],
                round(((mean(Total_mase[i,,])/mean(Total_mase[3,,]))+(mean(Total_smape[i,,])/mean(Total_smape[3,,])))/2, 3)))
  }
  
  
  
  # Save results
  saveRDS(Total_smape,paste('rds/',dataset,'-smape.rds',sep=''))
  saveRDS(Total_mase,paste('rds/',dataset,'-mase.rds',sep=''))
}


data_train <- read_data('data/m4/Train/Hourly-train.csv')
data_train <- lapply(data_train, function(x) ts(x[!is.na(x)], frequency=frq)) # Transform into TimeSeries object with frequency, drops na
data_test <- read_data('data/m4/Test/Hourly-test.csv')

sel <- 15
insample <- as.numeric(data_train[[sel]])
outsample <- as.numeric(data_test[[sel]])
mdl <- dtsf(ts=insample, poli=1, best=3, window=48, forecast=48)#$forecast



position <- order(mdl$windows$R2, decreasing=T)

setEPS(width=10, height=6)
postscript("img/example.eps")
eps('img/example.eps', width=10, height=6)

plot(append(insample,outsample), type='l', lwd=1, col='black', xlim=c(350,733), xlab="t", ylab="y",xaxt='n',yaxt='n')
#lines(append(rep(NA, length(insample)-fh),tail(insample,fh)), col='white', type='l', lty=3)
for (i in 1:3){
  st <- position[i]
  wd = window(insample, st, st+fh-1)
  lines(append(rep(NA, st-1),wd), col='black', lwd=2)
  text(x=st+24, y=8900, paste('Analogue', i))
}
abline(v=length(insample), lty=3)
abline(v=(length(insample)-fh), lty=3)
lines(append(rep(NA,length(insample)),mdl$forecast), type='l', col='blue', lwd=2, lty=2)
text(x=724,y=8900, 'Forecast')
text(x=677,y=8900, 'Query')

dev.off()