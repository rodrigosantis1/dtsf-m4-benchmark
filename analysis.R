setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

read_data <- function(path, frequency=F){
  df <- read.csv(path, header=T, row.names='V1')
  df <- t(df)
  df <- lapply(seq_len(ncol(df)), function(i) df[,i])
  return(df)
}



dataset <- 'Hourly'
frq <-1
data_train <- read_data(paste('data/m4/Train/',dataset,'-train.csv',sep=''))
data_train <- lapply(data_train, function(x) ts(x[!is.na(x)], frequency=frq)) # Transform into TimeSeries object with frequency, drops na


#### Hourly analysis

Names_benchmarks <- c("Naive", "sNaive", "Naive2", "SES", "Holt", "Damped", "Theta", "Comb", "DTSF")

Total_smape <- readRDS(paste('rds/',dataset,'-smape.rds',sep=''))
Total_mase  <- readRDS(paste('rds/',dataset,'-mase.rds',sep=''))


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

exclude <- c(-1,-2,-4)
Names_benchmarks <- Names_benchmarks[exclude]
Total_smape <- Total_smape[exclude,,]


#apply(Total_smape, c(1,2), mean)
err_h <- apply(Total_smape, c(2,1), mean)
#pdf('img/boxplot-step.pdf', width=8, height=6)
#setEPS(width=8, height=6)
#postscript("img/boxplot-step.eps")
boxplot(err_h, names=Names_benchmarks, ylab='sMAPE',cex.lab=1.25,cex.axis=1.25) #, main='Step error distribution by method'
#dev.off()

## Plot error through period
nn <- ncol(err_h)
#pdf('img/matplot-step.pdf', width=8, height=6)
setEPS(width=10, height=7)
postscript("img/matplot-step.eps")
matplot(err_h, type='l',col=seq_len(nn),lty=seq_len(nn), lwd=2,ylab='sMAPE', xlab='Step', cex.lab=1.25,cex.axis=1.25) #, main='Average error by step'
legend('topleft', Names_benchmarks,col=seq_len(nn),cex=1,lty=seq_len(nn))
dev.off()

## Error by total_smape
err_ts <- apply(Total_smape, c(3,1), mean)
err_ts <- err_ts[order(err_ts[,6]),]
#pdf('img/boxplot-ts.pdf', width=8, height=6)
#setEPS(width=8, height=6)
#postscript("img/boxplot-ts.eps")
boxplot(err_ts, names=Names_benchmarks, ylab='sMAPE',cex.lab=1.25,cex.axis=1.25) #main='Time series error distribution by method'
#dev.off()


#pdf('img/matplot-ts.pdf', width=8, height=6)
setEPS(width=10, height=7)
postscript("img/matplot-ts.eps")
matplot(log10(err_ts), type='p',col=seq_len(nn),lty=seq_len(nn), lwd=2, ylab='log10(sMAPE)', xlab='Time Series',cex.lab=1.25,cex.axis=1.25) #, main='Average error by time series'
legend('bottomright', Names_benchmarks,col=seq_len(nn),cex=1.0,lty=seq_len(nn))
dev.off()

