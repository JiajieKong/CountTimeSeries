rm(list=ls())
library(MASS)
library(mvtnorm)
library(data.table)
library(ellipse)
library(parallel)
library(matrixcalc)
library(psych)
library(MCMCpack)
library(simex)
library(Matrix)
library(SpatioTemporal)


setwd("/Users/jiajiekong/dropbox/Robert/" ) #Mac
# setwd("/home/jiajie/Dropbox/Research/Final_works/" ) #Linux

mydata.1 <- read.csv("storm_data_search_results_part1.csv", header=TRUE, sep=",")
mydata.2 <- read.csv("storm_data_search_results_part2.csv", header=TRUE, sep=",")

mydata = rbind(mydata.1, mydata.2)

time_series = paste(substr(mydata[,4], 7, 10), substr(mydata[,4], 1, 2), sep = "")

data_format = function(ts){
  
  ts_table = as.data.frame(table(ts))
  ts_table[,1] = as.numeric(as.character(ts_table[,1]))
  
  data_range = 201003:201012
  
  for (i in 1:9){data_range = c(data_range, (201001+i*100):(201012+i*100))}
  
  data_range = c(data_range, 202001:202003)
  
  M_info = cbind(1:length(data_range), data_range, 
                 matrix(rep(0, length(data_range)), ncol = 1))
  
  for (i in 1:dim(ts_table)[1]){
    
    aa = data_range == ts_table[i, 1]
    
    M_info[aa,3] = ts_table[i, 2]
    
  }
  
  
  
  
  colnames(M_info) = c("ID", "data", "freq")
  
  return(M_info)
}

t_data = data_format(time_series)

# save a numeric vector containing 121 monthly observations
# from Mar 2010 to Mar 2020 as a time series object
myts <- ts(t_data[,3], start=c(2010, 3), end=c(2020, 3), frequency=12)
# Seasonal decomposition
fit <- stl(myts, s.window="period")

s_para = matrix(rep(0, 2*12), ncol = 2)
m_data = list()
t_data = cbind(t_data, rep(0, dim(t_data)[1]))

for (m in 1:12){
  
  kk = t_data[seq(m,dim(t_data)[1], 12), 3]
  t_data[seq(m,dim(t_data)[1], 12), 4] = (kk - mean(kk)) / sd(kk)
  s_para[m, 1] = mean(kk)
  s_para[m, 2] = sd(kk)
  m_data[[m]] = kk
  
}

t_data[is.nan(t_data[,4]), 4] = 0  #standared 



###standardrize ##
myts_2 <- ts(t_data[,4], start=c(2010, 3), end=c(2020, 3), frequency=12)
par(mfrow = c(1, 1))


###Find correlation of time serise###
pdf(paste("plot_hist", '.pdf', sep = ""), width=10,height=13, paper='special')
par(mfrow = c(4,3 ))
mon = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", 
        "Sep", "Oct", "Nov", "Dec", "Jan", "Feb")

for (m in 1:12){
  hist(m_data[[m]], breaks = 20, main = mon[m], 
       xlab = paste("mean:", round(s_para[m,1],2), ",  sd:", round(s_para[m,2]),sep = ""))
}
dev.off()

pdf(paste("plot_Standard_TS", '.pdf', sep = ""), width=12,height=8, paper='special')
plot(myts_2, ylab = "Seasonal standardized data")
dev.off()

pdf(paste("plot_ACF", '.pdf', sep = ""), width=12,height=8, paper='special')
acf(t_data[,4])
dev.off()





pdf(paste("plot_data", '.pdf', sep = ""), width=27,height=8, paper='special')
plot(myts, ylab = "Tornado Count")
dev.off()

pdf(paste("plot_Seasonal_Analysis_1", '.pdf', sep = ""), width=8,height=8, paper='special')
plot(fit, main = "Seasonal Analysis")
dev.off()

pdf(paste("plot_Seasonal_Analysis_2", '.pdf', sep = ""), width=8,height=8, paper='special')
monthplot(myts, ylab = "Tornado Count", main = "Seasonal Analysis")
dev.off()


# additional plots
library(forecast)
seasonplot(myts)


