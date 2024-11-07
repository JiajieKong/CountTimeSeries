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
library(lubridate)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

mydata <- read.csv("data_2.csv", header=TRUE, sep=",")

tt = rep(0, dim(mydata)[1])

for (j in 1:dim(mydata)[1]){
  
  if (!is.na(mydata[j,4])){
    if (!(mydata[j, 4] == 0)){
      tt[j] = 1
    }
  }
}

mydata = cbind(mydata, tt)


week_sum = function(dataa){
  
  n = length(dataa)
  m = floor(n/7)
  
  s = rep(0, m)
  
  for (j in 1:m){
    
    s[j] = sum(dataa[(j*7-6):(j*7)])
    
  }
  
  return(s)
}


j = 1; i = 1


for (years in 2000:2019){
  
  sel = as.numeric(substr(mydata[,3], 1, 4)) == years
  
  z = week_sum(mydata[sel,7])
  
  beg_week = 53 - length(z)
  
  label = paste(years,"-",beg_week:52, sep = "")
  
  j = c(j, z)
  
  i = c(i, label)
  
}


result <- data.frame("count" = j[-1], "data" = i[-1], "normalized" = rep(0, 1040))

##normalize##

for (j in 1:52){
  
  select_item = result[seq(j, 1040, 52), 1]
  
  normallized_result = (select_item - mean(select_item)) / sd(select_item)
  
  result[seq(j, 1040, 52), 3] = normallized_result
  
}

X <- result[,1]

# ##histogram##
# pdf(paste("hist_weekly_data", '.pdf', sep = ""), width=16,height=30, paper='special')
# par(mfrow = c(9, 6))
# for (j in 1:52){
#   
#   hist(result[seq(j, 1040, 52),1], 
#        xlab = paste("week:", j),
#        xlim = c(0, 7))
#   
# }
# dev.off()
# 
# k = rep(0, 52)
# for (j in 1:52){
#   k[j] = mean(result[seq(j, 1040, 52), 1])
# }
# 
# 
# pdf(paste("weekly_mean_plot", '.pdf', sep = ""), width=8,height=6, paper='special')
# plot(ts(k))
# dev.off()
# 
# 
# 
# 
# 
# tes = ts(result[521:1040,1], freq=365.25/7, start=decimal_date(ymd("2000-01-01")))
# plot(tes)
# 
# acf(result[1:1040,3])
# 
# 
# ##save plot##
# pdf(paste("plot_new_weekly_data", '.pdf', sep = ""), width=12,height=8, paper='special')
# tes = ts(result[521:1040,1], freq=365.25/7, start=decimal_date(ymd("2010-01-01")))
# plot(tes, ylab = "weekly count")
# dev.off()
# 
# pdf(paste("plot_new_standardize_weekly_data", '.pdf', sep = ""), width=12,height=8, paper='special')
# tes = ts(result[521:1040,3], freq=365.25/7, start=decimal_date(ymd("2010-01-01")))
# plot(tes, ylab = "standardized weekly count")
# dev.off()
# 
# pdf(paste("plot_acf_standardize_weekly_data", '.pdf', sep = ""), width=12,height=8, paper='special')
# acf(result[521:1040,3], ylab = "autocorrelation")
# dev.off()
# 
# 
# 
# 
# time_series = paste(substr(mydata[,4], 7, 10), substr(mydata[,4], 1, 2), sep = "")
# 
# data_format = function(ts){
#   
#   ts_table = as.data.frame(table(ts))
#   ts_table[,1] = as.numeric(as.character(ts_table[,1]))
#   
#   data_range = 201003:201012
#   
#   for (i in 1:9){data_range = c(data_range, (201001+i*100):(201012+i*100))}
#   
#   data_range = c(data_range, 202001:202003)
#   
#   M_info = cbind(1:length(data_range), data_range, 
#                  matrix(rep(0, length(data_range)), ncol = 1))
#   
#   for (i in 1:dim(ts_table)[1]){
#     
#     aa = data_range == ts_table[i, 1]
#     
#     M_info[aa,3] = ts_table[i, 2]
#     
#   }
#   
#   
#   
#   
#   colnames(M_info) = c("ID", "data", "freq")
#   
#   return(M_info)
# }
# 
# t_data = data_format(time_series)
# 
# # save a numeric vector containing 121 monthly observations
# # from Mar 2010 to Mar 2020 as a time series object
# myts <- ts(t_data[,3], start=c(2010, 3), end=c(2020, 3), frequency=12)
# # Seasonal decomposition
# fit <- stl(myts, s.window="period")
# 
# s_para = matrix(rep(0, 2*12), ncol = 2)
# m_data = list()
# t_data = cbind(t_data, rep(0, dim(t_data)[1]))
# 
# for (m in 1:12){
#   
#   kk = t_data[seq(m,dim(t_data)[1], 12), 3]
#   t_data[seq(m,dim(t_data)[1], 12), 4] = (kk - mean(kk)) / sd(kk)
#   s_para[m, 1] = mean(kk)
#   s_para[m, 2] = sd(kk)
#   m_data[[m]] = kk
#   
# }
# 
# t_data[is.nan(t_data[,4]), 4] = 0  #standared 
# 
# 
# 
# ###standardrize ##
# myts_2 <- ts(t_data[,4], start=c(2010, 3), end=c(2020, 3), frequency=12)
# par(mfrow = c(1, 1))
# 
# 
# ###Find correlation of time serise###
# pdf(paste("plot_hist", '.pdf', sep = ""), width=10,height=13, paper='special')
# par(mfrow = c(4,3 ))
# mon = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", 
#         "Sep", "Oct", "Nov", "Dec", "Jan", "Feb")
# 
# for (m in 1:12){
#   hist(m_data[[m]], breaks = 20, main = mon[m], 
#        xlab = paste("mean:", round(s_para[m,1],2), ",  sd:", round(s_para[m,2]),sep = ""))
# }
# dev.off()
# 
# pdf(paste("plot_Standard_TS", '.pdf', sep = ""), width=12,height=8, paper='special')
# plot(myts_2, ylab = "Seasonal standardized data")
# dev.off()
# 
# pdf(paste("plot_ACF", '.pdf', sep = ""), width=12,height=8, paper='special')
# acf(t_data[,4])
# dev.off()
# 
# 
# 
# 
# 
# pdf(paste("plot_data", '.pdf', sep = ""), width=27,height=8, paper='special')
# plot(myts, ylab = "Tornado Count")
# dev.off()
# 
# pdf(paste("plot_Seasonal_Analysis_1", '.pdf', sep = ""), width=8,height=8, paper='special')
# plot(fit, main = "Seasonal Analysis")
# dev.off()
# 
# pdf(paste("plot_Seasonal_Analysis_2", '.pdf', sep = ""), width=8,height=8, paper='special')
# monthplot(myts, ylab = "Tornado Count", main = "Seasonal Analysis")
# dev.off()
# 
# 
# # additional plots
# seasonplot(myts)


