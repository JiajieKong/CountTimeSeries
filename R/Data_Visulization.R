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
library(truncnorm)
library(matrixStats)
library(numDeriv)
library(tscount)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(gridExtra)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# try(setwd(paste(getwd(),"/Dropbox/Robert", sep = "")), silent = TRUE)
# try(setwd(paste(getwd())), silent = TRUE)
# setwd("/Users/jiajiekong/dropbox/Robert/" ) #Mac
# setwd("/home/jiajie/Dropbox/Robert/" ) #Linux
# setwd("C:/Users/JiajieKong/Dropbox/Robert") #Windows

# X = dataa #from another file
source("Rainy_Day_Seattle.R")

# campyfit <- tsglm(ts=campy, model=list(past_obs=1:3, past_mean=1:3))
# pit(campyfit)

# options(scipen=999)  # turn off scientific notation like 1e+06
# data("midwest", package = "ggplot2")  # load the data
# # midwest <- read.csv("http://goo.gl/G1K41K") # alt source 
# 
# # Init Ggplot
# gg <- ggplot(midwest, aes(x=area, y=poptotal)) + 
#   geom_point(aes(col=state), size=1) +  # Set color to vary based on state categories.
#   geom_smooth(method="lm", col="firebrick", size=2) + 
#   coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 1000000)) + 
#   labs(title="Area Vs Population", subtitle="From midwest dataset", y="Population", x="Area", caption="Midwest Demographics") + 
#   scale_color_brewer(palette = "Spectral") +
#   scale_x_continuous(breaks=seq(0, 0.1, 0.01)) + 
#   theme_bw()
# plot(gg)


##Full Data
X_org = data.frame(1:1040, X)
colnames(X_org) = c("ind", "X")

ggplot(X_org, aes(x = ind, y = X)) + 
  geom_point(col = "blue", size = 1) + 
  geom_line(linetype = "dashed", size = .2) + 
  labs(title="Weekly Rainy Days in Seatle", y="# of Rainy Days", x="") 




##4 year data
X_new = data.frame(1:208, X[833:1040])
colnames(X_new) = c("ind", "X")

p1 = ggplot(X_new, aes(x = ind, y = X)) + 
  geom_point(col = "blue", size = 1) + 
  geom_line(linetype = "dashed", size = .2) + 
  labs(title="Weekly Rainy Days in Seatle", y="# of Rainy Days", x="Date") +
  scale_x_continuous(breaks=seq(0, 208, length.out = 9), 
                     labels = c("01-2016", "27-2016", "01-2017", "27-2017", "01-2018", "27-2018", "01-2019", "27-2019", "52-2019")) 
  



##acf
ts.acf <- acf(X, plot=FALSE)

alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.acf$n.used)

p2 = ts.acf$acf %>% 
  as_tibble() %>% mutate(lags = 0:(n()-1)) %>% 
  ggplot(aes(x=lags, y = V1)) + scale_x_continuous(breaks=seq(0,41,4)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  labs(y="Sample ACF", x="Lag", title= "") +
  geom_segment(aes(xend=lags, yend=0))


##pacf
ts.pacf <- pacf(X, main=NULL,plot=FALSE)

alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.pacf$n.used)

p3 = ts.pacf$acf %>% 
  as_tibble() %>% mutate(lags = 1:n()) %>%
  ggplot(aes(x=lags, y = V1)) + 
  geom_segment(aes(xend=lags, yend=0))  + 
  scale_x_continuous(breaks=seq(0,41,4))+ 
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  labs(y="Sample PACF", x="Lag", title= "")+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), size = .5)


##Seasonal Analysis
#wrap data#
X_W = matrix(rep(0, 1040), ncol = 52)
for (j in 1:52){X_W[,j] = X[seq(j, 1040, 52)]}

X_box = data.frame(rep(1:52, 20), X)
colnames(X_box) = c("week", "X")

ggplot(X_box, aes(x=as.factor(week), y=X)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  scale_x_discrete(breaks=seq(1, 52, 4)) +
  coord_cartesian(xlim = c(-1, 54)) +
  labs(y = "", x = "Week", title = "Box Plot of Weekly Rainy Days")
                   
sd_df = data.frame(1:52, apply(X_W, MARGIN = 2, var), apply(X_W, MARGIN = 2, mean))
colnames(sd_df) = c("week", "X", "m")
p4 = ggplot(sd_df, aes(x=week, y=X)) + 
  geom_point(col = "blue") +
  geom_line(linetype = "dashed", col = "blue") +
  geom_point(aes(x=week, y=m), col = "red") +
  geom_line(aes(x=week, y=m), linetype = "dashed", col = "red") +
  labs(y = "Red: Mean; Blue: Var")

##Overdisperse##
over_disperse_Bin = function(s, X){
  
  ind = seq(s, length(X), 52)
  pp = sum(X[ind]) / (7*length(X[ind]))
  
  expected = 7 * pp * (1 - pp)
  obs = var(X[ind])
  
  c(obs, expected )
}

over_disperse_Poi = function(s, X){
  
  ind = seq(s, length(X), 52)
  
  expected = mean(X[ind])
  obs = var(X[ind])
  
  c(obs, expected )
}

over_disperse_plot = function(over_disperse){
  
  M = matrix(rep(0, 2 * 52), ncol = 2)
  
  for (j in 1:52){
    
    M[j,] = over_disperse(j, X)
    
  }
  
  M_df = data.frame(M)
  colnames(M_df) = c("obs", "expected")
  
  ggplot(M_df, aes(x=expected, y=obs)) +
    geom_point(col = "blue") +
    geom_abline(intercept = 0, slope = 1, color="red", 
                linetype="dashed", size=1) +
    labs(x = "Expected Variance", y = "Obs Variane")
  
}

over_disperse_plot(over_disperse_Bin)
over_disperse_plot(over_disperse_Poi)

#pdf("Plot_Data_Analysis.pdf", width=10, height=10, paper='special')
grid.arrange(p1, p4, arrangeGrob(p2, p3, ncol = 2), nrow = 3)
#dev.off()






