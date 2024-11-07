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

# try(setwd(paste(getwd(),"/Dropbox/Robert/", sep = "")), silent = TRUE)
# try(setwd(paste(getwd())), silent = TRUE)
# setwd("/Users/jiajiekong/dropbox/Robert/" ) #Mac
# setwd("/home/jiajie/Dropbox/Robert/" ) #Linux
# setwd("C:/Users/JiajieKong/Dropbox/Robert/Wiley") #Windows

# X = dataa #from another file
source("Rainy_Day_Seattle.R")

Cov_function_SAR = function(i, j, TT = TT, a = a, phi = phi){
  
  h = abs(i - j)
  (a^h + phi*a^(TT - h))/(1 + phi*a^TT)
  
}

two_stat_MC = function(a, b){
  
  p_0_l = matrix(rep(0, 10*12), nrow = 10)
  p_1_l = matrix(rep(0, 10*12), nrow = 10)
  
  p_0_l[1,2] = a
  p_0_l[1,3] = 1 - a
  p_1_l[1,2] = 1 - b 
  p_1_l[1,3] = b
  
  for (l in 2:10){
    for (k in 0:l){
      
      p_0_l[l, k + 2] = a * p_0_l[l - 1, k + 2] + (1 - a) * p_1_l[l - 1, k + 1]
      p_1_l[l, k + 2] = b * p_1_l[l - 1, k + 1] + (1 - b) * p_0_l[l - 1, k + 2]
    }
  }
  
  return(list(p_0 = p_0_l, p_1 = p_1_l))
}


obj_function = function(par_star, X, model_select, Rp, org, TT){ 
  
  model_output = model_select(par_star, org, TT)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  Test_Pass = model_output$test_pass
  lag = model_output$lag
  
  if (Test_Pass == FALSE){
    
    FF = -1e22
    par = round(par, 3)
    print(paste(FF, par[1], par[2], par[3], par[4]))
    return(FF)
    
  }
  
  if (model_output$S == FALSE){
    
    N = dim(Rp)[1]
    Z = matrix(rep(0, N*length(X)), nrow = N)
    w = matrix(rep(0, N*length(X)), nrow = N)
    epsilon = matrix(rep(0, N*length(X)), nrow = N)
    log_w = matrix(rep(0, N*length(X)), nrow = N)
    
    for (t in 1:lag){
      P_C_1 = p(t)$PC_1
      phi_1 = phi(t)
      Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
      P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
      if (t == 1){
        log_w[,t] = log(P_x1) 
      }else{
        log_w[,t] = log_w[,t-1] + log(P_x1)
      }
    }
    
    for (t in (lag+1):length(X)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[,t] = log_w[,t - 1] + log(w_Z)
      
    }
    log_wt = log_w[,t]
    
  }else{
    
    phi = model_output$phi
    a = model_output$a
    
    N = dim(Rp)[1]
    Z = matrix(rep(0, N*length(X)), nrow = N)
    w = matrix(rep(0, N*length(X)), nrow = N)
    epsilon = matrix(rep(0, N*length(X)), nrow = N)
    log_w = matrix(rep(0, N*length(X)), nrow = N)
    
    P_C_1 = p(1)$PC_1
    phi_1 = phi
    Z[,1] = qtruncnorm(Rp[,1], a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[,1] = log(P_x1)
    
    
    for (t in 2:(lag+1)){
      
      CC = outer(1:t, 1:t, FUN = Cov_function_SAR, TT = TT, a = a, phi = phi)
      
      L = CC[1:(t-1),1:(t-1)]
      L_L = chol(L)
      
      D = forwardsolve(t(L_L), CC[1:(t-1),t])
      D_mu = apply(t(backsolve(L_L, D) * t(Z[,1:t-1])), MARGIN = 1, sum)
      
      
      D_sig = CC[t,t] - sum(D^2)
      
      P_C_1 = p(t %% TT)$PC_1
      phi_1 = phi
      Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2],
                         mean = D_mu, sd = sqrt(D_sig))
      P_x1 = pnorm(P_C_1[X[t] + 2], mean = D_mu, sd = sqrt(D_sig)) -
        pnorm(P_C_1[X[t] + 1], mean = D_mu, sd = sqrt(D_sig))
      
      log_w[, t] = log_w[, t-1] + log(P_x1)
      
    }
    
    
    for (t in (lag+2):length(X)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi
      sigma_v = sqrt((1-a^2)*(1-phi^2)*(1-phi*a^TT)/(1+phi*a^TT))
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v + a*Z[,t-1] - a*phi_v*Z[,t-lag-1]
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[,t] = log_w[,t - 1] + log(w_Z)
      
    }
    
    log_wt = log_w[,t]
    
  }
  
  par = round(par, 3)
  
  FF = - log(N) + logSumExp(log_wt)
  
  print(paste(FF, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]))
  
  FF
}


main = function(X, model_select, para_int_select, N, TT){
  
  Rp = matrix(runif(length(X) * N), nrow = N)
  para_int_select_tran = model_select(para_int_select, org = TRUE, TT)$par_tran
  
  system.time({
    MLE_1 = optim(para_int_select_tran, obj_function, 
                  X = X, model_select = model_select, Rp = Rp, org = FALSE, TT = TT,
                  method = "BFGS", control = list(factr = 1e7, fnscale=-1))
    
    MLE_1$par_est = model_select(MLE_1$par, org = FALSE, TT)$par
    
    MLE_1$hessian_2 = hessian(func = obj_function, MLE_1$par_est, org = TRUE,
                              X = X, model_select = model_select, Rp = Rp, TT = TT,
                              method.args=list(eps = 1e-4, d=0.01))
    # MLE_1$se = sqrt(diag(solve(-MLE_1$hessian_2)))
    
    MLE_1$AIC = -2*MLE_1$value + 2 * length(para_int_select_tran)
    MLE_1$BIC = -2*MLE_1$value + log(length(X)) * length(para_int_select_tran)
  })
  
  MLE_1
}


##simulation_begin_##
simulation_model = function(model_select, par, NN, TT, L){ #NN = 200
  
  model_output = model_select(par, org = TRUE, TT)
  p = model_output$p
  phi = model_output$phi
  lag = model_output$lag
  
  X_simulate = list()
  Z_simulate = list()
  
  if (model_output$S == FALSE){
    
    for (i in 1:NN){
      
      Z_simulate[[i]] = rep(0, L)
      X_simulate[[i]] = rep(0, L)
      
      for (j in 1:lag){
        Z_simulate[[i]][j] = rnorm(1)
        X_simulate[[i]][j] = sum(pnorm(Z_simulate[[i]][j]) > p(j)$PC_2 ) - 1
      }
      
      for (j in (lag+1):L){
        
        v = j %% TT
        phi_v = phi(v)
        sigma_v = sqrt(1 - phi_v^2)
        
        Z_simulate[[i]][j] = phi_v * Z_simulate[[i]][j-lag] + sigma_v * rnorm(1)
        X_simulate[[i]][j] = sum(pnorm(Z_simulate[[i]][j]) > p(v)$PC_2 ) - 1
      }
      
    }
  }else{
    
    for (i in 1:NN){
      
      lag = TT
      a = model_output$a
      phi = model_output$phi
      
      Z_simulate[[i]] = rep(0, L)
      X_simulate[[i]] = rep(0, L)
      
      CC = outer(1:(lag+1), 1:(lag+1), FUN = Cov_function_SAR, TT = TT, a = a, phi = phi)
      
      Z_simulate[[i]][1:(lag+1)] = rmvnorm(1, mean = rep(0, lag+1), sigma = CC)
      
      for (j in 1:(lag+1)){
        X_simulate[[i]][j] = sum(pnorm(Z_simulate[[i]][j]) > p(j)$PC_2 ) - 1
      }
      
      for (j in (lag+2):L){
        
        v = j %% TT
        phi_v = phi
        
        sigma_v = sqrt((1-a^2)*(1-phi^2)*(1-phi*a^TT)/(1+phi*a^TT))
        
        Z_simulate[[i]][j] = Z_simulate[[i]][j-lag] * phi_v + a*Z_simulate[[i]][j-1] - 
          a*phi_v*Z_simulate[[i]][j-lag-1] + sigma_v * rnorm(1)
        X_simulate[[i]][j] = sum(pnorm(Z_simulate[[i]][j]) > p(v)$PC_2 ) - 1
      }
      
    }
    
  }
  
  return(X_simulate)
}

simulation_learn = function(model_select, par_int, NN, TT, L, N, simulate_data){
  
  system.time({
    mc <- getOption("mc.cores", 9)
    simulation_MLE <- mclapply(mc.cores = mc, simulate_data, FUN = main, model_select = model_select, 
                               para_int_select = par_int, N = N, TT = TT)
    
  })
  
  return(simulation_MLE)
}

main_2 = function(model_select, par_int, NN, TT, L, N){
  
  simulation_data = simulation_model(model_select, par_int, NN, TT, L)
  simulation_MLE = simulation_learn(model_select, par_int, NN, TT, L, N, simulation_data)
  
  return(list(simulation_data = simulation_data, simulation_MLE = simulation_MLE))
}



simulation_data = list()
simulation_MLE = list()
model_simulate = list()
par_int_ss = list()
model = list()

##Define Model##
model_simulate_possion_Par1 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(log(par[1]),
                 log(par[2]),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi / 2)),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / TT - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[1] = exp(par[1])
    par[2] = exp(par[2])
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4])) / pi * 2
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * TT
    
    if (((par[1] - abs(par[2])) < 0) ||
        ((abs(par[4] + par[5]) > 1) || (abs(par[4] - par[5]) > 1))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] + par[5] * cos((2*pi*(v - par[6]))/ TT)
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1, S = FALSE))
  
} 
model_simulate_NB_Par1 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi / 2)),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / TT - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4])) / pi * 2
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * TT
    
    if (((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1)) ||
        ((abs(par[4] + par[5]) > 1) || (abs(par[4] - par[5]) > 1))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pnbinom(0:100, size = 10, prob = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] + par[5] * cos((2*pi*(v - par[6]))/ TT)
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1, S = FALSE))
  
} 
model_simulate_TM_Par1_10 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / TT - pi/2)),
                 tan((par[7] * pi / 2)),
                 tan((par[8] * pi / 2)),
                 tan((par[9] * pi / TT - pi/2)))
    test_pass = TRUE
  }else{
    par_tran = par = par_star
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * TT
    par[7] = (atan(par_star[7])) / pi * 2
    par[8] = (atan(par_star[8])) / pi * 2
    par[9] = (atan(par_star[9]) + pi/2) / pi * TT
    
    a = par[1]; b = par[2]; c = par[4]; d = par[5]; e = par[7]; f = par[8]
    
    if (((((a + b) < 0) || ((a - b) < 0)) || (((a + b) > 1) || ((a - b) > 1))) ||
        ((((c + d) < 0) || ((c - d) < 0)) || (((c + d) > 1) || ((c - d) > 1))) ||
        ((((e + f) < 0) || ((e - f) < 0)) || (((e + f) > 1) || ((e - f) > 1)))) {
      
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  
  p = function(v){
    
    a = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    b = par[4] + par[5] * cos((2*pi*(v - par[6]))/ TT)
    
    two_stat_MC_result = two_stat_MC(a, b)
    two_stat_MC_den = a/(a + b)*two_stat_MC_result$p_0[10, 2:12] + 
      b/(a + b)*two_stat_MC_result$p_1[10, 2:12]
    
    C = cumsum(two_stat_MC_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:10]), Inf)
    PC_2 = c(0, C[1:10], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[7] + par[8] * cos((2*pi*(v - par[9]))/ TT)
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1, S = FALSE))
} ##TM+PAR(1)

model_simulate_possion_Sar1 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(log(par[1]),
                 log(par[2]),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi / 1 - pi/2)),
                 tan((par[5] * pi / 1 - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[1] = exp(par[1])
    par[2] = exp(par[2])
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5]) + pi/2) / pi * 1
    
    if (((par[1] - abs(par[2])) < 0)){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = ppois(0:100, lambda = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[4]
  a = par[5]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, 
              par_tran=par_tran, lag = TT, S = TRUE, a = a))
  
} 
model_simulate_NB_Sar1 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi / 1 - pi/2)),
                 tan((par[5] * pi / 1 - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5]) + pi/2) / pi * 1
    
    if (((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pnbinom(0:100, size = 10, prob = probability)
    
    PC_1 = c(-Inf, qnorm(C), Inf)
    PC_2 = c(0, C, 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[4]
  a = par[5]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, 
              par_tran=par_tran, lag = TT, S = TRUE, a = a))  
} 
model_simulate_TM_Sar1_10 = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / TT - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / TT - pi/2)),
                 tan((par[7] * pi / 1 - pi/2)),
                 tan((par[8] * pi / 1 - pi/2)))
    test_pass = TRUE
  }else{
    par_tran = par = par_star
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * TT
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * TT
    par[7] = (atan(par_star[7]) + pi/2) / pi * 1
    par[8] = (atan(par_star[8]) + pi/2) / pi * 1
    
    a = par[1]; b = par[2]; c = par[4]; d = par[5]
    
    if (((((a + b) < 0) || ((a - b) < 0)) || (((a + b) > 1) || ((a - b) > 1))) ||
        ((((c + d) < 0) || ((c - d) < 0)) || (((c + d) > 1) || ((c - d) > 1)))){
      
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  
  
  p = function(v){
    
    a = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    b = par[4] + par[5] * cos((2*pi*(v - par[6]))/ TT)
    
    two_stat_MC_result = two_stat_MC(a, b)
    two_stat_MC_den = a/(a + b)*two_stat_MC_result$p_0[10, 2:12] + 
      b/(a + b)*two_stat_MC_result$p_1[10, 2:12]
    
    C = cumsum(two_stat_MC_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:10]), Inf)
    PC_2 = c(0, C[1:10], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = par[7]
  a = par[8]
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, 
              par_tran = par_tran, lag = TT, S = TRUE, a = a))
} ##TM+SAR_2(1)



model[[1]] = model_simulate_possion_Par1
model[[2]] = model_simulate_NB_Par1
model[[3]] = model_simulate_possion_Sar1
model[[4]] = model_simulate_NB_Sar1
model[[5]] = model_simulate_TM_Par1_10
model[[6]] = model_simulate_TM_Sar1_10

par_int_ss[[1]] = c(10, 5, 5, .5, .2, 5)
par_int_ss[[2]] = c(.5, .2, 5, .5, .2, 5)
par_int_ss[[3]] = c(10, 5, 5, .5, .3)
par_int_ss[[4]] = c(.5, .2, 5, .5, .3)
par_int_ss[[5]] = c(.4, .2, 5, .5, .3, 5, .2, .1, 5)
par_int_ss[[6]] = c(.4, .2, 5, .5, .3, 5, .5, .3)

# system.time({
#   
#   kk = list()
#   NN = 5e2  ## number of simulated sample paths
#   N = 5e2   ## number of sample paths for particle filetering approximation
#   Period = c(10, 50) ## length of 1 period
#   Time_Length = c(1e2, 3e2) ## length of the complete timestamp
# 
#   for (j in 1:4){
# 
#     kk[[j]] = list()
# 
#     kk[[j]][[1]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[1], Time_Length[1], N)
#     kk[[j]][[2]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[2], Time_Length[1], N)
#     kk[[j]][[3]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[1], Time_Length[2], N)
#     kk[[j]][[4]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[2], Time_Length[2], N)
# 
#     print(j)
# 
#   }
#   
# })

system.time({
  
  kk = list()
  NN = 500  ## number of simulated sample paths
  N = 500   ## number of sample paths for particle filtering approximation
  Period = c(10, 50) ## length of 1 period
  Time_Length = c(100, 300) #c(1e2, 3e2) ## length of the complete timestamp
  
  for (j in c(5)){
    
    kk[[j]] = list()
    
    #kk[[j]][[1]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[1], Time_Length[1], N)
    kk[[j]][[2]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[2], Time_Length[1], N)
    #kk[[j]][[3]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[1], Time_Length[2], N)
    #kk[[j]][[4]] = main_2(model[[j]], par_int_ss[[j]], NN, Period[2], Time_Length[2], N)
    
    print(j)
    
  }
  
})


#saveRDS(kk, file = "simulation_result_March152022")
kk1 = readRDS("simulation_result_March152022")
kk2 = readRDS("simulation_result_May4")
for (j in c(1,3,5,6)){
  kk2[[j]][[5]] = kk1[[j]][[1]]
  kk2[[j]][[6]] = kk1[[j]][[2]]
}
kk = kk2


result_output = function(s, kk, model_input, par_int_ss, par_num, title, subtitle){
  
  W = length(kk[[s]][[1]]$simulation_MLE[[1]]$par)
  NN = length(kk[[s]][[1]]$simulation_MLE)
  
  aa = list()
  aa2 = list()
  tt = list()
  tt2 = list()
  
  mu_est = list()
  sd_est = list()
  sd_Info_est =list()
  
  mle_output = list()
  model_input = model
  
  for (p in 1:6){
    
    aa[[p]] = matrix(rep(NaN, NN*W), nrow = NN)
    aa2[[p]] = matrix(rep(NaN, NN*W), nrow = NN)
    
    for (j in 1:NN){
      
      aa[[p]][j,] = kk[[s]][[p]]$simulation_MLE[[j]]$par_est
      
      try({aa2[[p]][j,] = 
        sqrt(diag(solve(-kk[[s]][[p]]$simulation_MLE[[j]]$hessian_2)))}, 
        silent = TRUE)
      
    }
    
    new_index = (apply(aa2[[p]], MARGIN = 1, function(x) length(which(!is.na(x))))) == W
    tt[[p]] = aa[[p]][new_index,]
    tt2[[p]] = aa2[[p]][new_index,]
    
    mu_est[[p]] = apply(tt[[p]], MARGIN = 2, mean, na.rm = TRUE)
    sd_est[[p]] = apply(tt[[p]], MARGIN = 2, sd, na.rm = TRUE)
    sd_Info_est[[p]] = apply(tt2[[p]], MARGIN = 2, mean, na.rm = TRUE)
    mle_output[[p]] = kk[[s]][[p]]$simulation_MLE
    
    print(p)
  }
  
  kk_new = as.data.frame(rbind(tt[[1]], tt[[2]],
                               tt[[3]], tt[[4]],
                               tt[[5]], tt[[6]]))
  
  kk_new$scheme = c(rep("A", dim(tt[[1]])[1]),rep("B", dim(tt[[2]])[1]),
                    rep("C", dim(tt[[3]])[1]),rep("D", dim(tt[[4]])[1]),
                    rep("E", dim(tt[[5]])[1]),rep("F", dim(tt[[6]])[1]))
  
  kk_new$size = c(rep(100, dim(tt[[1]])[1]),rep(100, dim(tt[[2]])[1]),
                  rep(300, dim(tt[[3]])[1]),rep(300, dim(tt[[4]])[1]),
                  rep(1000, dim(tt[[5]])[1]),rep(1000, dim(tt[[6]])[1]))
  
  kk_new$period = c(rep("Period=10", dim(tt[[1]])[1]),
                    rep("Period=50", dim(tt[[2]])[1]),
                    rep("Period=10", dim(tt[[3]])[1]),
                    rep("Period=50", dim(tt[[4]])[1]),
                    rep("Period=10", dim(tt[[5]])[1]),
                    rep("Period=50", dim(tt[[6]])[1]))
  
  kk_ff = data.frame(x_full = c(as.matrix(kk_new[,1:W])),
                     period = rep(kk_new$period, W),
                     size = rep(kk_new$size, W),
                     par_num = rep(letters[1:W], 1, each = dim(kk_new)[1]))
  
  levels(kk_ff$par_num) <- par_num[[s]]
  #levels(kk_ff$par_num) <- factor(kk_ff$par_num,levels=unique(kk_ff$par_num))
  
  pd.long$Site <- factor(pd.long$Site,levels=unique(pd.long$Site))
  
  a_mean = data.frame(mean_val = par_int_ss[[s]], par_num = levels(kk_ff$par_num))
  
  pp = ggplot(kk_ff, aes(x = factor(size), y = x_full))  + 
    geom_boxplot(aes(fill = period), position = "dodge", 
                 size = .3, outlier.size = .1) +
    facet_wrap(~par_num, scales = "free_y", nrow = 1,
               labeller = label_parsed)+
    geom_hline(data= a_mean, aes(yintercept=mean_val), linetype = 4, color = "red") +
    scale_fill_manual(values = c("white", "#2cb2aa", "green", "yellow")) +
    labs(title = title[[s]],
         subtitle = subtitle[[s]],
         x = "Length of Sample Path",
         y = "Parameter Estimate") +
    theme(legend.position="bottom") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # box_plot = function(model_input, mle_output, para_int, mu_est, sd_est, sd_Info_est, tt){
  #   
  #   par(mfrow=c(2,3))
  #   
  #   
  #   for (j in 1:dim(tt)[2]){
  #     
  #     boxplot(tt[,j], main = paste("mu_est=",round(mu_est[j],3),"; sd_est=", round(sd_est[j],3), sep = ""), 
  #             xlab = paste("sd_from_InfoMatrix: ", round(sd_Info_est[j],3), sep = ""))
  #     abline(h=para_int[j], col="red", lty = 2)
  #     
  #   }
  #   
  # }
  # 
  # gg_box_plot = function(model_input, mle_output, para_int, mu_est, sd_est, sd_Info_est, tt){
  # 
  # box_plot(model[[s]], kk[[s]][[p]]$simulation_MLE, 
  #          par_int_ss[[s]], mu_est, sd_est, sd_Info_est, tt)
  
  return(list(pp = pp, mu_est = mu_est, sd_est = sd_est, sd_Info_est = sd_Info_est))
}


par_num = list(c( expression(hat(a)[1]),
                  expression(hat(a)[2]),
                  expression(hat(a)[3]),
                  expression(hat(b)[1]),
                  expression(hat(b)[2]),
                  expression(hat(b)[3])),
               c( expression(hat(a)[1]),
                  expression(hat(a)[2]),
                  expression(hat(a)[3]),
                  expression(hat(b)[1]),
                  expression(hat(b)[2]),
                  expression(hat(b)[3])),
               c( expression(hat(a)[1]),
                  expression(hat(a)[2]),
                  expression(hat(a)[3]),
                  expression(hat(phi)),
                  expression(hat(alpha))),
               c( expression(hat(a)[1]),
                  expression(hat(a)[2]),
                  expression(hat(a)[3]),
                  expression(hat(phi)),
                  expression(hat(alpha))),
               c( expression(hat(c)[1]),
                  expression(hat(c)[2]),
                  expression(hat(c)[3]),
                  expression(hat(d)[1]),
                  expression(hat(d)[2]),
                  expression(hat(d)[3]),
                  expression(hat(b)[1]),
                  expression(hat(b)[2]),
                  expression(hat(b)[3])),
               c( expression(hat(c)[1]),
                  expression(hat(c)[2]),
                  expression(hat(c)[3]),
                  expression(hat(d)[1]),
                  expression(hat(d)[2]),
                  expression(hat(d)[3]),
                  expression(hat(phi)),
                  expression(hat(alpha))))

title = list(expression("Model 1: Poisson("*lambda*"(v))+"*plain(PAR(1))),
             expression("Model 3: NB("*p*"(v),1)"*plain(PAR(1))),
             expression("Model 2: Poisson("*lambda*"(v))+"*plain(SAR(1))),
             expression("Model 4: NB("*p*"(v),1)"*plain(SAR(1))),
             expression("Model 3: TSMC("*alpha*","*beta*")+"*plain(PAR(1))),
             expression("Model 4: TSMC("*alpha*","*beta*")+"*plain(SAR(1))))


subtitle = list(expression(paste(lambda*"("*v*")"==a[1]+a[2]*plain(cos)*"("*2*pi*"("*v-a[3]*")"/T*"), "," ",
                                 phi*"("*v*")"==b[1]+b[2]*plain(cos)*"("*2*pi*"("*v-b[3]*")"/T*")")),
                expression(paste(p*"("*v*")"==a[1]+a[2]*plain(cos)*"("*2*pi*"("*v-a[3]*")"/T*"), "," ",
                                 phi*"("*v*")"==b[1]+b[2]*plain(cos)*"("*2*pi*"("*v-b[3]*")"/T*")")),
                expression(paste(lambda*"("*v*")"==a[1]+a[2]*plain(cos)*"("*2*pi*"("*v-a[3]*")"/T*"), "," ",
                                 X[t]==phi*X[t-1]+eta[t]*","," ",
                                 eta[t]==alpha*eta[t-1]+epsilon[t])),
                expression(paste(p*"("*v*")"==a[1]+a[2]*plain(cos)*"("*2*pi*"("*v-a[3]*")"/T*"), "," ",
                                 X[t]==phi*X[t-1]+eta[t]*","," ",
                                 eta[t]==alpha*eta[t-1]+epsilon[t])),
                expression(paste(alpha*"("*v*")"==c[1]+c[2]*plain(cos)*"("*2*pi*"("*v-c[3]*")"/T*"), ",
                                 beta*"("*v*")"==d[1]+d[2]*plain(cos)*"("*2*pi*"("*v-d[3]*")"/T*"), "," ",
                                 phi*"("*v*")"==b[1]+b[2]*plain(cos)*"("*2*pi*"("*v-b[3]*")"/T*")")),
                expression(paste(alpha*"("*v*")"==c[1]+c[2]*plain(cos)*"("*2*pi*"("*v-c[3]*")"/T*"), ",
                                 beta*"("*v*")"==d[1]+d[2]*plain(cos)*"("*2*pi*"("*v-d[3]*")"/T*"), "," ",
                                 X[t]==phi*X[t-1]+eta[t]*","," ",
                                 eta[t]==alpha*eta[t-1]+epsilon[t])))

result_return = list()

result_return[[1]] = result_output(1, kk, model, par_int_ss, par_num, title, subtitle)
#result_return[[2]] = result_output(2, kk, model, par_int_ss, par_num, title, subtitle)
result_return[[3]] = result_output(3, kk, model, par_int_ss, par_num, title, subtitle)
#result_return[[4]] = result_output(4, kk, model, par_int_ss, par_num, title, subtitle)
result_return[[5]] = result_output(5, kk, model, par_int_ss, par_num, title, subtitle)
result_return[[6]] = result_output(6, kk, model, par_int_ss, par_num, title, subtitle)

# pdf("Plot_Simulation1_Analysis.pdf", width=10, height=8, paper='special')
# grid.arrange(result_return[[1]]$pp, 
#              result_return[[3]]$pp,
#              nrow = 2)
# dev.off()
# 
# pdf("Plot_Simulation2_Analysis.pdf", width=10, height=8, paper='special')
# grid.arrange(result_return[[2]]$pp, 
#              result_return[[4]]$pp,
#              nrow = 2)
# dev.off()
# 
# pdf("Plot_Simulation3_Analysis.pdf", width=10, height=8, paper='special')
# grid.arrange(result_return[[5]]$pp, 
#              result_return[[6]]$pp,
#              nrow = 2)
# dev.off()
w = 10
h = 4

pdf("boxPlot_Simulation1_Analysis_try.pdf", width=w, height=h, paper='special')
result_return[[1]]$pp
dev.off()

pdf("boxPlot_Simulation2_Analysis_try.pdf", width=w, height=h, paper='special')
result_return[[3]]$pp
dev.off()

pdf("boxPlot_Simulation3_Analysis_try.pdf", width=w, height=h, paper='special')
result_return[[5]]$pp
dev.off()

pdf("boxPlot_Simulation4_Analysis_try.pdf", width=w, height=h, paper='special')
result_return[[6]]$pp
dev.off()



print_table = function(s, result_return){
  
  for (j in 1:4){
    
    print(round(result_return[[s]]$mu_est[[j]], 5))
    print(round(result_return[[s]]$sd_est[[j]], 5))
    print(round(result_return[[s]]$sd_Info_est[[j]], 5))
    
  }
  
}

print_table(5, result_return)



M = 100000
sam_mcmc = rep(NA, M)
a = -3
v_cur = 0

for (i in 1:M){
  
  v_pro = v_cur + rnorm(1, 0, 1)
  
  if (v_pro < a){
    log_w = -Inf
  }else{
    log_w = dnorm(v_pro, log = TRUE) - dnorm(v_cur, log = TRUE) 
  }
  
  if (log(runif(1)) < log_w){
    v_cur = v_pro
  }
  
  sam_mcmc[i] = v_cur
  
}

hist(sam_mcmc, breaks = 100)

##acf
qqplot.acf <- function(x, name){
  
  ts.acf <- acf(x, plot=FALSE)
  
  alpha <- 0.95
  conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.acf$n.used)
  
  p_acf = ts.acf$acf %>% 
    as_tibble() %>% mutate(lags = 0:(n()-1)) %>% 
    ggplot(aes(x=lags, y = V1)) + 
    scale_x_continuous(breaks=seq(0,41,4)) +
    geom_hline(yintercept=conf.lims, lty=2, col='blue') +
    labs(y="Sample ACF", x="Lag", title= name) +
    geom_segment(aes(xend=lags, yend=0))+
    geom_hline(yintercept = 0)
  
  return(p_acf)
  
}
qqplot.pacf <- function(x){
  
  ts.pacf <- pacf(x, main=NULL,plot=FALSE)
  
  alpha <- 0.95
  conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.pacf$n.used)
  
  p_pacf = ts.pacf$acf %>% 
    as_tibble() %>% mutate(lags = 1:n()) %>%
    ggplot(aes(x=lags, y = V1)) + 
    geom_segment(aes(xend=lags, yend=0))  + 
    scale_x_continuous(breaks=seq(0,41,4))+ 
    geom_hline(yintercept=conf.lims, lty=2, col='blue') +
    labs(y="Sample PACF", x="Lag", title= "") +
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), size = .44)+
    geom_hline(yintercept = 0)
  
  return(p_pacf)
}
##4 year data
qqplot.ts <- function(x){
  X_new = data.frame(1:length(x), x)
  colnames(X_new) = c("ind", "X")
  
  p=ggplot(X_new, aes(x = ind, y = X)) + 
    geom_point(col = "blue", size = .6) + 
    geom_line(linetype = "dashed", size = .2) + 
    labs(title="", y="Value", x="Time Index") 
  
  return(p)
}

NN = 100
TT = 10
L = 100
model_select = model[[1]]
par_int = par_int_ss[[1]]
simulation_data = simulation_model(model_select, par_int, NN, TT, L)
dataa = simulation_data[[1]]
#standardrize
mu_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = mean)
sd_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = sd)
for (i in 1:length(dataa)){
  v = i %% TT
  
  if (v == 0){
    v = 10
  }
  
  dataa[i] = (dataa[i] - mu_sam[v])/sd_sam[v]
}
par(mfrow=c(2,3))
t_plot_1 = qqplot.ts(simulation_data[[1]])
p_acf_1 = qqplot.acf(dataa, "Poisson+PAR(1)")
p_pacf_1 = qqplot.pacf(dataa)

model_select = model[[3]]
par_int = par_int_ss[[3]]
simulation_data = simulation_model(model_select, par_int, NN, TT, L)
dataa = simulation_data[[1]]
mu_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = mean)
sd_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = sd)
for (i in 1:length(dataa)){
  v = i %% TT
  
  if (v == 0){
    v = 10
  }
  dataa[i] = (dataa[i] - mu_sam[v])/sd_sam[v]
}
t_plot_2 = qqplot.ts(simulation_data[[1]])
p_acf_2 = qqplot.acf(dataa, "Poisson+SAR(1)")
p_pacf_2 = qqplot.pacf(dataa)


pdf("boxPlot_Simulation_trajectory_Poisson.pdf", width=10, height=6, paper='special')
grid.arrange(arrangeGrob(t_plot_1, p_acf_1, p_pacf_1, ncol = 3), 
             arrangeGrob(t_plot_2, p_acf_2, p_pacf_2, ncol = 3), 
             nrow = 2)
dev.off()

TT = 10
L = 100
model_select = model[[5]]
par_int = par_int_ss[[5]]
simulation_data = simulation_model(model_select, par_int, NN, TT, L)
dataa = simulation_data[[1]]
#standardrize
mu_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = mean)
sd_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = sd)
for (i in 1:length(dataa)){
  v = i %% TT
  
  if (v == 0){
    v = 10
  }
  
  dataa[i] = (dataa[i] - mu_sam[v])/sd_sam[v]
}
par(mfrow=c(2,3))
t_plot_1 = qqplot.ts(simulation_data[[1]])
p_acf_1 = qqplot.acf(dataa, "TSMC+PAR(1)")
p_pacf_1 = qqplot.pacf(dataa)


model_select = model[[6]]
par_int = par_int_ss[[6]]
simulation_data = simulation_model(model_select, par_int, NN, TT, L)
dataa = simulation_data[[1]]
mu_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = mean)
sd_sam = apply(matrix(dataa, nrow = 10), MARGIN = 1, FUN = sd)
for (i in 1:length(dataa)){
  v = i %% TT
  
  if (v == 0){
    v = 10
  }
  
  dataa[i] = (dataa[i] - mu_sam[v])/sd_sam[v]
}
t_plot_2 = qqplot.ts(simulation_data[[1]])
p_acf_2 = qqplot.acf(dataa, "TSME+SAR(1)")
p_pacf_2 = qqplot.pacf(dataa)

pdf("boxPlot_Simulation_trajectory_TSMC.pdf", width=10, height=6, paper='special')
grid.arrange(arrangeGrob(t_plot_1, p_acf_1, p_pacf_1, ncol = 3), 
             arrangeGrob(t_plot_2, p_acf_2, p_pacf_2, ncol = 3), 
             nrow = 2)
dev.off()




simulation_learn_P = function(model_select, par_int, NN, TT, L, N, simulate_data){
  
  simulation_MLE <- main(simulate_data, model_select, par_int, N, TT)
  
  return(simulation_MLE)
}
simulation_P = function(par_star, X, model_select, org, TT, N, PP){ 
  
  Rp = matrix(runif((length(X)+PP) * N), nrow = N)
  rep.col<-function(x,n){matrix(rep(x,each=n), ncol=n, byrow=TRUE)}
  rep.row<-function(x,n){matrix(rep(x,each=n), nrow=n)}
  
  model_output = model_select(par_star, org, TT)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  Test_Pass = model_output$test_pass
  lag = model_output$lag
  
  if (model_output$S == FALSE){
    
    N = dim(Rp)[1]
    Z = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    w = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    epsilon = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    log_w = matrix(rep(0, N*length(X)), nrow = N)
    P = matrix(rep(0, PP*N), nrow = N)
    
    for (t in 1:lag){
      P_C_1 = p(t)$PC_1
      phi_1 = phi(t)
      Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
      P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
      if (t == 1){
        log_w[,t] = log(P_x1) 
      }else{
        log_w[,t] = log_w[,t-1] + log(P_x1)
      }
    }
    
    for (t in (lag+1):length(X)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[,t] = log_w[,t - 1] + log(w_Z)
      
    }
    log_wt = log_w[,t]
    
    for (t in (length(X)+1):(length(X)+PP)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v
      
      ##step 2, find epsolon##
      epsilon[,t] = qnorm(Rp[,t], mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      Z_tem = rep.col(pnorm(Z[,t]), length(p(v)$PC_2))
      P_tem = rep.row(p(v)$PC_2, length(Z[,t]))
      
      P[,t-length(X)] = apply(Z_tem > P_tem, MARGIN = 1, FUN = sum) - 1
    }
    
  }else{
    
    phi = model_output$phi
    a = model_output$a
    
    N = dim(Rp)[1]
    Z = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    w = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    epsilon = matrix(rep(0, N*(length(X)+PP)), nrow = N)
    log_w = matrix(rep(0, N*length(X)), nrow = N)
    P = matrix(rep(0, PP*N), nrow = N)
    
    P_C_1 = p(1)$PC_1
    phi_1 = phi
    Z[,1] = qtruncnorm(Rp[,1], a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[,1] = log(P_x1)
    
    
    for (t in 2:(lag+1)){
      
      CC = outer(1:t, 1:t, FUN = Cov_function_SAR, TT = TT, a = a, phi = phi)
      
      L = CC[1:(t-1),1:(t-1)]
      L_L = chol(L)
      
      D = forwardsolve(t(L_L), CC[1:(t-1),t])
      D_mu = apply(t(backsolve(L_L, D) * t(Z[,1:t-1])), MARGIN = 1, sum)
      
      
      D_sig = CC[t,t] - sum(D^2)
      
      P_C_1 = p(t %% TT)$PC_1
      phi_1 = phi
      Z[,t] = qtruncnorm(Rp[,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2],
                         mean = D_mu, sd = sqrt(D_sig))
      P_x1 = pnorm(P_C_1[X[t] + 2], mean = D_mu, sd = sqrt(D_sig)) -
        pnorm(P_C_1[X[t] + 1], mean = D_mu, sd = sqrt(D_sig))
      
      log_w[, t] = log_w[, t-1] + log(P_x1)
      
    }
    
    
    for (t in (lag+2):length(X)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi
      sigma_v = sqrt((1-a^2)*(1-phi^2)*(1-phi*a^TT)/(1+phi*a^TT))
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v + a*Z[,t-1] - a*phi_v*Z[,t-lag-1]
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[,t] = qtruncnorm(Rp[,t], a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[,t] = log_w[,t - 1] + log(w_Z)
      
    }
    
    log_wt = log_w[,t]
    
    for (t in (length(X)+1):(length(X)+PP)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi
      sigma_v = sqrt((1-a^2)*(1-phi^2)*(1-phi*a^TT)/(1+phi*a^TT))
      
      ##step 1, find Z_hat##
      Z_hat = Z[,t-lag] * phi_v + a*Z[,t-1] - a*phi_v*Z[,t-lag-1]
      
      ##step 2, find epsolon##
      epsilon[,t] = qnorm(Rp[,t], mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[,t] = Z_hat + sigma_v * epsilon[,t]
      
      Z_tem = rep.col(pnorm(Z[,t]), length(p(v)$PC_2))
      P_tem = rep.row(p(v)$PC_2, length(Z[,t]))
      
      P[,t-length(X)] = apply(Z_tem > P_tem, MARGIN = 1, FUN = sum) - 1
    }
    
  }
  
  D_tem = rep.col(exp(log_wt), PP)
  P_res = apply(D_tem*P, MARGIN = 2, FUN = mean)/
    apply(D_tem, MARGIN = 2, FUN = mean)
  
  return(list(P=P, log_wt = log_wt, P_res=P_res))
}
main_3 = function(model_select, par_int, NN, TT, L, N, PP){
  
  simulation_data = simulation_model(model_select, par_int, 1, TT, L+PP)
  simulation_MLE = simulation_learn_P(model_select, par_int, NN, TT, L, N, simulation_data[[1]][1:L])
  result_P = simulation_P(simulation_MLE$par_est, simulation_data[[1]][1:L], model_select, 
                          TRUE, TT, N, PP)
  
  return(list(simulation_data = simulation_data, 
              simulation_MLE = simulation_MLE, result_P = result_P))
}

dd = main_3(model[[1]], par_int_ss[[1]], 1, 10, 100, 500, 10)
dd$result_P$P_res
dd$simulation_data[[1]][101:110]

dd2 = main_3(model[[5]], par_int_ss[[5]], 1, 10, 100, 500, 10)
dd2$result_P$P_res
dd2$simulation_data[[1]][101:110]

pdf("boxPlot_Simulation_Prediction.pdf", width=10, height=8.5, paper='special')
par(mfrow = c(2,1))
plot(c(dd$simulation_data[[1]],max(dd$simulation_data[[1]]+2)), type = "n",
     xlab = "time_index",
     ylab = "Values",
     main = "Poisson+PAR(1)")
points(1:110, dd$simulation_data[[1]][1:110], col = "blue", pch = 19, cex = .7)
points(1:110, dd$simulation_data[[1]][1:110], type = 'l', lty = 2)
abline(v = 100, lty = 3)
#points(101:110, dd$simulation_data[[1]][101:110], col = "blue", pch = 3, cex = .7)
#points(101:110, dd$simulation_data[[1]][101:110], col = "blue", type = 'l', lty = 2)
points(101:110, dd$result_P$P_res, col = "red", pch = 4, cex = .7)
points(101:110, dd$result_P$P_res, col = "red", type = 'l', lty = 2)
legend("topleft", bty = "n",
       legend = c("Observation","Prediction"),
       pch = c(19, 4),
       col = c("blue", "red"))

plot(c(dd2$simulation_data[[1]],max(dd2$simulation_data[[1]]+2)), type = "n",
     xlab = "time_index",
     ylab = "Values",
     main = "TSMC+PAR(1)")
points(1:110, dd2$simulation_data[[1]][1:110], col = "blue", pch = 19, cex = .7)
points(1:110, dd2$simulation_data[[1]][1:110], type = 'l', lty = 2)
abline(v = 100, lty = 3)
#points(101:110, dd2$simulation_data[[1]][101:110], col = "green", pch = 3, cex = .7)
#points(101:110, dd2$simulation_data[[1]][101:110], col = "green", type = 'l', lty = 2)
points(101:110, dd2$result_P$P_res, col = "red", pch = 4, cex = .7)
points(101:110, dd2$result_P$P_res, col = "red", type = 'l', lty = 2)
legend("topleft", bty = "n",
       legend = c("Observation","Prediction"),
       pch = c(19, 4),
       col = c("blue", "red"))
dev.off()







