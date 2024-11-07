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

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# try(setwd(paste(getwd(),"/Dropbox/Robert", sep = "")), silent = TRUE)
# try(setwd(paste(getwd())), silent = TRUE)
# setwd("/Users/jiajiekong/dropbox/Robert/" ) #Mac
# setwd("/home/jiajie/Dropbox/Robert/" ) #Linux
# setwd("C:/Users/JiajieKong/Dropbox/Robert") #Windows

# X = dataa #from another file
source("Rainy_Day_Seattle.R")
X = X[1:104]

two_stat_MC = function(a, b){
  
  p_0_l = matrix(rep(0, 7*10), nrow = 7)
  p_1_l = matrix(rep(0, 7*10), nrow = 7)
  
  p_0_l[1,2] = a
  p_0_l[1,3] = 1 - a
  p_1_l[1,2] = 1 - b 
  p_1_l[1,3] = b
  
  for (l in 2:7){
    for (k in 0:l){
      
      p_0_l[l, k + 2] = a * p_0_l[l - 1, k + 2] + (1 - a) * p_1_l[l - 1, k + 1]
      p_1_l[l, k + 2] = b * p_1_l[l - 1, k + 1] + (1 - b) * p_0_l[l - 1, k + 2]
    }
  }
  
  return(list(p_0 = p_0_l, p_1 = p_1_l))
}

General_Possion_Truncated = function(mu, sig_2){
  
  lambda = 1 - sqrt(mu / sig_2)
  theta = mu * sqrt(mu / sig_2)
  
  prob_GP = function(x){
    
    theta * (theta + lambda * x)^(x - 1) * exp(-theta - lambda*x) / factorial(x)
    
  }
  
  prob_0_to_7 = prob_GP(0:7)
  
  const = sum(prob_0_to_7)
  
  prob_0_to_7_Truncated = prob_0_to_7 / const
  prob_0_to_7_Truncated[8] = 1 - sum(prob_0_to_7_Truncated[1:7])
  prob_0_to_7_Truncated
  
}

obj_function = function(par_star, X, model_select, Rp, org){ 
  
  model_output = model_select(par_star, org)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par
  Test_Pass = model_output$test_pass
  
  if (Test_Pass == FALSE){
    
    FF = -1e22
    par = round(par, 3)
    print(paste(FF, par[1], par[2], par[3], par[4]))
    return(FF)
    
  }
  
  TT = 52
  N = dim(Rp)[1]
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  for (i in 1:N){
    
    ##initiate 0##
    P_C_1 = p(1)$PC_1
    phi_1 = phi(1)
    Z[1] = qtruncnorm(Rp[i,1], a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    w[1] = 1
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[1] = 0
    
    for (t in 2:length(X)){
      
      v = t %% 52
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[t-1] * phi_v
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[t] = qtruncnorm(Rp[i,t], a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[t] = Z_hat + sigma_v * epsilon[t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[t] = log_w[t - 1] + log(w_Z)
      
      
    }
    
    log_wt[i] = log_w[length(X)]
    
  }
  
  par = round(par, 3)
  
  FF = log(P_x1) - log(N) + logSumExp(log_wt)
  
  print(paste(FF, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9]))
  
  FF
}

obj_function_Model_Check = function(par_star, X, model_select, org){ #par is vector of 6
  
  model_output = model_select(par_star, org)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par

  TT = 52
  N = 50
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  P_y = array(rep(0, N * length(X) * 8), dim = c(8, length(X), N))
  
  for (i in 1:N){
    
    ##initiate 0##
    P_C_1 = p(1)$PC_1
    phi_1 = phi(1)
    Z[1] = rtruncnorm(1, a = P_C_1[X[1] + 1], b = P_C_1[X[1] + 2], mean = 0, sd = 1)
    w[1] = 1
    P_x1 = pnorm(P_C_1[X[1] + 2]) - pnorm(P_C_1[X[1] + 1])
    log_w[1] = 0
    
    Z_hat = 0 ##marginal expectation for the first Obs
    P_y[, 1, i] = pnorm((P_C_1[2:9] - Z_hat)/ 1) - pnorm((P_C_1[1:8] - Z_hat)/ 1)
    
    for (t in 2:length(X)){
      
      v = t %% 52
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[t-1] * phi_v
      
      ##step 2, find epsolon##
      lowbound = (P_C[X[t] + 1] - Z_hat) / sigma_v
      upbound  = (P_C[X[t] + 2] - Z_hat) / sigma_v
      epsilon[t] = rtruncnorm(1, a = lowbound, b = upbound, mean = 0, sd = 1)
      
      ##step 3, update Z##
      Z[t] = Z_hat + sigma_v * epsilon[t]
      
      ##step 4, update w##
      w_Z = pnorm(upbound) - pnorm(lowbound)
      # w[t] = w[t - 1] * w_Z
      log_w[t] = log_w[t - 1] + log(w_Z)
      
      ##This part for PIT use##
      for (k in 0:7){
        
        lowbound = (P_C[k + 1] - Z_hat) / sigma_v
        upbound  = (P_C[k + 2] - Z_hat) / sigma_v
        P_y[k + 1, t, i] = pnorm(upbound) - pnorm(lowbound)
        
      }
    }
  }
  
  P_y_hat = matrix(rep(0, 8 * length(X)), nrow = 8)
  
  for (j in 1:N){
    
    P_y_hat = P_y_hat + P_y[,,j]
    
  }
  
  P_y_hat = apply(P_y_hat / N, MARGIN = 2, cumsum)
  P_y_hat[8,] = 1
  P_y_hat = rbind(rep(0, length(X)), P_y_hat)
  
  return(P_y_hat = P_y_hat)
}

res_calculation = function(X, par, model_select, org){
  
  TT = 52
  
  model_output = model_select(par, org)
  p = model_output$p
  phi = model_output$phi
  par = model_output$par

  z_hat_cal = function(X, par, p, phi){
    
    Z_t = rep(0, length(X))
    
    for (t in 1:length(X)){
      
      v = t %% 52
      P_C = p(v)$PC_1
      P_C_2 = p(v)$PC_2
      Z_t[t] = (exp(-P_C[X[t] + 1]^2/2) - exp(-P_C[X[t] + 2]^2/2))/
        (sqrt(2*pi) * (P_C_2[X[t] + 2] - P_C_2[X[t] + 1]))
      
    }
    
    epsilon_t = rep(0, length(X) - 1)

    for (t in 2:length(X)){
  
      epsilon_t[t] =  (Z_t[t] - Z_t[t-1] * phi(v)) / sqrt(1 - phi(v)^2)

    }
    
    return(list(Z_hat = Z_t, epsilon_t_hat = epsilon_t))
  }
  
  z_hat_cal(X, par, p, phi)
  
}

F_bar = function(u, P_y_hat, X){
  
  F_t = function(u, y, t, P_y_hat){
    
    if (u <= P_y_hat[y + 1,t]){
      F_t_u_y = 0
    }else if (u >= P_y_hat[y+2,t]){
      F_t_u_y = 1
    }else{
      F_t_u_y = (u - P_y_hat[y+1, t]) / (P_y_hat[y+2, t] - P_y_hat[y+1, t])
    }
    
    F_t_u_y
  }
  
  kk = 0
  
  for (t in 1:length(X)){
    
    kk = F_t(u, X[t], t, P_y_hat) + kk
    
  }
  
  kk / (length(X))
}

model = list()

model[[1]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  } else{
    par = rep(0, 4)
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52
    par[4] = (atan(par_star[4])) / pi * 2
    
    if ((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1)){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
      
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
      
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
      
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
    
  phi = function(v){
      
    par[4] 
      
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
  
}

model[[2]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  } else{
    par = rep(0, 4)
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52
    par[4] = (atan(par_star[4])) / pi * 2
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * 52
    
    if (((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1)) ||
        ((abs(par[4] + par[5]) > 1) || (abs(par[4] - par[5]) > 1))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
    
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] + par[5] * cos((2*pi*(v - par[6]))/ 52)
    
  }
  

  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
  
}

model[[3]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  } else {
    par = rep(0, 7)
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * 52
    par[7] = (atan(par_star[7])) / pi * 2
    
    if (((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1)) || 
        ((abs(par[4] + par[5]) > 1) || (abs(par[4] - par[5]) > 1))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  p = function(v){
    
    a = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
    b = par[4] + par[5] * cos((2*pi*(v - par[6]))/ 52)
    
    two_stat_MC_result = two_stat_MC(a, b)
    two_stat_MC_den = a/(a + b)*two_stat_MC_result$p_0[7, 2:9] + 
      b/(a + b)*two_stat_MC_result$p_1[7, 2:9]
    
    C = cumsum(two_stat_MC_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:7]), Inf)
    PC_2 = c(0, C[1:7], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[7] 
    
  }

  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
}

model[[4]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  }else{
    par = rep(0, 9)
    par[1] = (atan(par_star[1]) + pi/2) / pi * 1
    par[2] = (atan(par_star[2])) / pi * 2
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52
    par[4] = (atan(par_star[4]) + pi/2) / pi * 1
    par[5] = (atan(par_star[5])) / pi * 2
    par[6] = (atan(par_star[6]) + pi/2) / pi * 52
    par[7] = (atan(par_star[7])) / pi * 2
    par[8] = (atan(par_star[8])) / pi * 2
    par[9] = (atan(par_star[9]) + pi/2) / pi * 52
    
    if ((((abs(par[1] + par[2]) > 1) || (abs(par[1] - par[2]) > 1)) || 
         ((abs(par[4] + par[5]) > 1) || (abs(par[4] - par[5]) > 1))) ||
        ((abs(par[7] + par[8]) > 1) || (abs(par[7] - par[8]) > 1))){
      
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  
  p = function(v){
    
    a = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
    b = par[4] + par[5] * cos((2*pi*(v - par[6]))/ 52)
    
    two_stat_MC_result = two_stat_MC(a, b)
    two_stat_MC_den = a/(a + b)*two_stat_MC_result$p_0[7, 2:9] + 
      b/(a + b)*two_stat_MC_result$p_1[7, 2:9]
    
    C = cumsum(two_stat_MC_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:7]), Inf)
    PC_2 = c(0, C[1:7], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[7] + par[8] * cos((2*pi*(v - par[9]))/ 52)
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
}

model[[5]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  }else{
    par = rep(0, 7)
    par[1] = exp(par_star[1])
    par[2] = exp(par_star[2])
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52 #
    par[4] = exp(par_star[4])
    par[5] = exp(par_star[5])
    par[6] = (atan(par_star[6]) + pi/2) / pi * 52 #
    par[7] = (atan(par_star[7])) / pi * 2
    
    if ((((par[1] + par[2]) < 0) || ((par[1] - par[2]) < 0)) || 
        (((par[4] + par[5]) < 0) || ((par[4] - par[5]) < 0))){
      # (sum(GPT_test(0:51)) != 52)){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }

  
  p = function(v){
    
    mu = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
    sig_2 = par[4] + par[5] * cos((2*pi*(v - par[6]))/ 52)
    
    GPT_den = General_Possion_Truncated(mu, sig_2)
    
    C = cumsum(GPT_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:7]), Inf)
    PC_2 = c(0, C[1:7], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[7] 
    
  }

  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
}

model[[6]] = function(par_star, org){
  
  if (org == TRUE){
    par = par_star
    test_pass = TRUE
  }else{
    par = rep(0, 7)
    par[1] = exp(par_star[1])
    par[2] = exp(par_star[2])
    par[3] = (atan(par_star[3]) + pi/2) / pi * 52 #
    par[4] = exp(par_star[4])
    par[5] = exp(par_star[5])
    par[6] = (atan(par_star[6]) + pi/2) / pi * 52 #
    par[7] = (atan(par_star[7])) / pi * 2
    par[8] = (atan(par_star[8])) / pi * 2
    par[9] = (atan(par_star[9]) + pi/2) / pi * 52
    
    if ((((par[1] + par[2]) < 0) || ((par[1] - par[2]) < 0)) || 
        (((par[4] + par[5]) < 0) || ((par[4] - par[5]) < 0))){
      test_pass = FALSE
    }else{
      test_pass = TRUE
    }
  }
  
  
  p = function(v){
    
    mu = par[1] + par[2] * cos((2*pi*(v - par[3]))/ 52)
    sig_2 = par[4] + par[5] * cos((2*pi*(v - par[6]))/ 52)
    
    GPT_den = General_Possion_Truncated(mu, sig_2)
    
    C = cumsum(GPT_den) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C[1:7]), Inf)
    PC_2 = c(0, C[1:7], 1)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[7] + par[8] * cos((2*pi*(v - par[9]))/ 52)
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass))
}

para_int = list()

para_int[[1]] = c(tan((0.436 * pi) - pi/2),
             tan((0.225 * pi / 2)),
             tan((3.437 * pi / 52 - pi/2)),
             tan((0.1 * pi / 2)))

para_int[[2]] = c(tan((0.436 * pi) - pi/2),
                  tan((0.225 * pi / 2)),
                  tan((3.437 * pi / 52 - pi/2)),
                  tan((0.238 * pi / 2)),
                  tan((0 * pi / 2)),
                  tan((10 * pi / 52 - pi/2)))

para_int[[3]] = c(tan((0.436 * pi) - pi/2),
               tan((0.225 * pi / 2)),
               tan((3.437 * pi / 52 - pi/2)),
               tan((0.436 * pi) - pi/2),
               tan((0.225 * pi / 2)),
               tan((3.437 * pi / 52 - pi/2)),
               tan((0.1 * pi / 2)))

para_int[[4]] = c(tan((0.704 * pi) - pi/2),
               tan((-.122 * pi / 2)),
               tan((3.216 * pi / 52 - pi/2)),
               tan((0.645 * pi) - pi/2),
               tan((0.159 * pi / 2)),
               tan((3.246 * pi / 52 - pi/2)),
               tan((0.238 * pi / 2)),
               tan((0 * pi / 2)),
               tan((10 * pi / 52 - pi/2)))

para_int[[5]] = c(log(3),
               log(1),
               tan((3.437 * pi / 52 - pi/2)),
               log(4.5),
               log(1.2),
               tan((3.437 * pi / 52 - pi/2)),
               tan((0.1 * pi / 2)))

para_int[[6]] = c(log(3),
               log(1.5),
               tan((3.437 * pi / 52 - pi/2)),
               log(3.7),
               log(0.85),
               tan((3.6 * pi / 52 - pi/2)),
               tan((0.19 * pi / 2)),
               tan((0 * pi / 2)),
               tan((10 * pi / 52 - pi/2)))


# par[1] = (atan(par[1]) + pi/2) / pi * 1
# par[2] = (atan(par[2])) / pi * 2
# par[3] = (atan(par[3]) + pi/2) / pi * 52
# par[4] = (atan(par[4])) / pi * 2
# par[5] = (atan(par[5])) / pi * 2
# par[6] = (atan(par[6]) + pi/2) / pi * 52


main = function(X, model_select, para_int_select, N, TT){
 
  Rp = matrix(runif(TT * N), nrow = N)
  
  system.time({
    MLE_1 = optim(para_int_select, obj_function, 
                  X = X, model_select = model_select, Rp = Rp, org = FALSE,
                  method = "BFGS", control = list(factr = 1e7, fnscale=-1))
    
    MLE_1$par_est = model_select(MLE_1$par, org = FALSE)$par
    
    MLE_1$hessian_2 = hessian(func = obj_function, MLE_1$par_est, org = TRUE,
                              X = X, model_select = model_select, Rp = Rp)
  })
  
  MLE_1
}


Model_Checking = function(MLE_1, model_select){
  
  res = res_calculation(X, MLE_1$par, model_select, org = FALSE)
  res$P_y_hat = obj_function_Model_Check(MLE_1$par, X, model_select, org = FALSE)
  
  return(res)
}


plot_res = function(res){
  plot(res$epsilon_t_hat, type = 'l', main = 'Residual Series')
  hist(res$epsilon_t_hat, col = 'lightblue')
  acf(res$epsilon_t_hat)
  pacf(res$epsilon_t_hat)
  grid = seq(0, 1, length.out = 11)
  tt = rep(0, 11)
  for (j in 1:length(grid)){
    tt[j] = F_bar(grid[j], res$P_y_hat, X)
  }
  plot(grid[2:length(grid)], diff(tt), type = "S", ylim = c(0, .2), main = "Model Checking")
  abline(a = .1, b = 0, col = "red", lty = 2)
}


MLE = list()
res = list()


for (j in 1:6){
  
  MLE[[j]] = main(X, model[[j]], para_int[[j]], N = 5e2, TT = 104)
  res[[j]] = Model_Checking(MLE[[j]], model[[j]])
  
  MLE[[j]]$res = res[[j]]
  MLE[[j]]$se = sqrt(diag(solve(-MLE[[j]]$hessian_2)))
  
}


##begin to make comparison
pdf(paste("plot_result", '.pdf', sep = ""), 
    width=12,height=12, paper='special')
par(mfrow=c(6,5))
for (j in 1:6){
  plot_res(MLE[[j]]$res)
}
dev.off()

saveRDS(MLE, file = "result_2")


##simulation_begin_##
#Model_1_Bin and AR(1)
simulation_model = function(model_select, par, NN, TT){ #NN = 200
  
  model_output = model_select(par, org = TRUE)
  p = model_output$p
  phi = model_output$phi

  X_simulate = matrix(rep(0, NN * TT), nrow = NN)
  Z_simulate = matrix(rep(0, NN * TT), nrow = NN)
  
  for (i in 1:NN){
    
    Z_simulate[i,1] = rnorm(1)
    X_simulate[i,1] = sum(pnorm(Z_simulate[i,1]) > p(1)$PC_2 ) - 1
    
    for (j in 2:TT){
      
      v = j %% 52
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      Z_simulate[i,j] = phi_v * Z_simulate[i,j-1] + sigma_v * rnorm(1)
      X_simulate[i,j] = sum(pnorm(Z_simulate[i,j]) > p(v)$PC_2 ) - 1
    }

  }
  
  return(X_simulate)
}

simulation_learn = function(model_select, par_int, NN, TT, simulate_data){
  
  simulation_MLE = matrix(rep(0, NN * length(par_int)), ncol = length(par_int))
  
  for (j in 1:NN){
    
    cc = main(simulate_data[j,], model_select, par_int, N = 1e2, TT = TT)
    simulation_MLE[j,] = cc$par_est
    
  }
  
  return(simulation_MLE)
}

simulation_data = list()
simulation_MLE = list()

for (j in 1:6) {
  
  simulation_data[[j]] = simulation_model(model[[j]], par_int_ss[[j]], 1e2, 104)
  simulation_MLE[[j]] = simulation_learn(model[[j]], par_int_simu[[j]], 1e2, 104, 
                                         simulation_data[[j]])
  
}

simulation_data[[3]] = simulation_model(model[[3]], par_int_ss[[3]], 1e2, 104)
simulation_MLE[[3]] = simulation_learn(model[[3]], par_int_simu[[3]], 1e2, 104, 
                                       simulation_data[[3]])
simulation_data[[4]] = simulation_model(model[[4]], par_int_ss[[4]], 1e2, 104)
simulation_MLE[[4]] = simulation_learn(model[[4]], par_int_simu[[4]], 1e2, 104, 
                                       simulation_data[[4]])
simulation_data[[5]] = simulation_model(model[[5]], par_int_ss[[5]], 1e2, 104)
simulation_MLE[[5]] = simulation_learn(model[[5]], par_int_simu[[5]], 1e2, 104, 
                                       simulation_data[[5]])
simulation_data[[6]] = simulation_model(model[[6]], par_int_ss[[6]], 1e2, 104)
simulation_MLE[[6]] = simulation_learn(model[[6]], par_int_simu[[6]], 1e2, 104, 
                                       simulation_data[[6]])


par_int_simu = list()
par_int_ss = list()

par_int_ss[[1]] = c(.5, .2, 20, .3)
par_int_simu[[1]] = c(tan((0.5 * pi) - pi/2),
                      tan((0.2 * pi / 2)),
                      tan((20 * pi / 52 - pi/2)),
                      tan((0.3 * pi / 2)))

par_int_ss[[2]] = c(.5, .2, 20, .3, .2, 10)
par_int_simu[[2]] = c(tan((0.5 * pi) - pi/2),
                      tan((0.2 * pi / 2)),
                      tan((20 * pi / 52 - pi/2)),
                      tan((0.3 * pi / 2)),
                      tan((0.2 * pi / 2)),
                      tan((10 * pi / 52 - pi/2)))

par_int_ss[[3]] = c(.5, .2, 20, .4, .1, 10, .3)
par_int_simu[[3]] = c(tan((0.5 * pi) - pi/2),
                      tan((0.2 * pi / 2)),
                      tan((20 * pi / 52 - pi/2)),
                      tan((0.4 * pi) - pi/2),
                      tan((0.1 * pi / 2)),
                      tan((10 * pi / 52 - pi/2)),
                      tan((0.3 * pi / 2)))

par_int_ss[[4]] = c(.5, .2, 20, .4, .1, 10, .3, .15, 10)
par_int_simu[[4]] = c(tan((0.5 * pi) - pi/2),
                      tan((0.2 * pi / 2)),
                      tan((20 * pi / 52 - pi/2)),
                      tan((0.4 * pi) - pi/2),
                      tan((0.1 * pi / 2)),
                      tan((10 * pi / 52 - pi/2)),
                      tan((0.3 * pi / 2)),
                      tan((0.15 * pi / 2)),
                      tan((10 * pi / 52 - pi/2)))

par_int_ss[[5]] = c(3, 1, 10, 4.5, 1.5, 15, .2)
par_int_simu[[5]] = c(log(3),
                      log(1),
                      tan((10 * pi / 52 - pi/2)),
                      log(4.5),
                      log(1.5),
                      tan((15 * pi / 52 - pi/2)),
                      tan((0.2 * pi / 2)))

par_int_ss[[6]] = c(3, 1, 10, 4.5, 1.5, 15, .2, .1, 10)
par_int_simu[[6]] = c(log(3),
                      log(1.5),
                      tan((10 * pi / 52 - pi/2)),
                      log(4.5),
                      log(1.5),
                      tan((15 * pi / 52 - pi/2)),
                      tan((0.2 * pi / 2)),
                      tan((0.1 * pi / 2)),
                      tan((10 * pi / 52 - pi/2)))







