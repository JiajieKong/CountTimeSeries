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
library(lubridate)
library(truncnorm)
library(matrixStats)
library(numDeriv)


try(setwd(paste(getwd(),"/Dropbox/Robert/Wiley", sep = "")), silent = TRUE)
try(setwd(paste(getwd())), silent = TRUE)
# setwd("/Users/jiajiekong/dropbox/Robert/" ) #Mac
# setwd("/home/jiajie/Dropbox/Robert/" ) #Linux
# setwd("C:/Users/JiajieKong/Dropbox/Robert/Wiley") #Windows

# X = dataa #from another file
source("Rainy_Day_Seattle.R")

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

Possion_Reg = function(lambda){
  
  prob_0_to_100 = dpois(1:100)
  
  const = sum(prob_0_to_7)
  
  prob_0_to_100_Truncated = prob_0_to_7 / const
  prob_0_to_7_Truncated[8] = 1 - sum(prob_0_to_7_Truncated[1:7])
  prob_0_to_7_Truncated
  
}


model = list()

model[[1]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] 
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1))
  
} ##Bin+AR(1)

model[[2]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] 
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = TT))
  
} ##Bin+AR(1)

model[[3]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
  
} ##Bin + PAR(1)

model[[4]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)))
    test_pass = TRUE
  } else {
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
} ##TM+AR(1)

model[[5]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)))
    test_pass = TRUE
  } else {
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = TT))
} ##TM+SAR(1)

model[[6]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)),
                 tan((par[8] * pi / 2)),
                 tan((par[9] * pi / 52 - pi/2)))
  }else{
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
} ##TM+PAR(1)

para_int = list()
para_int[[1]] = c(0.436, 0.225, 3.437, .1)
para_int[[2]] = c(0.436, 0.225, 3.437, .1)
para_int[[3]] = c(.436, .225, 3.437, .238, 0, 10)
para_int[[4]] = c(.436, .225, 3.437, .436, .225, 3.437, .1)
para_int[[5]] = c(.436, .225, 3.437, .436, .225, 3.437, .1)
para_int[[6]] = c(.704, -.122, 3.216, .645, .159, 3.246, .238, 0, 10)


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
  
  N = dim(Rp)[1]
  W_I = rep(0, N)
  Z = rep(0, length(X))
  w = rep(0, length(X))
  epsilon = rep(0, length(X))
  log_w = rep(0, length(X))
  log_wt = rep(0, N)
  
  Path_PF = function(i){
    
    ##initiate 0##
    for (t in 1:lag){
      P_C_1 = p(t)$PC_1
      phi_1 = phi(t)
      Z[t] = qtruncnorm(Rp[i,t], a = P_C_1[X[t] + 1], b = P_C_1[X[t] + 2], mean = 0, sd = 1)
      P_x1 = pnorm(P_C_1[X[t] + 2]) - pnorm(P_C_1[X[t] + 1])
      log_w[t] = log(P_x1)
    }
    
    log_w[t] = sum(log_w)
    
    for (t in (lag+1):length(X)){
      
      v = t %% TT
      P_C = p(v)$PC_1
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      ##step 1, find Z_hat##
      Z_hat = Z[t-lag] * phi_v
      
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
    
    log_w[length(X)]
  }
  
  mc <- getOption("mc.cores", 4)
  log_wt = unlist(mclapply(1:100, Path_PF, mc.cores = mc))
  
  # log_wt = unlist(lapply(1:100, Path_PF))
  
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
                              X = X, model_select = model_select, Rp = Rp, TT = TT)
  })
  
  MLE_1
}


MLE = list()
res = list()


##simulation_begin_##
#Model_1_Bin and AR(1)
simulation_model = function(model_select, par, NN, TT, L){ #NN = 200
  
  model_output = model_select(par, org = TRUE, TT)
  p = model_output$p
  phi = model_output$phi
  lag = model_output$lag
  
  X_simulate = matrix(rep(0, NN * L), nrow = NN)
  Z_simulate = matrix(rep(0, NN * L), nrow = NN)
  
  for (i in 1:NN){
    
    for (j in 1:lag){
      Z_simulate[i,j] = rnorm(1)
      X_simulate[i,j] = sum(pnorm(Z_simulate[i,j]) > p(j)$PC_2 ) - 1
    }
    
    for (j in (lag+1):L){
      
      v = j %% TT
      phi_v = phi(v)
      sigma_v = sqrt(1 - phi_v^2)
      
      Z_simulate[i,j] = phi_v * Z_simulate[i,j-lag] + sigma_v * rnorm(1)
      X_simulate[i,j] = sum(pnorm(Z_simulate[i,j]) > p(v)$PC_2 ) - 1
    }
    
  }
  
  return(X_simulate)
}

simulation_learn = function(model_select, par_int, NN, TT, L, N, simulate_data){
  
  simulation_MLE = matrix(rep(0, 2 * NN * length(par_int)), ncol = 2 * length(par_int))
  
  for (j in 1:NN){
    
    cc = main(simulate_data[j,], model_select, par_int, N = 1e2, TT = TT)
    simulation_MLE[j,] = c(cc$par_est, sqrt(diag(solve(-cc$hessian_2))))
    
  }
  
  return(simulation_MLE)
}

simulation_data = list()
simulation_MLE = list()

for (j in 1:6) {
  
  simulation_data[[j]] = simulation_model(model[[j]], par_int_ss[[j]], 1e2, 52, 104)
  simulation_MLE[[j]] = simulation_learn(model[[j]], par_int_simu[[j]], 1e2, 52, 104, 1e2, 
                                         simulation_data[[j]])
  
}

model = list()

model[[1]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] 
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = 1))
  
} ##Bin+AR(1)

model[[2]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
    
    probability = par[1] + par[2] * cos((2*pi*(v - par[3]))/ TT)
    
    C = pbinom(0:7, size = 7, prob = probability) ##return C_0 to C_7
    
    PC_1 = c(-Inf, qnorm(C))
    PC_2 = c(0, C)
    return(list(PC_1 = PC_1, PC_2 = PC_2))
  }
  
  phi = function(v){
    
    par[4] 
    
  }
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran=par_tran, lag = TT))
  
} ##Bin+AR(1)

model[[3]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi / 2)),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)))
    test_pass = TRUE
  } else{
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
  
} ##Bin + PAR(1)

model[[4]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)))
    test_pass = TRUE
  } else {
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
} ##TM+AR(1)

model[[5]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)))
    test_pass = TRUE
  } else {
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = TT))
} ##TM+SAR(1)

model[[6]] = function(par_star, org, TT){
  
  if (org == TRUE){
    par_tran = par = par_star
    par_tran = c(tan((par[1] * pi) - pi/2),
                 tan((par[2] * pi / 2)),
                 tan((par[3] * pi / 52 - pi/2)),
                 tan((par[4] * pi) - pi/2),
                 tan((par[5] * pi / 2)),
                 tan((par[6] * pi / 52 - pi/2)),
                 tan((par[7] * pi / 2)),
                 tan((par[8] * pi / 2)),
                 tan((par[9] * pi / 52 - pi/2)))
  }else{
    par_tran = par = par_star
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
  
  return(list(par = par, p = p, phi = phi, test_pass = test_pass, par_tran = par_tran, lag = 1))
} ##TM+PAR(1)

para_int = list()
para_int[[1]] = c(0.436, 0.225, 3.437, .1)
para_int[[2]] = c(0.436, 0.225, 3.437, .1)
para_int[[3]] = c(.436, .225, 3.437, .238, 0, 10)
para_int[[4]] = c(.436, .225, 3.437, .436, .225, 3.437, .1)
para_int[[5]] = c(.436, .225, 3.437, .436, .225, 3.437, .1)
para_int[[6]] = c(.704, -.122, 3.216, .645, .159, 3.246, .238, 0, 10)


model_simulate = list()
model_simulate[[1]] = function(par_star, org, TT){  ##Poi+SAR(1)
  
}



par_int_ss = list()

par_int_ss[[1]] = c(.5, .2, 20, .3)

par_int_ss[[2]] = c(.5, .2, 20, .3, .2, 10)

par_int_ss[[3]] = c(.5, .2, 20, .4, .1, 10, .3)

par_int_ss[[4]] = c(.5, .2, 20, .4, .1, 10, .3, .15, 10)

par_int_ss[[5]] = c(3, 1, 10, 4.5, 1.5, 15, .2)

par_int_ss[[6]] = c(3, 1, 10, 4.5, 1.5, 15, .2, .1, 10)

par_int_ss[[7]] = c(10, 4, 20, .3, .2, 10)

par_int_ss[[8]] = c(10, -4, 20, -.3, .2, 10)

tt = readRDS("simulation_poisson")


pdf(paste("plot_simulation_1", '.pdf', sep = ""), width=8,height=2, paper='special')
par(mfrow=c(1,6))
dk = c("a1", "a2", "a3", "b1", "b2", "b3")
for (j in 1:6){
  
  hist(tt[[7]][,j], col = "lightblue", probability = TRUE, breaks = 10,
       main = dk[j], xlab = paste("mean=",round(mean(tt[[7]][,j]),2), "; sd=", round(sd(tt[[7]][,j]),2), sep = ""))
  abline(v = mean(tt[[7]][,j]), col = "red", lty = 2)
  
}
dev.off()


pdf(paste("plot_simulation_2", '.pdf', sep = ""), width=8,height=2, paper='special')
par(mfrow=c(1,6))
dk = c("a1", "a2", "a3", "b1", "b2", "b3")
for (j in 1:6){
  
  hist(tt[[8]][,j], col = "lightblue", probability = TRUE, breaks = 10,
       main = dk[j], xlab = paste("mean=",round(mean(tt[[8]][,j]),2), "; sd=", round(sd(tt[[8]][,j]),2), sep = ""))
  abline(v = mean(tt[[8]][,j]), col = "red", lty = 2)
  
}
dev.off()


