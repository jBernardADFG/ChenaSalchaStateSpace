# PROJECT OUT A FEW YEARS #
# SOME THOUGHTS -- WE CAN PROJECT THE NUMBER OF SPAWNERS/RECRUITS OUT A FEW YEARS #
# THINK ABOUT TOSSING IN A TERM FOR SEA-SURFACE TEMPERATURE #
# ADDITIVE VS LINEAR EFFECTS #

# HIERARCHICAL TERMS FOR ALPHA AND BETA ? 

# READ IN THE AND ORGANIZE DATA FOR USE IN JAGS MODEL #
{
  rm(list=objects())
  library(readxl)
  jordyData <- as.matrix(read_excel("D:/Jordy/ChenaSalchaSpawnerRecruit/Data/jordyData.xlsx"))
  n_years <- nrow(jordyData)
  n_ages <- 8
  log_mr_est <- log(jordyData[,2:3])
  log_tow_est <- log(jordyData[,6:7])
  mr_cv <- jordyData[,4:5]
  tow_cv <- jordyData[,8:9]
  for (y in 1:n_years){  
    for (r in 1:2){  
      if(is.na(mr_cv[y,r])){
        mr_cv[y,r] <- -999
      }
      if(is.na(tow_cv[y,r])){
        tow_cv[y,r] <- -999
      }
    }
  }
  log_h_riv_est <- log(jordyData[,10:11])
  log_h_my_est <- log(jordyData[,12]) 
  n_my_riv <- jordyData[,13:14]
  n_my_tot <- jordyData[,15]
  for (y in 1:length(n_my_tot)){
    if(is.na(n_my_tot[y])){
      n_my_tot[y] <- 0
    }
  }
  age_samples <- array(NA, dim=c(2,n_years,6))
  age_samples[1,,] <- jordyData[,16:21]
  age_samples[2,,] <- jordyData[,22:27]
  n_age_samples <- apply(age_samples[,,], c(2,1), sum, na.rm=T)
  # MEAN MONTHLY DISCHARGE #
  chena_water_data <- read.delim("D:/Jordy/ChenaSalchaSpawnerRecruit/Data/chena-at-fairbanks")
  salcha_water_data <- read.delim("D:/Jordy/ChenaSalchaSpawnerRecruit/Data/salcha-near-salchaket")
  get.monthly.discharge <- function(water_data, month, start_year, end_year){
    water_data <- water_data[,5:7]
    names(water_data) <- c("year", "month", "discharge")
    water_data <- water_data[water_data$month==month,]
    water_data <- water_data[,c(1,3)]
    water_data <- water_data[water_data$year>=start_year & water_data$year<=end_year,]
    water_data <- water_data[,2]
    water_data <- water_data - median(water_data)
    return(water_data)
  }
  chena_water_data <- get.monthly.discharge(chena_water_data, 8, 1986, 2007)
  salcha_water_data <- get.monthly.discharge(salcha_water_data, 8, 1986, 2007)
  water_level <- cbind(chena_water_data, salcha_water_data)
  s_max=5000
  data <- list(n_years=n_years,
               n_ages=n_ages,
               log_mr_est=log_mr_est,
               log_tow_est=log_tow_est,
               mr_cv=mr_cv,
               tow_cv=tow_cv,
               log_h_riv_est=log_h_riv_est,
               log_h_my_est=log_h_my_est,
               n_my_riv=n_my_riv,
               n_my_tot=n_my_tot,
               age_samples=age_samples,
               n_age_samples=n_age_samples,
               water_level=water_level)
  
}



# WRITE JAGS MODEL #


mod <-
  "model{
  
  ########################################################################
  ############################ LATENT PROCESS ############################
  ########################################################################
  
  # --------------
  # IN-RIVER-RUN-ABUNDANCE ON THE CHENA AND SALCHA DURING THE INITIAL YEARS #
  for (r in 1:2){
    for (y in 1:n_ages){
      IRRA[y,r] ~ dnorm(mu_e[r], tau_e[r])T(0,)
    }
    mu_e[r] ~ dunif(0,25000)
    tau_e[r] <- pow(1/sig_e[r], 2)
    sig_e[r] ~ dexp(1E-4)
  }
  
  # --------------
  # HARVEST ON THE CHENA AND SALCHA #
  for (r in 1:2){
    for (y in 1:n_years){                                                    
      H_riv[y,r] ~ dnorm(mu_h_riv[r], tau_h_riv[r])T(0,)
    }
    mu_h_riv[r] ~ dunif(0, 2000)
    tau_h_riv[r] <- pow(1/sig_h_riv[r], 2)
    sig_h_riv[r] ~ dexp(1E-4)
  }
  
  # --------------
  # SPAWNERS GIVEN IN-RIVER-RUN ABUNDANCE AND HARVEST ON THE CHENA AND SALCHA #
  for (r in 1:2){
    for (y in 1:n_years){                               
      S[y,r] <- max(IRRA[y,r]-H_riv[y,r], 1)
    }
  }
  
  
  # ------------------------------------------------- #
  # ----------------- RS PROCESSES  ----------------- #
  # ------------------------------------------------- #
  
  # ---------------- USE ONE OF THE FOLLOWING TEN OPTIONS ---------------- #
  
# --------------
# SIMPLE RICKER RS PROCESS #
for (r in 1:2){
  for (y in 1:n_years){
    log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
    mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - beta[r]*S[y,r]
    R[y,r] <- exp(log_R[y,r])
  }
  tau_w[r] <- pow(1/sig_w[r], 2)
  sig_w[r] ~ dexp(0.1)
  alpha[r] ~ dexp(1E-2)T(1,)
  log_alpha[r] <- log(alpha[r])
  beta[r] ~ dexp(1E2)
}

#  # --------------  
#  # RICKER RS PROCESS WITH AN AR(1) TERM #
#  for (r in 1:2){
#    log_R[1,r] ~ dnorm(mu_sr[1,r], tau_w[r])
#    mu_sr[1,r] <- log(alpha[r]) + log(S[1,r]) - beta[r]*S[1,r]
#    resid[1, r] <- 0
#    R[1,r] <- exp(log_R[1,r])
#    for (y in 2:n_years){
#      log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#      mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - beta[r]*S[y,r] + phi[r]*resid[y-1,r]
#      resid[y, r] <- log_R[y,r]-log(alpha[r])-log(S[y,r])+beta[r]*S[y,r]
#      R[y,r] <- exp(log_R[y,r])
#    }
#    tau_w[r] <- pow(1/sig_w[r], 2)
#    sig_w[r] ~ dexp(0.1)
#    alpha[r] ~ dexp(1E-2)T(1,)
#    log_alpha[r] <- log(alpha[r])
#    phi[r] ~ dunif(-1,1)
#    beta[r] ~ dexp(1E2)
#  }

# # --------------
# # RICKER RS PROCESS WITH A TIME VARYING PRODUCTIVITY PARAMETER #
#   for (r in 1:2){
#     for (y in 1:n_years){
#       log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#       mu_sr[y,r] <- log_alpha[y,r] + log(S[y,r]) - beta[r]*S[y,r]
#       log_alpha[y,r] <- c[r] + d[r]*(y-1)
#       alpha[y,r] <- exp(log_alpha[y,r])
#       R[y,r] <- exp(log_R[y,r])
#     }
#     tau_w[r] <- pow(1/sig_w[r], 2)
#     sig_w[r] ~ dexp(0.1)
#     c[r] ~ dexp(1E-2)
#     d[r] ~ dnorm(0, 1E-1)T(-c[r]/(n_years-1),)
#     beta[r] ~ dexp(1E2)
#   }

# # --------------
# # RICKER RS PROCESS WITH AN AR(1) TERM AND A TIME VARYING PRODUCTIVITY PARAMETER #
#   for (r in 1:2){
#     log_R[1,r] ~ dnorm(mu_sr[1,r], tau_w[r])
#     mu_sr[1,r] <- log_alpha[1,r] + log(S[1,r]) - beta[r]*S[1,r]
#     log_alpha[1,r] <- c[r]
#     alpha[1,r] <- exp(log_alpha[1,r])
#     R[1,r] <- exp(log_R[1,r])
#     resid[1, r] <- 0
#     for (y in 2:n_years){
#       log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#       mu_sr[y,r] <- log(alpha[y,r]) + log(S[y,r]) - beta[r]*S[y,r] + phi[r]*resid[y-1,r]
#       log_alpha[y,r] <- c[r] + d[r]*(y-1)
#       alpha[y,r] <- exp(log_alpha[y,r])
#       resid[y, r] <- log_R[y,r]-log(alpha[y, r])-log(S[y,r])+beta[r]*S[y,r]
#       R[y,r] <- exp(log_R[y,r])
#     }
#     tau_w[r] <- pow(1/sig_w[r], 2)
#     sig_w[r] ~ dexp(0.1)
#     c[r] ~ dexp(1E-2)
#     d[r] ~ dnorm(0, 1E-1)T(-c[r]/(n_years-1),)
#     phi[r] ~ dunif(-1,1)
#     beta[r] ~ dexp(1E2)
#   }

# # --------------
# # RICKER RS PROCESS WITH A TERM FOR THE WATER LEVEL #
# for (r in 1:2){
#   for (y in 1:n_years){
#     log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#     mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - beta[r]*S[y,r] + omega[r]*water_level[y,r]
#     R[y,r] <- exp(log_R[y,r])
#   }
#   tau_w[r] <- pow(1/sig_w[r], 2)
#   sig_w[r] ~ dexp(0.1)
#   alpha[r] ~ dexp(1E-2)T(1,)
#   omega[r] ~ dnorm(0, 1E-4)
#   log_alpha[r] <- log(alpha[r])
#   beta[r] ~ dexp(1E2)
# }

# # --------------
# # RICKER RS PROCESS WITH WATER LEVEL AND AR(1) TERMS #
# for (r in 1:2){
#   log_R[1,r] ~ dnorm(mu_sr[1,r], tau_w[r])
#   mu_sr[1,r] <- log(alpha[r]) + log(S[1,r]) - beta[r]*S[1,r] + omega[r]*water_level[1,r]
#   resid[1, r] <- 0
#   R[1,r] <- exp(log_R[1,r])
#     for (y in 2:n_years){
#     log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#     mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - beta[r]*S[y,r] + omega[r]*water_level[y,r] + phi[r]*resid[y-1,r]
#     resid[y, r] <- log_R[y,r]-log(alpha[r])-log(S[y,r])+beta[r]*S[y,r]-omega[r]*water_level[1,r]
#     R[y,r] <- exp(log_R[y,r])
#   }
#   tau_w[r] <- pow(1/sig_w[r], 2)
#   sig_w[r] ~ dexp(0.1)
#   alpha[r] ~ dexp(1E-2)T(1,)
#   omega[r] ~ dnorm(0, 1E-4)
#   phi[r] ~ dunif(-1,1)
#   log_alpha[r] <- log(alpha[r])
#   beta[r] ~ dexp(1E2)
# }

#  # --------------
#  # RICKER RS PROCESS WITH TIME VARYING PRODUCTIVITY AND A TERM FOR WATER LEVEL #
#  for (r in 1:2){
#    for (y in 1:n_years){
#      log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#      mu_sr[y,r] <- log_alpha[y,r] + log(S[y,r]) - beta[r]*S[y,r] + omega[r]*water_level[y,r]
#      log_alpha[y,r] <- c[r] + d[r]*(y-1)
#      alpha[y,r] <- exp(log_alpha[y,r])
#      R[y,r] <- exp(log_R[y,r])
#    }
#    tau_w[r] <- pow(1/sig_w[r], 2)
#    sig_w[r] ~ dexp(0.1)
#    c[r] ~ dexp(1E-2)
#    d[r] ~ dnorm(0, 1E-1)T(-c[r]/(n_years-1),)
#    omega[r] ~ dnorm(0, 1E-4)
#    beta[r] ~ dexp(1E2)
#  }

#  # --------------
#  # RICKER RS PROCESS WITH TIME VARYING PRODUCTIVITY AND AR(1) AND WATER LEVEL TERMS #
#  for (r in 1:2){
#    log_R[1,r] ~ dnorm(mu_sr[1,r], tau_w[r])
#    mu_sr[1,r] <- log_alpha[1,r] + log(S[1,r]) - beta[r]*S[1,r]  + omega[r]*water_level[1,r]
#    log_alpha[1,r] <- c[r]
#    alpha[1,r] <- exp(log_alpha[1,r])
#    R[1,r] <- exp(log_R[1,r])
#    resid[1, r] <- 0
#    for (y in 2:n_years){
#      log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#      mu_sr[y,r] <- log(alpha[y,r]) + log(S[y,r]) - beta[r]*S[y,r] + phi[r]*resid[y-1,r] + omega[r]*water_level[y,r]
#      log_alpha[y,r] <- c[r] + d[r]*(y-1)
#      alpha[y,r] <- exp(log_alpha[y,r])
#      resid[y, r] <- log_R[y,r] - log(alpha[y, r]) - log(S[y,r]) + beta[r]*S[y,r] - omega[r]*water_level[y,r]
#      R[y,r] <- exp(log_R[y,r])
#    }
#    tau_w[r] <- pow(1/sig_w[r], 2)
#    sig_w[r] ~ dexp(0.1)
#    c[r] ~ dexp(1E-2)
#    d[r] ~ dnorm(0, 1E-1)T(-c[r]/(n_years-1),)
#    phi[r] ~ dunif(-1,1)
#    omega[r] ~ dnorm(0, 1E-4)
#    beta[r] ~ dexp(1E2)
#  }

#  # --------------
#  # SIMPLE BEVERTON-HOLT RS PROCESS #
#  for (r in 1:2){
#    for (y in 1:n_years){
#      log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#      mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - log(1+beta[r]*S[y,r])
#      R[y,r] <- exp(log_R[y,r])
#    }
#    tau_w[r] <- pow(1/sig_w[r], 2)
#    sig_w[r] ~ dexp(0.1)
#    alpha[r] ~ dexp(1E-4)T(1,)
#    log_alpha[r] <- log(alpha[r])
#    beta[r] ~ dexp(1E-4)
#  }

#  # --------------
#  # BEVERTON-HOLT RS PROCESS WITH AN AR(1) TERM #
#  for (r in 1:2){
#    log_R[1,r] ~ dnorm(mu_sr[1,r], tau_w[r])
#    mu_sr[1,r] <- log(alpha[r]) + log(S[1,r]) - beta[r]*S[1,r]
#    resid[1, r] <- 0
#    for (y in 1:n_years){
#      log_R[y,r] ~ dnorm(mu_sr[y,r], tau_w[r])
#      mu_sr[y,r] <- log(alpha[r]) + log(S[y,r]) - log(1+beta[r]*S[y,r]) + phi[r]*resid[y-1,r]
#      resid[y, r] <- log_R[y,r]-log(alpha[r])-log(S[y,r])+beta[r]*S[y,r]
#      R[y,r] <- exp(log_R[y,r])
#    }
#    tau_w[r] <- pow(1/sig_w[r], 2)
#    sig_w[r] ~ dexp(0.1)
#    alpha[r] ~ dexp(1E-4)T(1,)
#    log_alpha[r] <- log(alpha[r])
#    phi[r] ~ dunif(-1,1)
#    beta[r] ~ dexp(1E-4)
#  }
  
  # ------------------------------------------------- #

  # --------------
  # RETURNERS GIVEN RECRUITS #
  for (r in 1:2){
    for (y in (n_ages+1):n_years){
      for (a in 1:6){
        A1[y,r,a] <- R[(y-9+a),r]*p_maturity[r,y,7-a]
      }
      Returners[y,r] <- sum(A1[y,r,1:6])
    }
  }
  
  
#  # ---------------- CHOOSE ONE OF THE FOLLOWING TWO OPTIONS  ---------------- #
#  # ----------------
#  # ----------------
#  # WITHOUT TIME VARYING AGE-AT-MATURITY VERSION 1 #
#  for (r in 1:2){
#    for (y in 1:n_years){
#      p_maturity[r,y,1:6] ~ ddirch(gamma[1:6]+0.1)
#    }
#  }
#  for (a in 1:n_ages){
#    gamma[a] ~ dexp(0.1)
#  }

#  # ----------------
#  # WITHOUT TIME VARYING AGE-AT-MATURITY VERSION 2 #
#  for (r in 1:2){
#    for (y in 1:n_years){
#      p_maturity[r,y,1:6] ~ ddirch(gamma[r,1:6]+0.1)
#    }
#  }
#  for (r in 1:2){
#    for (a in 1:n_ages){
#      gamma[r,a] ~ dexp(0.1)
#    }  
#  }

  # ----------------
  # WITH TIME VARYING AGE AT MATURITY #  
  for (r in 1:2){
    for (y in 1:n_years){
    p_maturity[r,y,1:6] ~ ddirch(gamma[r,y,1:6]+0.1)
      for (a in 1:6){
        gamma[r,y,a] <- pi[r,y,a]*D[r]
        pi[r,y,a] <- logistic[r,y,a]/sum(logistic[r,y,1:6])
        logistic[r,y,a] <- exp(n[1,r,a]+n[2,r,a]*y)
      }
    }
    D[r] ~ dexp(0.001)T(1,)
    for (a in 1:6){
      n[1,r,a] ~ dunif(-100, 100)
      n[2,r,a] ~ dunif((-100-n[1,r,a])/n_years, (100-n[1,r,a])/n_years)
    }
  }
  
  # ------------------------------------------------- #
  
  
  # --------------
  # PROBABILITY OF CHENA AND SALCHA HARVEST ON THE MIDDLE YUKON #
  for (r in 1:2){
    for (y in 1:n_years){
      p_my_harvest[y,r] ~ dbeta(a[r], b[r])
    }
    a[r] <- xi[r] + 1                     
    b[r] <- nu[r] - xi[r] + 1
    nu[r] ~ dexp(1E-4)
    xi[r] ~ dunif(0, nu[r])
  }
  
  # --------------
  # MIDDLE YUKON HARVEST #
  for (y in 1:n_years){
    H_my[y] ~ dnorm(mu_h_my, tau_h_my)T(0,)
  }
  mu_h_my ~ dunif(0, 50000)
  tau_h_my <- pow(1/sig_h_my, 2)
  sig_h_my ~ dexp(1E-5)
  
  # --------------
  # IRRA GIVEN RETURNERS AND MIDDLE YUKON HARVEST #
  for (r in 1:2){
     for (y in (n_ages+1):n_years){
       IRRA[y,r] <- max(Returners[y,r]-p_my_harvest[y,r]*H_my[y], 1) 
     }
  }
  
  
  #############################################################################
  ############################ OBSERVATION PROCESS ############################
  #############################################################################
  
  for(r in 1:2){
    for (y in 1:n_years){
    
      # --------------
      # MARK-RECAPTURE ABUNDANCE ESTIMATES #
      log_mr_est[y, r] ~ dnorm(log(IRRA[y,r]), tau_mr[y,r])
      tau_mr[y,r] <- 1/var_mr[y,r]
      var_mr[y,r] <- log(pow(mr_cv[y,r], 2)+1)
    
      # --------------
      # TOWER COUNTS #
      log_tow_est[y, r] ~ dnorm(log(IRRA[y,r]), tau_tow[y,r])
      tau_tow[y,r] <- 1/var_tow[y,r]
      var_tow[y,r] <- log(pow(tow_cv[y,r],2)+1)
      
      # --------------
      # CHENA AND SALCHA HARVEST #
      log_h_riv_est[y,r] ~ dnorm(log(H_riv[y,r]), tau_h_riv_est[r])
      
      # --------------
      # AGE DATA FROM THE CHENA AND SALCHA #
      age_samples[r, y, 1:6] ~ dmulti(p_maturity[r,y,1:6], n_age_samples[y,r])
      
      # --------------
      # MOVEMENT BETWEEN THE MIDDLE YUKON AND THE CHENA AND SALCHA #
      n_my_riv[y,r] ~ dbin(p_my_harvest[y,r], n_my_tot[y])
      
    }
    
    # --------------
    # HYPERPRIOR FOR THE CHENA AND SALCHA HARVEST #
    tau_h_riv_est[r] <- pow(1/sig_h_riv_est[r], 2)
    sig_h_riv_est[r] ~ dexp(1E-2)
    
  }
  
  # --------------
  # HARVEST IN THE MIDDLE YUKON #
  for (y in 1:n_years){
    log_h_my_est[y] ~ dnorm(log(H_my[y]), tau_h_my_est)
  }
  tau_h_my_est <- pow(1/sig_h_my_est, 2)
  sig_h_my_est ~ dexp(1E-2)
  
  
  ############################################################################################
  ############################ CALCULATING SOME USEFUL STATISTICS ############################
  ############################################################################################
  
  # ----------------  TOGGLE THE COMMENTS ACCORDING TO THE RS PROCESS ---------------- #
  
  
  for (r in 1:2){
    
    # ---------------- USE ONE OF THE FOLLOWING FOUR OPTIONS ---------------- # 
    
    # WITHOUT THE AR(1) TERM #
    alpha_prime[r] <- alpha[r]*exp(pow(sig_w[r], 2)/2)
  
    # # WITH THE AR(1) TERM #
    # alpha_prime[r] <- alpha[r]*exp(pow(sig_w[r], 2)/(2*(1-pow(phi[r], 2))))
  
    # # TIME VARYING PRODUCTIVITY WITHOUT THE AR(1) TERM #
    # for (y in 1:n_years){
    #  alpha_prime[y,r] <- alpha[y,r]*exp(pow(sig_w[r], 2)/2)
    # }
    
    # # TIME VARYING PRODUCTIVITY WITH THE AR(1) TERM #
    # for (y in 1:n_years){
    #   alpha_prime[y,r] <- alpha[y,r]*exp(pow(sig_w[r], 2)/(2*(1-pow(phi[r], 2))))
    # }
    
    # ---------------- USE ONE OF THE FOLLOWING THREE OPTIONS ---------------- #
    
    # FOR THE RICKER RS RELATIONSHIP WITHOUT TIME VARYING PRODUCTIVITY #
    S_msy[r] <- log(alpha_prime[r])/beta[r]*(0.5-0.07*log(alpha_prime[r]))
    S_max[r] <- 1/beta[r]
    S_eq[r] <- log(alpha_prime[r])/beta[r]
    U_msy[r] <- log(alpha_prime[r])*(0.5-0.07*log(alpha_prime[r]))
    
    # RICKER RS RELATIONSHIP WITH TIME VARYING PRODUCTIVITY PARAMETER #
    # for (y in 1:n_years){
    #   S_msy[y,r] <- log(alpha_prime[y,r])/beta[r]*(0.5-0.07*log(alpha_prime[y,r]))
    #   S_max[y,r] <- 1/beta[r]
    #   S_eq[y,r] <- log(alpha_prime[y,r])/beta[r]
    #   U_msy[y,r] <- log(alpha_prime[y,r])*(0.5-0.07*log(alpha_prime[y,r]))
    # }
    
    # # FOR THE BEVERTON-HOLT RS RELATIONSHIP #
    # S_msy[r] <- (sqrt(alpha_prime[r])-1)/beta[r]
    # S_max[r] <- 1/beta[r]
    # S_eq[r] <- (alpha_prime[r]-1)/beta[r]
    # U_msy[r] <- 1-sqrt(alpha_prime[r])/alpha_prime[r]
  
  }
  
}"
  
  fmod = "D:/Jordy/ChenaSalchaSpawnerRecruit/R/JAGS/SR3.R" 
  writeLines(mod,con=fmod)
  
  
  # RUN JAGS MODEL #
  library(rjags)
  jags_model = jags.model(fmod,
                          data = data,
                          n.chains = 4,
                          n.adapt = 20000)
  
  samples = coda.samples(jags_model,
                         variable.names = c("pi"),                              
                         n.burnin = 25000,                              
                         n.iter = 50000,                              
                         thin = 50)
  
  dic <- dic.samples(jags_model, 
                     n.burnin = 25000,                              
                     n.iter = 50000,
                     thin = 50)
  
  
  
  
  
  # CONVERGENCE DIAGNOSTICS #
  MCMCvis::MCMCtrace(samples)
  # GET SOME GELMAN RUBIN BROOKS PLOTS
  summary(samples)
  
  
  
  
  
  
  
  
  
  