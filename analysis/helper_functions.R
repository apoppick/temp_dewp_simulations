# This file contains the functions needed to estimate the statistical models and produce the observation-based simulations described in the paper.

get_GMTSmooth <- function(TGlobal, start.year = 1920, end.year = 2100, ref.years = 1990:2005, 
                          span = 0.05){
  # Estimates the forced GMT trend using a lowess smoother. 
  # Inputs:
  # TGlobal is the ensemble mean annual GMT from start.year to end.year
  # ref.years defines the period for calculating anomalies (i.e., GMTA averaged over ref.years is zero)
  # span is the span used in the lowess smoother
  #
  # Outputs:
  # GMT_df is a data frame consisting of the raw GMT (TGlobal), the smoothed estimated GMTA, and the corresponding years.
  
  nYears <- end.year - start.year + 1
  GMT_df <- data.frame(GMT_years = start.year:end.year,
                       GMT = rowMeans(TGlobal))
  loessGMT <- loess(GMT ~ GMT_years, data = GMT_df,
                    span = span, control = loess.control(surface = "direct"))
  GMT_df$GMTSmoother_year <- predict(loessGMT)
  GMT_df$GMTSmoother_year <- GMT_df$GMTSmoother_year - 
    mean(GMT_df$GMTSmoother_year[GMT_df$GMT_years %in% ref.years])
  
  return(GMT_df)
}

make_combinedData <- function(TREFHT, DEWPREFHT, GMT_df, dates, 
                              season_ind, num_harmonics){
  # Creates a tidy dataset for estimating the models used in the paper
  # Inputs:
  # TREFHT is a matrix consisting of the temperature values for each run at a location (nDates x nRuns).
  # DEWPREFHT is a matrix consisting of the dew point values for each run at a locatioj(nDates x nRuns).
  # GMT_df  is the output of get_GMTSmooth()
  # dates are the dates corresponding to TREFHT and DEWPREFHT
  # season_ind defines which days of the years are included in dates
  # num_harmonics defines how many harmonics are going to be used in the statistical models
  #
  # Outputs:
  # local_df is a data frame with columns correspond to temperature, dew point, dew point depression, smoothed GMTA, run, date, and the harmonics corresponding to the date. We also adjust any dew point and dew point depression values where the dew point depression is less than zero to the smallest positive value.
  
  nSeason <- sum(season_ind)
  nYears <- length(unique(year(dates)))
  nRuns <- dim(as.matrix(TREFHT))[2]
  
  X_c <- (cos((2*pi/365) * outer(0:(365 - 1), 1:num_harmonics)))[season_ind, ]
  X_s <- (sin((2*pi/365) * outer(0:(365 - 1), 1:num_harmonics)))[season_ind, ]
  
  GMTSmoother_subs <- GMT_df$GMTSmoother_year[GMT_df$GMT_years %in% unique(year(dates))]
  
  local_df <- data.frame(TREFHT = c(TREFHT),
                         DEWP = c(DEWPREFHT),
                         GMTSmoother = rep(rep(GMTSmoother_subs, rep(nSeason, nYears)),
                                           nRuns),
                         X_c = repmat(X_c, nYears*nRuns, 1),
                         X_s = repmat(X_s, nYears*nRuns, 1),
                         date = rep(dates, nRuns),
                         run = rep(1:nRuns, rep(nYears*nSeason, nRuns)))
  
  ## Using data.table
  # local_df <- data.frame(date = dates,
  #                        TREFHT = TREFHT,
  #                        DEWPREFHT = DEWPREFHT,
  #                        X_c = repmat(X_c, nYears, 1),
  #                        X_s = repmat(X_s, nYears, 1),
  #                       GMTSmoother = rep(GMTSmoother_year, rep(nSeason, nYears)))
  # 
  # harmonic_cols <- which(grepl("X_c", colnames(local_df), fixed = TRUE) |
  #                       grepl("X_s", colnames(local_df), fixed = TRUE))
  # 
  # 
  # local_df <- melt.data.table(as.data.table(local_df), id.vars = c(1, harmonic_cols, 86),
  #                                measure.vars = 2:(harmonic_cols[1]-1))
  # local_df[, c("variable", "run") := tstrsplit(variable, ".", fixed=TRUE)]
  # local_df$run <- as.numeric(local_df$run)
  # local_df <- dcast.data.table(local_df, ... ~ variable)
  # setorder(local_df, run, date)
  
  ## Using tidyr
  # local_df <- local_df %>%
  #               gather(key, value, -date, -GMTSmoother, -harmonic_cols) %>%
  #               extract(key, c("var", "run"), "(.+)\\.(.+)") %>%
  #               spread(var, value) %>% 
  #               mutate(run = as.numeric(run)) %>%
  #               arrange(run, date)
  
  local_df$DPS <- local_df$TREFHT - local_df$DEWP
  
  bad_dewp <- (local_df$DPS <= 0)
  min_pos_DPS <- min(local_df$DPS[!bad_dewp])
  local_df$DEWP[bad_dewp] <- local_df$TREFHT - min_pos_DPS
  local_df$DPS[bad_dewp] <- min_pos_DPS
  
  return(local_df)
}

repmat <- function(X,m,n){
  # R version of Matlab's repmat() function. X is a matrix; m is the number of row repetitions and n is the number of column repetitiions.
  
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow = TRUE)
}

estimate_meanTemp <- function(combinedData, I.scores  = TRUE, boot = FALSE){
  # estimates the mean temperature model from the paper
  # inputs:
  # combined Data is the output from make_combinedData()
  # I.scores is a logical indicating whether the diagnostic scores I_0 and I_1 should be calculated
  # boot is a logical indicating whether to produce bootstrap coefficient estimates. This is based on blocking by run, so only works with the CESM-LE data (not observations)
  #
  # outputs:
  # meanTemp_model is an lm object
  # I0 and I1 are the diagnostic scores
  # bootCoefs is a matrix of bootstrap coefficient estimates
  
  # this assumes 2 harmonics. It would be better to re-write this using the number of harmonic terms available in combinedData.
  meanTemp_model <- lm(TREFHT ~  GMTSmoother*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                       data = combinedData)
  if(!I.scores){
    return(list(model = meanTemp_model))
  } else {
    nRuns <- length(unique(combinedData$run))
    TREFHT_mat <- matrix(combinedData$TREFHT, ncol = nRuns, byrow = FALSE)
    ens_mean <- rowMeans(TREFHT_mat)
    time_mean <- colMeans(TREFHT_mat)
    nDays_transient <- dim(TREFHT_mat)[1]
    
    I1 <- sum((resid(meanTemp_model)[!is.na(combinedData$run)])^2) / 
      (nRuns / (nRuns - 1) * sum(sweep(TREFHT_mat, 1, ens_mean, FUN = "-")^2))
    I0 <- (nDays_transient / (nDays_transient - 1) * 
             sum(sweep(TREFHT_mat, 2, time_mean, FUN = "-")^2)) / 
      (nRuns / (nRuns - 1) * sum(sweep(TREFHT_mat, 1, ens_mean, FUN = "-")^2))
    
    if(boot){
      B <- 1e3
      bootCoefs <- matrix(NA, nrow = length(coef(meanTemp_model)), ncol = B)
      combinedData$TempAnom <- resid(meanTemp_model)
      for(b in 1:B){
        bootRuns <- sample(1:nRuns, replace = TRUE)
        bootTempAnom <- numeric(length(combinedData$TempAnom))
        for(r in 1:nRuns){
          bootTempAnom[combinedData$run == r] <- combinedData$TempAnom[combinedData$run == bootRuns[r]]
        }
        bootTemp <- combinedData$TREFHT + bootTempAnom
        # this assumes 2 harmonics. It would be better to re-write this using the number of harmonic terms available in combinedData.
        boot_meanModel <- lm(bootTemp ~  GMTSmoother*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                             data = combinedData)
        bootCoefs[,b] <-coef(boot_meanModel)
      }
    } else {
      bootCoefs = NULL
    }
    
    return(list(model = meanTemp_model, I1 = I1, I0 = I0, bootCoefs = bootCoefs))
  }
}

estimate_logRho <- function(combinedData){
  # estimates the variability change model from the paper. This function is a wrapper for estimate_logRho_old()
  # inputs:
  # combinedData is the output from make_combinedData()
  #
  # outputs:
  # delta_hat_smoothed is the estimated function delta(omega) from the model (change in log spectrum per degree increase in GMT)
  # V_deltaSmooth is the estimated variance of delta_hat_smoothed(omega)
  # loga0_hat is the estimate of the spectrum in the first year of the data record (not smoothed)
  # deviance is the deviance score for diagnostics
  
  nRuns <- length(unique(combinedData$run))
  
  local_temp <- matrix(combinedData$TREFHT, ncol = nRuns, byrow = FALSE)
  
  GMT_year <- unique(combinedData$GMTSmoother) #make this nicer
  GMT_year <- GMT_year - GMT_year[1]
  
  seasonInd <- (month(seq.Date(ymd("1990-01-01"), ymd("1990-12-31"), 1)) %in% 
                  unique(month(combinedData$date))) # make this better
  
  estimate_logRho_old(local_temp, local_temp_PI = NULL, GMT_year, seasonInd, 
                  smoother_bandwidth = NA)
}

estimate_logRho_old <- function(local_temp, local_temp_PI = NULL, GMT_year, seasonInd, smoother_bandwidth = NA){
  # estimates the model logRho(omega, y) = delta(omega) * GMTA_Year(y), where logRho(omega, y) is the ratio of 
  # spectra in year y vs. year zero. If local_temp_PI is NULL, year zero is interpretted as the first year of 
  # local_temp.
  #
  #INPUTS:
  # local_temp: an nDays x nRuns matrix of temperatures from one gridcell
  # local_temp_PI: an nPI x 1 vector of PI temperatures from the same gridcell
  # GMT_year: an nDays x 1 vector of global mean temperature anomalies (relative to year zero). If local_temp_PI is                 NULL, the first value should be zero; otherwise, the difference vs. the PI mean temperature
  # seasonInd: a 365 x 1 logical vector indicating the days in the year being supplied in local_temp and local_temp_PI
  # smoother_bandwidth: bandwidth for smoothed estimate of delta. If NA, this is chosen automatically.
  
  nSeason <- sum(seasonInd)
  nDays <- dim(local_temp)[1]
  nRuns <- dim(local_temp)[2]
  nYears <- nDays / nSeason
  
  Tanom <- sweep(local_temp, 1, rowMeans(local_temp))
  Tanom_reshape <- array(Tanom, dim = c(nSeason, nYears, nRuns))
  I_reshape <- Mod(apply(Tanom_reshape, c(2,3), fft)/sqrt(nSeason))^2
  I_bar <- rowMeans(I_reshape, dims = 2)
  #I_bar[1,] <- I_bar[2,] #think about this
  
  if(!is.null(local_temp_PI)){
    nDays_PI <- length(local_temp_PI)
    nYears_PI <- nDays_PI / nSeason
    fullDateInd <- 1:(365*nYears_PI)
    SeasonDateInd <- fullDateInd[rep(seasonInd, nYears_PI)]
    num_harmonics <- 2 # this should be sufficient for summer, but might need to be bigger if doing the whole year
    X_c <- cos((2*pi/365) * outer(SeasonDateInd, 1:num_harmonics))
    X_s <- sin((2*pi/365) * outer(SeasonDateInd, 1:num_harmonics))
    mean_PI_lm <- lm(local_temp_PI ~ X_c + X_s)
    Tanom_PI <- resid(mean_PI_lm)
    
    Tanom_PI_reshape <- array(Tanom_PI, dim = c(nSeason, nYears_PI))
    I_PI_reshape <- Mod(mvfft(Tanom_PI_reshape)/sqrt(nSeason))^2
    I_bar_PI <- rowMeans(I_PI_reshape)
  } else { #is.null(local_temp_PI)
    nYears_PI <- 0
    I_bar_PI <- 0
  }
  
  if(!is.null(local_temp_PI)){
    loga0_hat <- log(I_bar_PI)
  } else {
    loga0_hat <- log(nRuns/(nRuns - 1) * I_bar[,1])
  }
  delta_hat <- rowSums(sweep(log(I_bar), 1, loga0_hat)) / sum(GMT_year)
  stopCond <- 100
  while(stopCond > 1e-8){
    GDelta <- outer(delta_hat, GMT_year)
    loga0_hat_new <- log((nYears_PI*I_bar_PI  + nRuns* rowSums(I_bar*exp(-GDelta)))/(nYears_PI + nYears*(nRuns-1)))
    
    a_omega <- exp(sweep(GDelta, 1, loga0_hat_new, FUN = "+"))
    g_term <- (nRuns - 1) * sum(GMT_year) - nRuns * rowSums(sweep(I_bar/a_omega, 2, GMT_year, FUN = "*"))
    h_term <- nRuns * rowSums(sweep(I_bar/a_omega, 2, GMT_year^2, FUN = "*"))
    delta_hat_new <- delta_hat - g_term / h_term
    stopCond <- sum((g_term / h_term)^2)
    loga0_hat <- loga0_hat_new
    delta_hat <- delta_hat_new
  }
  
  delta_hat_smoothed <- autoSmoothBP(delta_hat, smoother_bandwidth) # check this
  delta_hat_smoothed$delta_smooth[1] <- delta_hat[1]
  
  V_deltaRough <- (2*nYears*(nRuns-1) + nYears_PI) / ((nRuns-1)*sum(GMT_year^2)*(nYears*(nRuns-1) + nYears_PI) - 
                                                        ((nRuns-1)*sum(GMT_year))^2)
  V_deltaSmooth <- V_deltaRough * delta_hat_smoothed$V_multiplier
  V_deltaSmooth[1] <- V_deltaRough[1]
  
  a_omega <- exp(sweep(outer(c(delta_hat_smoothed$delta_smooth, delta_hat_smoothed$delta_smooth[46:2]),
                             GMT_year), 1, loga0_hat, FUN = "+"))
  
  # Deviance does not include PI terms here because they seem irrelevant to the question of whether the model captures the changes..
  a_omega_rough <- exp(sweep(outer(delta_hat,GMT_year), 1, loga0_hat, FUN = "+"))
  dev <-   sum((nRuns - 1) * log(a_omega_rough) + nRuns * I_bar/a_omega_rough) - 
    sum((nRuns - 1) * log(I_bar) + nRuns)
  
  
  return(list(delta_hat_smoothed = delta_hat_smoothed$delta_smooth, V_deltaSmooth = V_deltaSmooth,
              loga0_hat = loga0_hat, deviance = dev))
}

fftshift <- function(A){ #make this better
  # version of matlab's fftshift applied to the rows of the matrix A
  nRows <- dim(A)[1]
  nRows_half <- ceiling(nRows/2)
  
  Ashift <- rbind(A[((nRows_half+1):nRows), , drop = FALSE],
                  A[(1:nRows_half), , drop = FALSE])
  
  return(Ashift)
}

autoSmoothBP <- function(delta_rough, M = NA){
  # Smooths the columns of delta_rough using an epanechnikov kernel, with the bandwidth M chosen
  # via the function RCV(), same bandwidth for each column. The secondary output, V_multiplier, is 
  # the number such that the variance of the smoothed estimator is equal to V_multiplier*Var(delta_rough)
  
  delta_rough <- as.matrix(delta_rough)
  N <- dim(delta_rough)[1]
  Ks <- -floor(N/2):(N-1)
  Js <- 0:(N-1)
  
  extTheta <- rbind(fftshift(delta_rough), delta_rough[(ceiling(N/2)+1):N, , drop = FALSE])
  
  S0 <- as.matrix(cumsum(extTheta))
  S1 <- as.matrix(cumsum(sweep(extTheta, 1 , Ks, FUN = "*")))
  S2 <- as.matrix(cumsum(sweep(extTheta, 1, Ks^2, FUN = "*")))
  
  if(is.na(M)){
    RCVFun <- function(x) RCV(delta_rough,S0,S1,S2,x);
    M <- round(optimize(RCVFun, c(1, round(N/4) - 1))$minimum)
  }
  
  if(M == round(N/4)-1){
    warning('Final estimate on boundary')
  }
  
  l <- ((-floor(N/2)):(N-1))
  j <- (0:(N/2))
  
  bigWeightMat <- (1-(outer(l, j, FUN = "-")/(M+1))^2) *
    (abs(outer(l, j, FUN = "-"))<= M) * (3*M + 3)/(4*M^2 + 8*M + 3)
  #each column j of this matrix sums of 1, 
  #with delta_smooth[j] = sum_l W[l,j] * delta_rough[l]
  
  realWeightMat <- bigWeightMat[l >= 0 & l <= N/2, ]
  numRealWeights <- N/2+1
  #add periodic weights to the left:
  realWeightMat[2:(numRealWeights-1),] = realWeightMat[2:(numRealWeights-1),] + 
    bigWeightMat[seq(N/2, 2, by = -1), 1:(N/2+1)]
  
  #add periodic weights to the right:                              
  realWeightMat[2:(numRealWeights-1),] = realWeightMat[2:(numRealWeights-1),] + 
    bigWeightMat[seq(N+N/2, N+2, by = -1), 1:(N/2+1)]; 
  
  V_multiplier <- colSums(realWeightMat^2) #(3/5)*(4*M^2 + 8*M + 5)/(4*M^3+12*M^2+11*M+3)
  delta_smooth <- SmoothBP(S0,S1,S2,M,N)
  
  return(list(delta_smooth = delta_smooth[1:(floor(N/2)+1)],
              V_multiplier = V_multiplier))
}

RCV <- function(lograt,S0,S1,S2,M){
  # Calculates crossvalidation mean squared error of estimates delta(omega) function for variability change model (see notes for derivation of the score).
  # lograt is the original unsmoothed estimate. 
  #S0, S1, and S2 are cumulative sums (see code for autoSmoothBP)
  # M is the bandwidth of the kernel
  
  M <- round(M)
  N <- dim(lograt)[1]
  nLocs <- dim(lograt)[2]
  
  weight0 <- (3*M + 3)/(4*M^2 + 8*M + 3)
  sIrat <- SmoothBP(S0,S1,S2,M,N)
  
  colSums(as.matrix(rowSums((sIrat - lograt)^2)/(1-weight0)^2))/(N*nLocs)
}

SmoothBP <- function(S0,S1,S2,M,N){
  # Calculates kernel smoothed estimate of delta(omega) using an Epanechnikov kernel. This code is written so that the smoother can be calculated essentially instantly after computing three comulative sums once. This is helpful when trying to optimize the bandwidth of the kernel.
  #S0, S1, and S2 are cumulative sums of the unsmoothed estimate and products involving the kernel (see code for autoSmoothBP)
  # M is the bandwidth of the kernel
  # N is the length of the time series
  
  
  M <- round(M)
  Js <- 0:(N-1)
  ind1 <- (1:(floor(N/2)+1)) + M + ceiling(N/2)
  ind2 <- (ceiling(N/2):N) - M
  
  sItemp <- (S0[ind1, ] - S0[ind2, ]) - 
    (
      (S2[ind1,] - S2[ind2,]) - 
        sweep(as.matrix(S1[ind1,]-S1[ind2,]), 1, 2*Js[1:length(ind1)], FUN = "*") +
        sweep(as.matrix(S0[ind1,]-S0[ind2,]), 1, Js[1:length(ind1)]^2, FUN = "*")
    )/(M+1)^2
  
  sItemp <- rbind(as.matrix(sItemp[1:(floor(N/2)+1),]),
                  as.matrix(sItemp[ceiling(N/2):2,])) * (3*M+3)/(4*M^2+8*M+3)
  return(sItemp)
}

estimate_quantModel <- function(combinedData, obs = FALSE, pearson_stat = TRUE,
                                taus = c(0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 
                                         0.75, 0.9, 0.95, 0.99, 0.995)){
  # estimates the quantile regression model used in the paper
  #
  # inputs:
  # combinedData is the output of make_combinedData()
  # obs is a logical indicating whether combinedData is from observation (TRUE) or CESM-LE (FALSE)
  # pearson_stat is a logical indicating whether the Pearson statistic described in the paper supplement should be calculated.
  # taus are the quantile levels at which to estimate the model
  
  N <- dim(combinedData)[1]
  first.method <- ifelse(N > 1e4, "pfn", "br")
  
  err <- 0
  if(!obs){
    quantModel <- tryCatch({
      rq(log(DPS) ~ GMTSmoother +  
           (ns(TempAnom, df = 10) + 
              (X_c.1 + X_c.2 + X_s.1 + X_s.2)), 
         data = combinedData, 
         tau = taus, method = first.method)
    }, error = function(e) {
      tryCatch({
        err[1] <<- 1
        return(rq(log(DPS) ~ GMTSmoother +  
                    (ns(TempAnom, df = 10) + 
                       (X_c.1 + X_c.2 + X_s.1 + X_s.2)), 
                  data = combinedData,
                  tau = taus))
      }, error = function(e2){
        err[1] <<- 2
        return(NULL)
      }
      )
    }
    )
    
    if(pearson_stat){
      ## Calculate Pearson statistic for pooled data
      e <- apply(quantModel$resid, 
                 MARGIN = 1, FUN = sort, decreasing = TRUE)
      
      # get upper and lower bounds of quantile level for each observation
      e_minus <- ifelse(e < 0, -e, NA)
      e_plus <- ifelse(e > 0, e, NA)
      tau_U_ind <- apply(e_minus, 2, which.min)
      tau_U_ind <- unlist(lapply(tau_U_ind, function(x) ifelse(length(x) == 1, x, NA)))
      tau_U_ind[is.na(tau_U_ind)] <- length(taus)
      
      combinedData$tau_U_ind <- as.factor(tau_U_ind)
      combinedData$period <- ifelse(year(combinedData$date) <= 2005, "historical",
                                    ifelse(year(combinedData$date) <= 2035, "f1", "f2"))
      combinedData$TempAnomGroup <- cut(combinedData$TempAnom, 
                                        breaks = quantile(combinedData$TempAnom, c(0, 1/4, 1/2, 3/4, 1)))
      count_table <- table(combinedData$tau_U_ind, combinedData$period, month(combinedData$date), 
                           combinedData$TempAnomGroup)
      theor_props <- c(taus[1], diff(taus))
      theor_counts <- array(dim = dim(count_table))
      for(p in 1:3){
        for(m in 1:3){
          for(t in 1:4){
            theor_counts[,p,m,t] <- theor_props * sum(count_table[,p,m,t])
          }
        }
      }
      
      pearson_stat <- sum((count_table - theor_counts)^2/theor_counts)
    } else {
      pearson_stat <- NULL
    }
  } else { #if observations, add jitter and don't include GMT in model
    combinedData$DPS <- combinedData$DPS + runif(dim(combinedData)[1], min = -0.05, max = 0.05) + 
      runif(dim(combinedData)[1], min = -0.05, max = 0.05) * 5/9 - 
      runif(dim(combinedData)[1], min = -0.05, max = 0.05) - 
      runif(dim(combinedData)[1], min = -0.05, max = 0.05) * 5/9 
    combinedData$DPS[combinedData$DPS <= 0] <- min(combinedData$DPS[combinedData$DPS > 0])
    quantModel <- rq(log(DPS) ~ (ns(TempAnom, df = 10) + 
                                   (X_c.1 + X_c.2 + X_s.1 + X_s.2)), 
                     data = combinedData,
                     tau = taus)
    pearson_stat <- NULL
  }

  
  return(list(model = quantModel, err = err, pearson_stat = pearson_stat))
}

calculate_riskProbs <- function(combinedData, quantModel, meanTemp){
  # calculates changes in the risk of high and low humidity events given historically high temperature in CESM-LE
  #
  # inputs:
  # combinedData is the output of make_combinedData()
  # quantModel is the output of estimate_quantModel()
  # meanTemp is the output of estimate_meanTemp()
  #
  # output: a data frame consisting of
  # TempRisk_future_CESM: the 2071-2081 quantile level associated with the 1990-2005 95th percentile of temperature.
  # DewpRisk_future_high: the 2071-2080 quantile level associated with the 1990-2005 95th percentile of dew point, given the 1990-2005 95th percentile of temperature
  # DewpRisk_future_low: the 2071-2080 quantile level associated with the 1990-2005 5th percentile of dew point, given the 1990-2005 95th percentile of temperature
  
  num_harmonics <- sum(grepl("X_c", names(combinedData)))
  month_starts <- c(152, 182, 213) #6/1, 7/1, and 8/1
  
  X_c <- (cos((2*pi/365) * outer(month_starts-1, 1:num_harmonics)))
  X_s <- (sin((2*pi/365) * outer(month_starts-1, 1:num_harmonics)))
  harmonics_df <- data.frame(months = 6:8, X_c = X_c, X_s = X_s)
  
  simulDates <- unique(combinedData$date[year(combinedData$date) >= 2071])
  #GMT_year <- tapply(GMT_JJA$GMTSmoother, year(GMT_JJA$date), mean)
  
  taus <- quantModel$model$tau
  
  TempRisk_future <- numeric(3)
  DewpRisk_future_low <- numeric(3)
  DewpRisk_future_high<- numeric(3)
  for(m in 6:8){
    tempThresh <- quantile(combinedData$TREFHT[month(combinedData$date) == m & 
                                                    year(combinedData$date) %in% 1990:2005], 
                                0.95)
    
    newData_hist <- harmonics_df %>% filter(months == m) %>% mutate(GMTSmoother = 0)
    tempAnomThresh <- tempThresh - predict(meanTemp$model, newData_hist)
    newData_hist <- newData_hist %>% mutate(TempAnom = tempAnomThresh)
    lDPD_preds_hist <- predict(quantModel$model, newdata = newData_hist)
    dewpThresh_low <- tempThresh  - exp(lDPD_preds_hist[taus == 0.95])
    dewpThresh_high <- tempThresh  - exp(lDPD_preds_hist[taus == 0.05])
    
    GMT_future <- mean(combinedData$GMTSmoother[year(combinedData$date) >= 2071])
    newData_future <- harmonics_df %>% filter(months == m) %>% mutate(GMTSmoother = GMT_future)
    tempAnom_future <- tempThresh - predict(meanTemp$model, newData_future) # check this
    newData_future <- newData_future %>% mutate(TempAnom = tempAnom_future)
    lDPD_preds_future <- predict(quantModel$model, newdata = newData_future)
    dewp_preds_future <- tempThresh  - exp(lDPD_preds_future)
    
    TempRisk_future[m-5] <- mean(combinedData$TREFHT[year(combinedData$date) >= 2071 &
                                                       month(combinedData$date) == m] > tempThresh)
    DewpRisk_future_high[m-5] <- 1 - approx(dewp_preds_future, taus, dewpThresh_high, rule = 2)$y
    DewpRisk_future_low[m-5]  <- 1 - approx(dewp_preds_future, taus, dewpThresh_low, rule = 2)$y
  }
  riskProbs <- data.frame(month = 6:8, 
                          TempRisk_future_CESM = TempRisk_future,
                          DewpRisk_future_high_CESM = DewpRisk_future_high,
                          DewpRisk_future_low_CESM = DewpRisk_future_low)
  
  return(riskProbs)
}

simulTemp <- function(obs_data, obs_meanModel, GMT, CESM_meanCoefs, delta, obsYears, simulYears){
  # Simulate future temperatures according to the method described in the paper. This is a wrapper for simulTemp_old().
  #
  # inputs:
  # obs_data are observations output from make_combinedData()
  # obs_meanModel is output of estimate_meanModel(obs_data)
  # GMT is output of get_GMTSmooth
  # CESM_meanCoefs are the mean model coefficients for CESM
  # delta is the estimated delta(omega) function from estimate_logRho()
  # obsYears is a vector of years defining the subset of obs_data to use for the future simulation
  # simulYears defines the years for the future simulation (should be same length as obsYears)
  #
  # outputs:
  # simul is the temperature simulation
  # anom is simul - the estimated forced mean
  
  nSeason <- dim(obs_data[year(obs_data$date) == first(year(obs_data$date)),])[1]
  histData <- subset(obs_data, year(date) %in% obsYears)
  futureData <- histData
  futureData$GMTSmoother <- rep(GMT$GMTSmoother[GMT$GMT_years %in% simulYears], 
                                rep(nSeason, length(simulYears)))
  
  X_hist <- model.matrix(obs_meanModel$model, data = histData)
  X_future <- model.matrix(obs_meanModel$model, data = futureData)
  
  Delta <- X_future %*% CESM_meanCoefs - X_hist %*% CESM_meanCoefs
  histMean <- X_hist %*% coef(obs_meanModel$model)
  
  GMT_year_change <- GMT$GMTSmoother[GMT$GMT_years %in% simulYears] - 
                    GMT$GMTSmoother[GMT$GMT_years %in% obsYears]
  logRho <- outer(c(delta,delta[46:2]), GMT_year_change)
  
  simul <- simulTemp_old(histData$TREFHT, obsMeanTemp = histMean,
                                 meanTempChange = Delta,
                                 logRho = logRho,
                                 nSeason = 92)
  anom <- simul - (histMean  + Delta)
  return(list(simul = simul, anom = anom))
}

simulTemp_old <- function(obsTemp, obsMeanTemp, meanTempChange, logRho, nSeason){
  # Simulate future temperatures according to the method described in the paper. This is an older version of the function that is more computationally efficient but accepts inputs in a form that is inconsistent with the rest of my code.
  #
  # inputs:
  # obsTemp is a vector of observed temperatures
  # obsMeanTemp is a vector of estimated mean temperatures
  # meanTempChange is a vector of the change in mean temperature from CESM-LE
  # logRho is a matrix containing the values delta(omega)*Delta(y) from the variability change model
  # nSeason is the number of days in the season (92 for JJA).
  #
  # outputs a vector of simulated temperatures
  
  nDays <- length(obsTemp)
  nYears <- nDays/nSeason
  obsAnom <- obsTemp - obsMeanTemp
  obsAnom_reshape <- array(obsAnom, dim = c(nSeason, nYears))
  F_reshape <- mvfft(obsAnom_reshape)
  newAnom_reshape <- Re(mvfft(exp(0.5*logRho) * F_reshape, inverse = TRUE)/nSeason)
  return(obsMeanTemp + meanTempChange + c(newAnom_reshape))
}

simulDEWP <- function(obs_data, quantModel_OBS, GMT, quantModel_CESM_coefs, simulTemp_fromObs,
                      obsYears, simulYears){
  # Simulate future dew points given a future temperature simulation
  # inputs:
  # obs_data is observed data output from make_combinedData()
  # quantModel_OBS is output from estimate_quantModel(obs_data)
  # quantModel_CESM_coefs are quantile regresison model coefficients estimated from CESM
  # simulTemp_fromObs is output from simulTemp()
  # obsYears are the years of obs_data to use for the simulation
  # simulYears are the years for the simulation
  #
  # output:
  # vector of simulated dew points
  
  nSeason <- dim(obs_data[year(obs_data$date) == first(year(obs_data$date)),])[1]
  
  GMT_coefs <- quantModel_CESM_coefs[rownames(quantModel_CESM_coefs) == "GMTSmoother", ]
  histData <- subset(obs_data, year(date) %in% obsYears)
  
  # calculate observed residuals and fitted values
  changeTerm_hist <- outer(histData$GMTSmoother, GMT_coefs)
  
  e_obs <- apply(quantModel_OBS$model$resid[year(obs_data$date) %in% obsYears,] - changeTerm_hist, 
                 MARGIN = 1, FUN = sort, decreasing = TRUE)
  
  
  yhat_obs <- apply(quantModel_OBS$model$fitted[year(obs_data$date) %in% obsYears,] + changeTerm_hist, 
                    MARGIN = 1, FUN = sort, decreasing = FALSE)
  
  nDPSimul <- dim(yhat_obs)[2]
  
  # calculate fitted quantiles for simulation period
  futureData <- subset(histData, year(date) %in% obsYears)
  futureData$TempAnom <- simulTemp_fromObs$anom
  
  GMT_year_change <- GMT$GMTSmoother[GMT$GMT_years %in% simulYears] - 
                          GMT$GMTSmoother[GMT$GMT_years %in% obsYears]
  changeTerm <- outer(rep(GMT_year_change, rep(nSeason, length(simulYears))), GMT_coefs)
  
  logDPS_hats <- predict(quantModel_OBS$model, newdata = futureData)
  logDPS_hats <- t(apply(logDPS_hats, 1, sort)) + changeTerm #+ t(changeTerm_hist)
  
  # get upper and lower bounds of quantile level for each observation
  taus <- quantModel_OBS$model$tau
  e_obs_minus <- ifelse(e_obs < 0, -e_obs, NA)
  e_obs_plus <- ifelse(e_obs > 0, e_obs, NA)
  tau_U_ind <- apply(e_obs_minus, 2, which.min)
  tau_U_ind <- unlist(lapply(tau_U_ind, function(x) ifelse(length(x) == 1, x, NA)))
  tau_U_ind[is.na(tau_U_ind)] <- length(taus)
  
  tau_L_ind <- apply(e_obs_plus, 2, which.min)
  tau_L_ind <- unlist(lapply(tau_L_ind, function(x) ifelse(length(x) == 1, x, NA)))
  tau_L_ind[is.na(tau_L_ind)] <- 1
  
  
  tau_star <- numeric(nDPSimul)
  logDPS_hat <- numeric(nDPSimul)
  DPSObs_forSimul <- histData$DPS
  for(k in 1:nDPSimul){
    if(tau_L_ind[k] == tau_U_ind[k]){
      tau_star[k] <- taus[tau_L_ind[k]]
      logDPS_hat[k] <- logDPS_hats[k,tau_L_ind[k]]
    } else {
      tau_star[k] <- taus[tau_L_ind[k]] + 
        (log(DPSObs_forSimul[k]) - yhat_obs[tau_L_ind[k],k]) *
        (taus[tau_U_ind[k]] - taus[tau_L_ind[k]]) /
        (yhat_obs[tau_U_ind[k],k] - yhat_obs[tau_L_ind[k],k])
      
      logDPS_hat[k] <- logDPS_hats[k,tau_L_ind[k]] + 
        (tau_star[k] - taus[tau_L_ind[k]]) *
        (logDPS_hats[k,tau_U_ind[k]] - logDPS_hats[k,tau_L_ind[k]]) /
        (taus[tau_U_ind[k]] - taus[tau_L_ind[k]])
    }
  }
  DEWP_hat <- simulTemp_fromObs$simul - exp(logDPS_hat)
  return(DEWP_hat)
}
  
calculate_riskProbs_forSimul <- function(obs_data, meanTempModel_OBS, quantModel_OBS, GMT, 
                                         quantModel_CESM_coefs, meanTempModel_CESM,
                                         simulTemp_fromObs){
  # this is the analogue of calculate_riskProbs() but for the observation-based simulated temperature and dew point values
  
  num_harmonics <- sum(grepl("X_c", names(obs_data)))
  month_starts <- c(152, 182, 213) #6/1, 7/1, and 8/1
  simulDates <- obs_data$date[1:length(simulTemp_fromObs$simul)] # just need the days here, not the years
  
  X_c <- (cos((2*pi/365) * outer(month_starts-1, 1:num_harmonics)))
  X_s <- (sin((2*pi/365) * outer(month_starts-1, 1:num_harmonics)))
  harmonics_df <- data.frame(months = 6:8, X_c = X_c, X_s = X_s)
  
  taus <- quantModel_OBS$model$tau
  GMT_coefs <- quantModel_CESM_coefs[rownames(quantModel_CESM_coefs) == "GMTSmoother", ]
  
  TempRisk_future <- numeric(3)
  DewpRisk_future_low <- numeric(3)
  DewpRisk_future_high<- numeric(3)
  
  for(m in 6:8){
    tempThresh <- quantile(obs_data$TREFHT[month(obs_data$date) == m & 
                                       year(obs_data$date) %in% 1990:2005], 0.95)
    obsNewData <- harmonics_df %>% filter(months == m) %>% mutate(GMTSmoother = 0)
    tempAnomThresh <- tempThresh - predict(meanTempModel_OBS$model, obsNewData)
    obsNewData <- obsNewData %>% mutate(TempAnom = tempAnomThresh, TREFHT = tempThresh)
    lDPD_preds_hist <- predict(quantModel_OBS$model, newdata = obsNewData) + 
      obsNewData$GMTSmoother * GMT_coefs
    dewpThresh_low <- tempThresh  - exp(lDPD_preds_hist[taus == 0.95])
    dewpThresh_high <- tempThresh  - exp(lDPD_preds_hist[taus == 0.05])
    
    simulNewData <- harmonics_df %>% filter(months == m) %>%
      mutate(GMTSmoother = mean(GMT$GMTSmoother_year[GMT$GMT_years %in% 2071:2080]),
             TREFHT = tempThresh)
    
    X_hist <- model.matrix(meanTempModel_OBS$model, data = obsNewData)
    X_future <- model.matrix(meanTempModel_OBS$model, data = simulNewData)
    Delta_forSimulObs <- X_future %*% meanTempModel_CESM - X_hist %*% meanTempModel_CESM
    
    tempAnom_future <- tempThresh - (predict(meanTempModel_OBS$model, obsNewData) + Delta_forSimulObs) 
    simulNewData <- simulNewData %>% mutate(TempAnom = tempAnom_future)
    lDPD_preds_future <- predict(quantModel_OBS$model, newdata = simulNewData) + 
        mean(simulNewData$GMTSmoother - obsNewData$GMTSmoother)*GMT_coefs
    
    dewp_preds_future <- tempThresh  - exp(lDPD_preds_future)
    
    TempRisk_future[m-5] <- mean(simulTemp_fromObs$simul[month(simulDates) == m] > tempThresh)
    DewpRisk_future_high[m-5] <- 1 - approx(dewp_preds_future, taus, dewpThresh_high, rule = 2)$y
    DewpRisk_future_low[m-5]  <- 1 - approx(dewp_preds_future, taus, dewpThresh_low, rule = 2)$y
  }
  
  riskProbs <- data.frame(month = 6:8, 
                          TempRisk_future_OBS = TempRisk_future,
                          DewpRisk_future_high_OBS = DewpRisk_future_high,
                          DewpRisk_future_low_OBS = DewpRisk_future_low)
  return(riskProbs)
}