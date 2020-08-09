# This code contains functions that create figures showing the behavior at individual locations.

meanDiagnostics <- function(combinedData, meanTemp, lon, lat){
  # plot (a) ensemble mean temperature and modeled mean temperature, and (b) difference between the two.
  
  # combinedData is the output of make_combinedData()
  # meanTemp is the output of estimate_meanTemp()
  # lon is the longitude of the location (for plot title and file name)
  # lat is the latitude of the location being plotted (for plot title and file name)
  
  ens_mean <- tapply(combinedData$TREFHT, combinedData$date, mean)
  JJADates <- combinedData$date[combinedData$run == 1]
  
  
  validation_df <- data.frame(ens_mean = tapply(combinedData$TREFHT, combinedData$date, mean), 
                              JJADates = combinedData$date[combinedData$run == 1],
                              mean_fitted =fitted(meanTemp$model)[combinedData$run == 1])
  
  plot_dir <- "/research/Temp_Dewpt/Results/figures/validation/"
  f <- paste0(plot_dir, paste0("MeanModel_Lon", -round(lon), 
                               "Lat", round(lat),".pdf"))
  pdf(f, width = 10, height = 7)
  par(mfrow = c(2, 1), mar = c(5,5,4,4))
  plot(ens_mean ~ JJADates, data = validation_df, pch = 16, col = 'gray', xlab = "Date", 
       ylab = expression(paste("Mean Temperature (", degree, "C)")))
  lines(mean_fitted ~ JJADates, data = validation_df, col = 'red', 
        subset = (year(JJADates_6hourly) %in% 1990:2005))
  lines(mean_fitted ~ JJADates, data = validation_df, col = 'red', 
        subset = (year(JJADates_6hourly) %in% 2026:2035))
  lines(mean_fitted ~ JJADates, data = validation_df, col = 'red', 
        subset = (year(JJADates_6hourly) %in% 2071:2080))
  legend("topleft", pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), col = c('gray', 'red'),
         legend = c("Ensemble mean", "Estimated by model"))
  
  plot(I(ens_mean - mean_fitted) ~ JJADates, data = validation_df, pch = 16, col = 'black',
       xlab = "Date", 
       ylab = "Ensemble mean - Modeled mean")
  abline(h=0, col = 'gray', lwd = 2, lty = 'dashed')
  
  mtext(text = paste0("Lon ", round(lon, 1), 
                      ", Lat ", round(lat, 1), 
                      " (I1 = ", round(meanTemp$I1,3),")"), outer = TRUE, line = -1.5, cex = 1)
  dev.off()
}


varDiagnostics <- function(combinedData, logRhoEst, lon, lat){
  # plot ensemble JJA average temperature standard deviation and modeled estimated
  
  # combinedData is the output of make_combinedData()
  # logRhoEst is the output of estimate_logRho()
  # lon is the longitude of the location (for plot title and file name)
  # lat is the latitude of the location being plotted (for plot title and file name)
  
  anAvg <- tapply(combinedData$TREFHT, list(year(combinedData$date), combinedData$run), mean)
  var_anAvg <- apply(anAvg, 1, var)
  GMT_year <- unique(combinedData$GMTSmoother)
  GMT_year_change <- GMT_year - GMT_year[1]
  var_anAvg_hat <-  exp(logRhoEst$loga0_hat[1] + logRhoEst$delta_hat_smoothed[1]*GMT_year_change)/92
  years <- unique(year(combinedData$date))
  
  plot_dir <- "/research/Temp_Dewpt/Results/figures/validation/"
  f <- paste0(plot_dir, paste0("VarModel_Lon", -round(lon), 
                               "Lat", round(lat),".pdf"))
  
  pdf(f)
  par(mfrow = c(1,1))
  plot(sqrt(var_anAvg) ~ years, 
       ylab = expression(paste("SD of JJA Average (", degree, "C)")), 
       xlab = "Year", main = "", pch = 16, col = 'gray')
  lines(sqrt(var_anAvg_hat) ~ years, col = 'red', 
        subset = (years %in% 1990:2005), lwd = 2)
  lines(sqrt(var_anAvg_hat) ~ years, col = 'red', 
        subset = (years %in% 2026:2035), lwd = 2)
  lines(sqrt(var_anAvg_hat) ~ years, col = 'red', 
        subset = (years %in% 2071:2080), lwd = 2)
  legend("bottomright", pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), col = c("gray", "red"),
         legend = c("Ensemble SD", "Estimated by model"))
  mtext(text = paste0("Lon ", round(lon, 1), 
                      ", Lat ", round(lat, 1), 
                      " (Deviance = ", round(logRhoEst$deviance, 1),")"), outer = FALSE, line = 0.5, cex = 1.5)
  dev.off()
}

quantModelDiagnostics <- function(combinedData, quantModel, meanTemp, lon, lat){
  # plot July dew point vs. temperature in 1990-2005 and 2071-2080, along with mid-July estimate of quantile functions when GMT is the average value in those two time periods, and empirical estimates of the quantiles based on pooling by time period and temperature deviation value.
  
  # combinedData is the output of make_combinedData()
  # quantModel is the output of estimate_quantModel()
  # meanTemp is the output of estimate_meanTemp()
  # lon is the longitude of the location (for plot title and file name)
  # lat is the latitude of the location being plotted (for plot title and file name)
  
  
  numNewTemps <- 100
  newTempAnom1 <- seq(quantile(combinedData$TempAnom, 0.005), quantile(combinedData$TempAnom, 0.995), 
                      len = numNewTemps)
  avgGMT1 <- mean(combinedData$GMTSmoother[year(combinedData$date) %in%  1990:2005])
  newData1 <- data.frame(TempAnom = newTempAnom1, 
                         GMTSmoother = rep(avgGMT1, numNewTemps),
                         X_c.1 = rep(combinedData$X_c.1[31+15], numNewTemps),
                         X_c.2 = rep(combinedData$X_c.2[31+15], numNewTemps),
                         X_s.1 = rep(combinedData$X_s.1[31+15], numNewTemps),
                         X_s.2 = rep(combinedData$X_s.2[31+15], numNewTemps))
  newX1 <- model.matrix(~ GMTSmoother*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                        data = newData1)
  meanTemps1 <- newX1[1,] %*% meanTemp$model$coef
  newTREFHT1 <- newTempAnom1 + meanTemps1
  newData1$TREFHT <- newTREFHT1
  quantModelPreds1 <- predict(quantModel$model, newdata = newData1)
  quantModelPreds1 <-  t(apply(quantModelPreds1, 1, sort))
  
  newTempAnom2 <- newTempAnom1
  avgGMT2 <- mean(combinedData$GMTSmoother[year(combinedData$date) %in%  2071:2080])
  newData2 <- data.frame(TempAnom = newTempAnom2, 
                         GMTSmoother = rep(avgGMT2, numNewTemps),
                         X_c.1 = rep(combinedData$X_c.1[31+15], numNewTemps),
                         X_c.2 = rep(combinedData$X_c.2[31+15], numNewTemps),
                         X_s.1 = rep(combinedData$X_s.1[31+15], numNewTemps),
                         X_s.2 = rep(combinedData$X_s.2[31+15], numNewTemps))
  newX2 <- model.matrix(~ GMTSmoother*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                        data = newData2)
  meanTemps2 <- newX2[1,] %*% meanTemp$model$coef
  newTREFHT2 <- newTempAnom2 + meanTemps2
  newData2$TREFHT <- newTREFHT2
  quantModelPreds2 <- predict(quantModel$model, newdata = newData2)
  quantModelPreds2 <-  t(apply(quantModelPreds2, 1, sort))
  
  ## calculate empirical quantiles
  combinedData$period <- factor(ifelse(year(combinedData$date) <= 2005, "historical",
                                ifelse(year(combinedData$date) <= 2035, "f1", "f2")),
                                levels = c("historical", "f1", "f2"))
  combinedData$TempAnomGroup <- cut(combinedData$TempAnom, 
                                    breaks = quantile(combinedData$TempAnom, c(0, 1/4, 1/2, 3/4, 1)))
  
  TempAnomGroup_midpoints <- quantile(combinedData$TempAnom, c(1/8, 3/8, 5/8, 7/8))
  empirical_quantiles <- tapply(combinedData$DEWP, list(combinedData$period, combinedData$TempAnomGroup, 
                                                        month(combinedData$date)),
                                quantile, probs = quantModel$model$tau[c(3, 6, 9)]) 
  
  plot_dir <- "/research/Temp_Dewpt/Results/figures/validation/"
  f <- paste0(plot_dir, paste0("quantModel_Lon", -round(lon), 
                               "Lat", round(lat),".pdf"))
  pdf(f)
  par(mfrow = c(1,1))
  plot(DEWP ~ TREFHT, data = combinedData, cex = 0.5, subset = (year(date) %in% 1990:2005 & 
                                                                  month(date) == 7),
       xlim = c(min(TREFHT[month(date) == 7]), max(TREFHT[month(date) == 7])),
       ylim = c(min(DEWP[month(date) == 7]), max(DEWP[month(date) == 7])), col = rgb(0,0,0, alpha = 0.05),
       main = "",
       xlab = expression(paste("Temperature (", degree,"C)")), ylab = expression(paste("Dew Point (", degree,"C)")),
       axes = F)
  box()
  axis(1, lwd  = 0, lwd.ticks = 1)
  axis(2, lwd  = 0, lwd.ticks = 1)
  lines(DEWP ~ TREFHT, data = combinedData, cex = 0.5, subset = (year(date) >= 2071 & 
                                                                   year(date) <= 2080 & 
                                                                   month(date) == 7),
        type = "p", col = rgb(1,0,0, alpha = 0.05))
  abline(a = 0, b = 1, col = 'gray', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,3]), col = 'black', lty = 'dashed', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,6]), col = 'black', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,9]), col = 'black',lty = 'dashed',  lwd = 2)
  
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,3]), col = 'red', lty = 'dashed', lwd = 2)
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,6]), col = 'red', lwd = 2)
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,9]), col = 'red',lty = 'dashed',  lwd = 2)
  
  points(TempAnomGroup_midpoints + rep(meanTemps1, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[1, ,2])[1,], col = 'black', pch = 0, cex = 2)
  points(TempAnomGroup_midpoints + rep(meanTemps1, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[1, ,2])[2,], col = 'black', pch = 15, cex = 2)
  points(TempAnomGroup_midpoints + rep(meanTemps1, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[1, ,2])[3,], col = 'black', pch = 0, cex = 2)
  
  points(TempAnomGroup_midpoints + rep(meanTemps2, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[3, ,2])[1,], col = 'red', pch = 0, cex = 2)
  points(TempAnomGroup_midpoints + rep(meanTemps2, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[3, ,2])[2,], col = 'red', pch = 15, cex = 2)
  points(TempAnomGroup_midpoints + rep(meanTemps2, length(TempAnomGroup_midpoints)), 
         simplify2array(empirical_quantiles[3, ,2])[3,], col = 'red', pch = 0, cex = 2)
  
  
  legend("bottomright", col = c("black", "red", "black", "black", "black"), 
         pch = c(1,1,NA, NA,15), lty = c(NA, NA, 'solid', 'dashed',NA),
         legend = c("1990 - 2005, CESM", "2071 - 2080, CESM", "median", 
                    "5th & 95th percentiles", "empirical quantile"),
         bty = 'n')
  
  mtext(text = paste0("Lon ", round(lon, 1), 
                      ", Lat ", round(lat, 1), 
                      " (Pearson stat = ", 
                      round(quantModel$pearson_stat),")"), 
        outer = TRUE, line = -2.5, cex = 1.5)
  dev.off() 
}


plot_combinedTemp <- function(obs_data, simulTemp_fromObs, CESMmeanTemp_coef, CESM_meanCoefs_boot, 
                              delta_hat, var_delta_hat,
                              obsYearsForSimul, simulYears){
  # plot (a) projected changes in mean temperature per degree warming in GMT, (b) projected changes in temperature variability per degree warming in GMT, (c) the first year of the temperature observations and future  simulation, (d) the same but for temperature deviations, and (e) & (f) the same but for JJA annual average temperatures and deviations.
  
  # obs_data is the output of make_combinedData() for observations
  # simulTemp_fromObs is the output of simulTemp()
  # CESMmeanTemp_coef is estimate_meanTemp()$model$coef for the CESM data
  # CESM_meanCoefs_boot is estimate_meanTemp()$bootCoefs for the CESM data
  # delta_hat is estimate_logRho()$delta_hat_smoothed
  # var_delta_hat is estimate_logRho()$V_deltaSmooth
  # obsYearsForSimul are the years of obs_data used for the future simulation
  # simulYears are the years for the future simulation
  
  num_harmonics <- sum(grepl("X_c", names(obs_data)))
  nSeason <- dim(obs_data[year(obs_data$date) == first(year(obs_data$date)),])[1]
  me_ind <- which(names(CESMmeanTemp_coef) == "GMTSmoother")
  interaction_inds <- which(grepl("GMTSmoother:X_c", names(CESMmeanTemp_coef)) |
                              grepl("GMTSmoother:X_s", names(CESMmeanTemp_coef)))
  harmonics <- as.matrix(obs_data[year(obs_data$date) == first(year(obs_data$date)),
                                  grepl("X_", colnames(obs_data))])
  
  # change in mean temperature for 1 degree increase in GMT:
  Delta <- cbind(rep(1, nSeason), harmonics) %*% CESMmeanTemp_coef[c(me_ind, interaction_inds)]
  bootDelta <- cbind(rep(1, nSeason), harmonics) %*% CESM_meanCoefs_boot[c(me_ind, interaction_inds),]
  Delta_CI <- apply(bootDelta, 1, quantile, c(0.025, 0.975))
  
  # JJA annual averages
  anMeanObs <- tapply(obs_data$TREFHT[year(obs_data$date) %in% obsYearsForSimul],
                      year(obs_data$date[year(obs_data$date) %in% obsYearsForSimul]), mean)
  anMeanObsAnom <- tapply(obs_data$TempAnom[year(obs_data$date) %in% obsYearsForSimul],
                          year(obs_data$date[year(obs_data$date) %in% obsYearsForSimul]),mean)
  anMeanSimul <- tapply(simulTemp_fromObs$simul, 
                        year(obs_data$date[year(obs_data$date) %in% obsYearsForSimul]), mean)
  anMeanSimulAnom <- tapply(simulTemp_fromObs$anom, 
                            year(obs_data$date[year(obs_data$date) %in% obsYearsForSimul]), mean)
  
  
  summerDays <- seq.Date(mdy("06-01-1990"), mdy("8-31-1990"), by = 1)
  save_file <- "/research/Temp_Dewpt/Results/figures/combined_temp_plotMN.pdf"
  
  pdf(save_file, height = 10, width = 16)
  par(mfrow = c(2,3), mar = c(5,5,4,4), cex.axis = 2, cex.lab = 2)
  
  plot(summerDays, Delta, type='l',
       ylab = expression(paste("Change in local mean / ", degree, "C")), xlab = "", axes = F,
       ylim = c(min(Delta_CI), max(Delta_CI)))
  box()
  lines(summerDays, Delta_CI[1,], lty = 'dashed')
  lines(summerDays, Delta_CI[2,], lty = 'dashed')
  axis(1, lwd  = 0, lwd.ticks = 1, at = c(mdy("6-1-1990"),mdy("7-1-1990"),mdy("8-1-1990"),mdy("9-1-1990")),
       labels = c("June 1", "July 1", "August 1", "September 1"))
  axis(2, lwd  = 0, lwd.ticks = 1)
  text(par("usr")[1],par("usr")[4],"(a)",cex=2.5,adj=c(-0.15,1.15))
  
  plot(TREFHT ~ date, data = obs_data, type = 'l',
       ylim = c(min(TREFHT, simulTemp_fromObs$simul)-1.25, max(TREFHT, simulTemp_fromObs$simul)),
       subset = (year(date) == 1996), axes = F,
       ylab = expression(paste("Temperature (", degree, "C)")), xlab = "")
  lines(simulTemp_fromObs$simul[1:92] ~ obs_data$date[year(obs_data$date) == 1996], col = 'red')
  box()
  axis(1, lwd  = 0, lwd.ticks = 1,
       at = c(mdy("6-1-1996"), mdy("7-1-1996"), mdy("8-1-1996"), mdy("9-1-1996")),
       labels = c("June 1", "July 1", "August 1", "September 1"))
  axis(2, lwd  = 0, lwd.ticks = 1)
  legend("bottomleft", col = c('black', 'red'), lty = c(1, 1),
         legend = c("Observations, 1996", "Simulation, 2071"),
         cex = 1.75)
  text(par("usr")[1],par("usr")[4],"(c)",cex=2.5,adj=c(-0.15,1.15))
  
  plot(1:10, anMeanObs,
       ylim = c(min(anMeanObs, anMeanSimul)-1, max(anMeanObs, anMeanSimul)),
       type = 'l', ylab = expression(paste("Temperature (", degree, "C)")), 
       xlab = "",
       axes = FALSE)
  mtext("Year", side = 1, line = 4, cex = 1.75)
  lines(1:10, anMeanSimul, col = 'red')
  abline(h = 0, col = 'gray', lty = 'dashed')
  box()
  axis(1, seq(1, 10, by = 2), labels = seq(1996, 2005, by = 2), lwd  = 0, lwd.ticks = 1)
  axis(1, seq(1, 10, by = 2), labels = seq(2071, 2080, by = 2), col.axis = 'red', line = 1.25,
       lwd  = 0, lwd.ticks = 0)
  axis(2, lwd  = 0, lwd.ticks = 1)
  legend("bottomleft", col = c('black', 'red'), lty = c(1, 1),
         legend = c("Observations, JJA Average", "Simulation, JJA Average"),
         cex = 1.75)
  text(par("usr")[1],par("usr")[4],"(e)",cex=2.5,adj=c(-0.15,1.15))
  
  plot((1:46)/92, exp(0.5*delta_hat[2:47])-1, type='l', ylim = c(-0.01, 0.1), axes = F,
       xlab = "Period (days)", ylab = expression(paste("rel. change in square root spectral density / ", degree, "C")),
       main = "",
       xlim = c(-0.01, 0.5))
  box()
  axis(1, lwd  = 0, lwd.ticks = 1, at = c(1/92, 1/14, 1/7, 1/2), labels = c(92, 14, 7, 2))
  axis(2, lwd  = 0, lwd.ticks = 1,)
  lines((1:46)/92, exp(0.5*(delta_hat[2:47] + 2*(sqrt(var_delta_hat[2:47]))))-1, lty = 'dashed')
  lines((1:46)/92, exp(0.5*(delta_hat[2:47] - 2*(sqrt(var_delta_hat[2:47]))))-1, lty = 'dashed')
  abline(h=0, col = 'gray', lwd = 2, lty = 'dashed')
  points(0, exp(0.5*delta_hat[1])-1, col = 'black', pch  = 16)
  text(0, exp(0.5*delta_hat[1])-1, "JJA annual average", col = 'black', pos = 4, cex = 1.75)
  arrows(0, exp(0.5*(delta_hat[1] - 2*(sqrt(var_delta_hat[1]))))-1,
         0, exp(0.5*(delta_hat[1] + 2*(sqrt(var_delta_hat[1]))))-1,
         angle = 90, col = 'black', length = 0)
  text(par("usr")[1],par("usr")[4],"(b)",cex=2.5,adj=c(-0.15,1.15))
  
  plot(TempAnom ~ date, data = obs_data, type = 'l',
       ylim = c(min(TempAnom, simulTemp_fromObs$anom), max(TempAnom, simulTemp_fromObs$anom)),
       subset = (year(date) == 1996), axes = F,
       ylab = expression(paste("Temperature Deviation (", degree, "C)")), xlab = "")
  lines(simulTemp_fromObs$anom[1:92] ~ obs_data$date[year(GSOD$date) == 1996], col = 'red')
  box()
  axis(1, lwd  = 0, lwd.ticks = 1,
       at = c(mdy("6-1-1996"), mdy("7-1-1996"), mdy("8-1-1996"), mdy("9-1-1996")),
       labels = c("June 1", "July 1", "August 1", "September 1"))
  axis(2, lwd  = 0, lwd.ticks = 1)
  abline(h = 0, col = 'gray', lty = 'dashed')
  legend("bottomleft", col = c('black', 'red'), lty = c(1, 1),
         legend = c("Observations, 1996", "Simulation, 2071"), cex = 1.75)
  text(par("usr")[1],par("usr")[4],"(d)",cex=2.5,adj=c(-0.15,1.15))
  
  plot(1:10, anMeanObsAnom,
       ylim = c(min(anMeanObsAnom, anMeanSimulAnom), max(anMeanObsAnom, anMeanSimulAnom)),
       type = 'l', ylab = expression(paste("Temperature Deviation (", degree, "C)")), 
       xlab = "",
       axes = FALSE)
  mtext("Year", side = 1, line = 4, cex = 1.75)
  lines(1:10, anMeanSimulAnom, col = 'red')
  abline(h = 0, col = 'gray', lty = 'dashed')
  box()
  axis(1, seq(1, 10, by = 2), labels = seq(1996, 2005, by = 2), lwd  = 0, lwd.ticks = 1)
  axis(1, seq(1, 10, by = 2), labels = seq(2071, 2080, by = 2), col.axis = 'red', line = 1.25,
       lwd  = 0, lwd.ticks = 0)
  axis(2, lwd  = 0, lwd.ticks = 1)
  legend("bottomleft", col = c('black', 'red'), lty = c(1, 1),
         legend = c("Observations, JJA Average", "Simulation, JJA Average"),
         cex = 1.75)
  text(par("usr")[1],par("usr")[4],"(f)",cex=2.5,adj=c(-0.15,1.15))
  dev.off()
}

plot_DewpVsTemp <- function(obs_data, CESM_data, quantModel_OBS, quantModel_CESM, 
                            meanTemp_obs, meanTemp_CESM, simulTemp_fromObs, 
                            simulDewp_fromObs, obsYearsForSimul, simulYears){
  # plot dew point vs. temperature and dew point vs. temperature deviation for CESM and for observations and observation-based simulation, along with quantile curves. Plots are for July data in 1996-2005 and 2071-2080, and quantile curves are for mid-July when GMT is equal to the average value in the time period.
  
  # obs_data is the output of make_combinedData() for observations
  # CESM_data is the output of make_combinedData() for CESM-LE
  # quantModel_OBS is the output of estimate_quantModel(obs_data)
  # quantModel_CESM is the output of estimate_quantModel(CESM_data)
  # meanTemp_obs is the output of estimate_meanTemp(obs_data)
  # meanTemp_CESM is the output of estimate_meanTemp(CESM_data)
  # simulDewp_fromObs is the output of simulTemp(...)
  # simulTemp_fromObs is the output of simulDEWP(...)
  # obsYearsForSimul are the years of obs_data used for the future simulation
  # simulYears are the years for the future simulation
  
  
  GMT_coefs <- coef(quantModel_CESM$model)[rownames(coef(quantModel_CESM$model)) == "GMTSmoother", ]
  
  numNewTemps <- 100
  newTempAnom1 <- seq(quantile(CESM_data$TempAnom, 0.01), quantile(CESM_data$TempAnom, 0.99), 
                      len = numNewTemps)
  avgGMT1 <- mean(CESM_data$GMTSmoother[year(CESM_data$date) %in%  obsYearsForSimul])
  newData1 <- data.frame(TempAnom = newTempAnom1, 
                         GMTSmoother = rep(avgGMT1, numNewTemps),
                         X_c.1 = rep(CESM_data$X_c.1[31+15], numNewTemps),
                         X_c.2 = rep(CESM_data$X_c.2[31+15], numNewTemps),
                         X_s.1 = rep(CESM_data$X_s.1[31+15], numNewTemps),
                         X_s.2 = rep(CESM_data$X_s.2[31+15], numNewTemps))
  newX1 <- model.matrix(~ GMTSmoother*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                        data = newData1)
  newTREFHT1 <- newTempAnom1 + newX1 %*% meanTemp_CESM$model$coef
  newData1$TREFHT <- newTREFHT1
  quantModelPreds1 <- predict(quantModel_CESM$model, newdata = newData1)
  quantModelPreds1 <-  t(apply(quantModelPreds1, 1, sort))
  
  newTempAnom2 <- newTempAnom1
  avgGMT2 <- mean(CESM_data$GMTSmoother[year(CESM_data$date) %in%  simulYears])
  newData2 <- data.frame(TempAnom = newTempAnom2, 
                         GMTSmoother = rep(avgGMT2, numNewTemps),
                         X_c.1 = rep(CESM_data$X_c.1[31+15], numNewTemps),
                         X_c.2 = rep(CESM_data$X_c.2[31+15], numNewTemps),
                         X_s.1 = rep(CESM_data$X_s.1[31+15], numNewTemps),
                         X_s.2 = rep(CESM_data$X_s.2[31+15], numNewTemps))
  newX2 <- model.matrix(~ I(GMTSmoother)*(X_c.1 + X_c.2 + X_s.1 + X_s.2), 
                        data = newData2)
  
  newTREFHT2 <- newTempAnom2 + newX2 %*% meanTemp_CESM$model$coef
  newData2$TREFHT <- newTREFHT2
  quantModelPreds2 <- predict(quantModel_CESM$model, newdata = newData2)
  quantModelPreds2 <-  t(apply(quantModelPreds2, 1, sort))
  
  
  newData1_obs <- data.frame(TempAnom = newTempAnom1, 
                             GMTSmoother = rep(avgGMT1, numNewTemps),
                             X_c.1 = rep(CESM_data$X_c.1[31+15], numNewTemps),
                             X_c.2 = rep(CESM_data$X_c.2[31+15], numNewTemps),
                             X_s.1 = rep(CESM_data$X_s.1[31+15], numNewTemps),
                             X_s.2 = rep(CESM_data$X_s.2[31+15], numNewTemps))
  newTREFHT1_obs <- newTempAnom1 + newX1 %*% meanTemp_obs$model$coef
  newData1_obs$TREFHT <- newTREFHT1_obs
  quantModelPreds1_obs <- predict(quantModel_OBS$model, newdata = newData1_obs) + 
    outer(newData1_obs$GMTSmoother, GMT_coefs)
  quantModelPreds1_obs <-  t(apply(quantModelPreds1_obs, 1, sort))
  
  newData2_obs <- data.frame(TempAnom = newTempAnom2, 
                             GMTSmoother = rep(avgGMT2, numNewTemps),
                             X_c.1 = rep(CESM_data$X_c.1[31+15], numNewTemps),
                             X_c.2 = rep(CESM_data$X_c.2[31+15], numNewTemps),
                             X_s.1 = rep(CESM_data$X_s.1[31+15], numNewTemps),
                             X_s.2 = rep(CESM_data$X_s.2[31+15], numNewTemps))
  
  Delta_forIllustObs <- mean(meanTemp_CESM$model$fitted[year(CESM_data$date) %in% simulYears & 
                                                          month(CESM_data$date) == 7 & 
                                                          day(CESM_data$date) == 1 &
                                                          CESM_data$run == 1]) - 
    mean(meanTemp_CESM$model$fitted[year(CESM_data$date) %in% obsYearsForSimul & 
                                      month(CESM_data$date) == 7 & 
                                      day(CESM_data$date) == 1 &
                                      CESM_data$run == 1])
  
  newTREFHT2_obs <- newTempAnom2 +  newX1 %*% meanTemp_obs$model$coef + Delta_forIllustObs
  newData2_obs$TREFHT <- newTREFHT2_obs
  quantModelPreds2_obs <- predict(quantModel_OBS$model, newdata = newData2_obs) + 
    outer(newData2_obs$GMTSmoother, GMT_coefs)
  quantModelPreds2_obs <-  t(apply(quantModelPreds2_obs, 1, sort))
  
  pdf('/research/Temp_Dewpt/Results/figures/DEWPvsTemp_MN.pdf', width = 10, height = 10)
  par(mfrow = c(2,2), cex.axis = 1.75, cex.lab = 1.75, mar = c(5,5,4,4))
  plot(DEWP ~ TREFHT, data = CESM_data, cex = 0.5, subset = (year(date) %in% obsYearsForSimul &
                                                               month(date) == 7),
       xlim = c(min(TREFHT), max(TREFHT)),
       ylim = c(min(DEWP), max(DEWP)), col = rgb(0,0,0, alpha = 0.05),
       main = "",
       xlab = expression(paste("Temperature (", degree,"C)")), ylab = expression(paste("Dew Point (", degree,"C)")),
       axes = F)
  box()
  axis(1, lwd  = 0, lwd.ticks = 1)
  axis(2, lwd  = 0, lwd.ticks = 1)
  lines(DEWP ~ TREFHT, data = CESM_data, cex = 0.5, subset = (year(date) %in% simulYears &
                                                                month(date) == 7),
        type = "p", col = rgb(1,0,0, alpha = 0.05))
  abline(a = 0, b = 1, col = 'gray', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,3]), col = 'black', lty = 'dashed', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,6]), col = 'black', lwd = 2)
  lines(newTREFHT1, newTREFHT1 - exp(quantModelPreds1[,9]), col = 'black',lty = 'dashed',  lwd = 2)
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,3]), col = 'red', lty = 'dashed', lwd = 2)
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,6]), col = 'red', lwd = 2)
  lines(newTREFHT2, newTREFHT2 - exp(quantModelPreds2[,9]), col = 'red',lty = 'dashed',  lwd = 2)
  legend("bottomright", col = c("black", "red", "black", "black"), 
         pch = c(1,1,NA, NA), lty = c(NA, NA, 'solid', 'dashed'),
         legend = c("1996 - 2005, CESM", "2071 - 2080, CESM", "median", "5th & 95th percentiles"), cex = 1.25)
  
  plot(DPS ~ TempAnom, data = CESM_data, cex = 0.5, subset = (year(date) %in% obsYearsForSimul &
                                                                month(date) == 7),
       xlim = c(min(TempAnom), max(TempAnom)),
       ylim = (c(min(DPS[DPS != 0]), max(DPS))), col = rgb(0,0,0, alpha = 0.05),
       main = "",
       xlab = expression(paste("Temperature - mean trend (", degree, "C)")), 
       ylab = expression(paste("Temperature - Dew Point (", degree,"C)")), log = "y", axes = FALSE)
  box()
  axis(1, lwd  = 0, lwd.ticks = 1)
  axis(2, lwd  = 0, lwd.ticks = 1)
  lines(DPS ~ TempAnom, data = CESM_data, cex = 0.5, subset = (year(date) %in% simulYears &
                                                                 month(date) == 7),
        type = "p", col = rgb(1,0,0, alpha = 0.05))
  lines(newTempAnom1, exp(quantModelPreds1[,3]), col = 'black', lty = 'dashed', lwd = 2)
  lines(newTempAnom1, exp(quantModelPreds1[,6]), col = 'black', lwd = 2)
  lines(newTempAnom1, exp(quantModelPreds1[,9]), col = 'black',lty = 'dashed',  lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2[,3]), col = 'red', lty = 'dashed', lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2[,6]), col = 'red', lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2[,9]), col = 'red',lty = 'dashed',  lwd = 2)
  legend("bottomright", col = c("black", "red", "black", "black"), 
         pch = c(1,1,NA, NA), lty = c(NA, NA, 'solid', 'dashed'),
         legend = c("1996 - 2005, CESM", "2071 - 2080, CESM", "median", "5th & 95th percentiles"), cex = 1.25)
  
  plot(DEWP ~ TREFHT, data = obs_data, cex = 0.5, subset = (year(date) %in% obsYearsForSimul &
                                                              month(date) == 7),
       xlim = c(min(TREFHT, simulTemp_fromObs$simul), max(TREFHT, simulTemp_fromObs$simul)),
       ylim = c(min(DEWP, simulDewp_fromObs), max(DEWP, simulDewp_fromObs)), 
       col = rgb(0,0,0, alpha = 0.5),
       main = "",
       xlab = expression(paste("Temperature (", degree,"C)")), 
       ylab = expression(paste("Dew Point (", degree,"C)")),
       axes = F)
  box()
  axis(1, lwd  = 0, lwd.ticks = 1)
  axis(2, lwd  = 0, lwd.ticks = 1)
  simulDates <- CESM_data$date[year(CESM_data$date) %in% simulYears & CESM_data$run == 1]
  lines(simulDewp_fromObs ~ simulTemp_fromObs$simul, cex = 0.5, subset = (month(simulDates) == 7),
        type = "p", col = rgb(1,0,0, alpha = 0.5))
  abline(a = 0, b = 1, col = 'gray', lwd = 2)
  lines(newTREFHT1_obs, newTREFHT1_obs - exp(quantModelPreds1_obs[,3]), col = 'black', lty = 'dashed', lwd = 2)
  lines(newTREFHT1_obs, newTREFHT1_obs - exp(quantModelPreds1_obs[,6]), col = 'black', lwd = 2)
  lines(newTREFHT1_obs, newTREFHT1_obs - exp(quantModelPreds1_obs[,9]), col = 'black',lty = 'dashed',  lwd = 2)
  lines(newTREFHT2_obs, newTREFHT2_obs - exp(quantModelPreds2_obs[,3]), col = 'red', lty = 'dashed', lwd = 2)
  lines(newTREFHT2_obs, newTREFHT2_obs - exp(quantModelPreds2_obs[,6]), col = 'red', lwd = 2)
  lines(newTREFHT2_obs, newTREFHT2_obs - exp(quantModelPreds2_obs[,9]), col = 'red',lty = 'dashed',  lwd = 2)
  legend("bottomright", col = c("black", "red", "black", "black"), 
         pch = c(1,1,NA, NA), lty = c(NA, NA, 'solid', 'dashed'),
         legend = c("1996 - 2005, GSOD", "2071 - 2080, simulation", "median", "5th & 95th percentiles"), cex = 1.25)
  
  plot(DPS ~ TempAnom, data = obs_data, cex = 0.5,  subset = (year(date) %in% obsYearsForSimul &
                                                                month(date) == 7),
       xlim = c(min(TempAnom, simulTemp_fromObs$anom), max(TempAnom, simulTemp_fromObs$anom)),
       ylim = (c(min(DPS[DPS != 0], simulDewp_fromObs[simulDewp_fromObs != 0]),
                 max(DPS[DPS != 0], simulDewp_fromObs[simulDewp_fromObs != 0]))),
       col = rgb(0,0,0, alpha = 0.5),
       main = "",
       xlab = expression(paste("Temperature - mean trend (", degree, "C)")), 
       ylab = expression(paste("Temperature - Dew Point (", degree,"C)")), log = "y", 
       axes = FALSE)
  box()
  axis(1, lwd  = 0, lwd.ticks = 1)
  axis(2, lwd  = 0, lwd.ticks = 1)
  lines(I(simulTemp_fromObs$simul - simulDewp_fromObs) ~ simulTemp_fromObs$anom, 
        subset = (month(simulDates) == 7),
        cex = 0.5,
        type = "p", col = rgb(1,0,0, alpha = 0.5))
  lines(newTempAnom1, exp(quantModelPreds1_obs[,3]), col = 'black', lty = 'dashed', lwd = 2)
  lines(newTempAnom1, exp(quantModelPreds1_obs[,6]), col = 'black', lwd = 2)
  lines(newTempAnom1, exp(quantModelPreds1_obs[,9]), col = 'black',lty = 'dashed',  lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2_obs[,3]), col = 'red', lty = 'dashed', lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2_obs[,6]), col = 'red', lwd = 2)
  lines(newTempAnom2, exp(quantModelPreds2_obs[,9]), col = 'red',lty = 'dashed',  lwd = 2)
  legend("bottomright", col = c("black", "red", "black", "black"), 
         pch = c(1,1,NA, NA), lty = c(NA, NA, 'solid', 'dashed'),
         legend = c("1996 - 2005, GSOD", "2071 - 2080, simulation", "median", "5th & 95th percentiles"), cex = 1.25)
  dev.off()
}

plotGMTSmoother <- function(GMT){
  # plot ensemble average annual average GMT along with smoothed lowess estimate of forced trend.
  
  GMT$GMT <- GMT$GMT - mean(GMT$GMT[GMT$GMT_years %in% 1990:2005])
  pdf("/research/Temp_Dewpt/Results/figures/GMT_smoother_CESM.pdf")
  plot(GMT ~ GMT_years, data = GMT, pch = 16, 
       col = ifelse(GMT_years %in% c(1990:2005, 2026:2035, 2071:2080), 'pink', 'gray'),
       xlab = "Year", ylab = expression(paste("GMT Anomaly (", degree, "C)")), 
       main = "CESM-LE", axes = F, xlim = c(1920, 2100))
  box()
  axis(1, lwd  = 0, lwd.ticks = 1, 
       at = seq(from = 1950, to = 2100, length.out = 4))
  axis(2, lwd  = 0, lwd.ticks = 1)
  
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 1990:2005, col = 'red')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 2026:2035, col = 'red')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 2071:2080, col = 'red')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 1920:1989, col = 'black')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 2006:2025, col = 'black')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 2036:2070, col = 'black')
  lines(GMTSmoother_year ~ GMT_years, data = GMT, lwd = 2, 
        subset = GMT_years %in% 2081:2100, col = 'black')
  legend("topleft", pch = c(16, NA), lty = c(NA, 1), legend = c("Ensemble mean", "LOWESS smoother"))
  dev.off()
}