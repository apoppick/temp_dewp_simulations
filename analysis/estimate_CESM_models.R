#nohup R CMD BATCH estimate_CESM_models.R estimate_CESM_models.Rout &

# This code estimates changes in mean temperature, temperature variability, and dew point quantiles needed to produce the observation-based simulation. It also outputs several figures illustrating the behavior of these models at an example location near Minneapolis and at several locations of varying qualities of model fit.


rm(list=ls())

library(MASS)
library(quantreg)
library(splines)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(MBC)

source("/research/Temp_Dewpt/Code/analysis/helper_functions.R")
source("/research/Temp_Dewpt/Code/analysis/plotting_functions.R")

load('/research/comps1617/processed_data/RCP85RunsSmall.RData') #Global Mean Temperatures
TGlobal <- TGlobal - 273.15
TGlobal <- TGlobal[,1:35] # remove UToronto runs
rm(list = c("fname", "sname"))

load('/research/Temp_Dewpt/LENS/DEWP_TREFHT_calculated6Hour_LENS_CONUS.RData') #CESM-LE data
TREFHT_6hourly <- TREFHT - 273.15
DEWPREFHT <- DEWPREFHT - 273.15
rm("TREFHT")

dates_6hourlyData <- c(seq.Date(mdy("01-01-1990"), mdy("12-31-2005"), by = 1),
                       seq.Date(mdy("01-01-2026"), mdy("12-31-2035"), by = 1),
                       seq.Date(mdy("01-01-2071"), mdy("12-31-2080"), by = 1))
dates_6hourlyData <- dates_6hourlyData[!(month(dates_6hourlyData) == 2 & day(dates_6hourlyData) == 29)]
JJADates_6hourly <- dates_6hourlyData[month(dates_6hourlyData) >= 6 & month(dates_6hourlyData) <= 8]
season_ind <-(month(dates_6hourlyData[1:365]) %in% 6:8)

nLon <- length(lon_convert_subs)
nLat <- length(lat_subs)
nRuns <- dim(TREFHT_6hourly)[4]
nSeason <- sum(season_ind)
nYears_6hourly <- length(dates_6hourlyData)/365

## Get GMT smoother
GMT <- get_GMTSmooth(TGlobal)
plotGMTSmoother(GMT)

MSP_loc <- data.frame(L = 22, l = 26) #for Section 3 and 4 plots
meanDiagnostics_ind <- c(634, 159, 996) #worst I1, best I1, median I1
varDiagnostics_ind <- c(643, 1053, 964) # worst deviance, best deviance, median deviance
quantModelDiagnostics_ind <- c(301, 1283, 896) #worst pearson, best pearson, median pearson

numCores <- 4 #detectCores()
cl <- registerDoParallel(numCores) 
tic <- Sys.time()
CESM_models <- foreach(loc = 1:(nLon * nLat)) %dopar% {
  lonlat <- arrayInd(loc, c(nLon, nLat))
  l <- lonlat[1]
  L <- lonlat[2]
  TREFHT_JJA_loc <- TREFHT_6hourly[l,L,month(dates_6hourlyData) >= 6 & 
                                     month(dates_6hourlyData) <= 8,]
  DEWPREFHT_JJA_loc <- DEWPREFHT[l,L,month(dates_6hourlyData) >= 6 & 
                                   month(dates_6hourlyData) <= 8,]
  
  combined_data <- make_combinedData(TREFHT_JJA_loc, DEWPREFHT_JJA_loc, GMT,
                                     JJADates_6hourly, 
                                     season_ind, num_harmonics = 2)
  
  ## Estimate mean temperature model and add anomaly to dataset
  meanModel_boot <- ifelse(l == MSP_loc$l & L == MSP_loc$L, TRUE, FALSE)
  meanTemp <- estimate_meanTemp(combined_data, boot = meanModel_boot)
  combined_data$TempAnom <- resid(meanTemp$model)
  if(loc %in% meanDiagnostics_ind){
    meanDiagnostics(combined_data, meanTemp, lon = lon_convert_subs[l], lat = lat_subs[L])
  }
  
  ## Estimate temperature variability change model
  logRhoEst <- estimate_logRho(combined_data)
  if(loc %in% varDiagnostics_ind){
    varDiagnostics(combined_data, logRhoEst, lon = lon_convert_subs[l], lat = lat_subs[L])
  }
  
  ## Estimate quantile regression model
  quantModel <- estimate_quantModel(combined_data)
  if(loc %in% quantModelDiagnostics_ind){
    quantModelDiagnostics(combined_data, quantModel, meanTemp, lon = lon_convert_subs[l], lat = lat_subs[L])
  }
  
  ## Calculate changes in risk probabilities
  riskProbs <- calculate_riskProbs(combined_data, quantModel, meanTemp)
  
  list(meanTemp_coef = meanTemp$model$coef, 
       meanTemp_bootCoefs = meanTemp$bootCoefs,
       R2 = summary(meanTemp$model)$r.squared, 
       I0 = meanTemp$I0, 
       I1 = meanTemp$I1,
       delta_hat = logRhoEst$delta_hat, 
       V_deltaSmooth = logRhoEst$V_deltaSmooth, 
       deviance = logRhoEst$deviance,
       quantModelCESMCoef = coef(quantModel$model),
       err = quantModel$err,
       riskProbs = riskProbs,
       quantModel_pearson_stat = quantModel$pearson_stat)
}
stopImplicitCluster()
Sys.time() - tic

save(file = "/research/Temp_Dewpt/Results/data_output/estimate_CESM_models.RData", 
     list = c("CESM_models", "lon_convert_subs", "lat_subs"))
q(save = "no")