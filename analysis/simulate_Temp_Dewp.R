# this code takes the estimated changes in temperature and dew point distributions from estimate_CESM_models and produces observation-based future simulations of temperature and dew point for years 2071-2080 of RCP8.5, using GSOD observations from 1996-2005 for the simulation.

#nohup R CMD BATCH simulate_Temp_Dewp.R simulate_Temp_Dewp.Rout &
rm(list=ls())

set.seed(309456354)

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
library(fields)

load('/research/Temp_Dewpt/Results/data_output/estimate_CESM_models.RData')
source("/research/Temp_Dewpt/Code/analysis/helper_functions.R")
source("/research/Temp_Dewpt/Code/analysis/plotting_functions.R")

load('/research/Temp_Dewpt/LENS/landfrac_CONUS.RData')

load('/research/comps1617/processed_data/RCP85RunsSmall.RData') #Global Mean Temperatures
TGlobal <- TGlobal - 273.15
TGlobal <- TGlobal[,1:35] # remove UToronto runs
rm(list = c("fname", "sname"))

lonlat <- expand.grid(lon_convert_subs, lat_subs)
landfrac_CONUS_vec <- c(landfrac_CONUS)

GMT <- get_GMTSmooth(TGlobal)

simulYears <- 2071:2080
obsYearsForSimul <- 1996:2005

GMT_year_change_forSimulObs <- GMT$GMTSmoother_year[GMT$GMT_years %in% simulYears] - 
                                GMT$GMTSmoother_year[GMT$GMT_years %in% obsYearsForSimul]

GSOD_metadata <- read.csv("/research/Temp_Dewpt/GSOD_USA/newJJA/metadata_new.csv", header = TRUE)
GSOD_metadata <- subset(GSOD_metadata, n_missing == 0 & LON >= - 125 & LON <= -66 & LAT >= 25 & LAT <= 50)

nStations <- dim(GSOD_metadata)[1]
nYearsSimul <- length(simulYears)
nLat <- length(lat_subs)
nLon <- length(lon_convert_subs)

season_ind <-(month(seq.Date(ymd("1990-1-1"), ymd("1990-12-31"), by = 1)) %in% 6:8)


## Load MN CESM data (needed to make plots for Sections 3 and 4)
load('/research/Temp_Dewpt/LENS/DEWP_TREFHT_calculated6Hour_LENS_CONUS.RData')
TREFHT_6hourly <- TREFHT - 273.15
DEWPREFHT <- DEWPREFHT - 273.15
dates_6hourlyData <- c(seq.Date(mdy("01-01-1990"), mdy("12-31-2005"), by = 1),
                       seq.Date(mdy("01-01-2026"), mdy("12-31-2035"), by = 1),
                       seq.Date(mdy("01-01-2071"), mdy("12-31-2080"), by = 1))
dates_6hourlyData <- dates_6hourlyData[!(month(dates_6hourlyData) == 2 & day(dates_6hourlyData) == 29)]
JJADates_6hourly <- dates_6hourlyData[month(dates_6hourlyData) >= 6 & month(dates_6hourlyData) <= 8]


MSP_loc <- data.frame(L = 22, l = 26) #ind = 117
TREFHT_JJA_MSP <- TREFHT_6hourly[MSP_loc$l,MSP_loc$L,month(dates_6hourlyData) >= 6 & 
                                   month(dates_6hourlyData) <= 8,]
DEWPREFHT_JJA_MSP <- DEWPREFHT[MSP_loc$l,MSP_loc$L,month(dates_6hourlyData) >= 6 & 
                                 month(dates_6hourlyData) <= 8,]

MSP_data <- make_combinedData(TREFHT_JJA_MSP, DEWPREFHT_JJA_MSP, GMT,
                                   JJADates_6hourly, 
                                   season_ind, num_harmonics = 2)
rm("TREFHT", "TREFHT_6hourly", "DEWPREFHT", "dates_6hourlyData", "JJADates_6hourly")


# numCores <- 10 #detectCores()
# cl <- registerDoParallel(numCores)
tic <- Sys.time()
GSOD_simul <- foreach(s = 1:nStations) %do% {
  dists <- rdist.earth(x1 = cbind(GSOD_metadata$LON[s], GSOD_metadata$LAT[s]), 
                       x2 = lonlat)
  dists_land <- dists * ifelse(landfrac_CONUS_vec >= 0.5, 1, NA)
  closest_ind <- which.min(dists_land)
  closest_ind_arrayInd <- arrayInd(closest_ind, c(nLon, nLat))
  l <- closest_ind[1]
  L <- closest_ind[2]
  
  ## Load GSOD data
  GSOD <- read.csv(paste0("/research/Temp_Dewpt/GSOD_USA/newJJA/", as.character(GSOD_metadata$STNID[s]), '.csv'),
                   header = TRUE)
  
  GSOD <- make_combinedData(GSOD$TEMP, GSOD$DEWP, GMT, ymd(GSOD$YEARMODA), 
                                season_ind, num_harmonics = 2)
  
  ## Estimate mean temperature model and add anomaly to dataset
  meanTemp_GSOD <- estimate_meanTemp(GSOD, I.scores = FALSE)
  GSOD$TempAnom <- resid(meanTemp_GSOD$model)
  
  ## Estimate quantile regression model
  quantModel_OBS <- estimate_quantModel(GSOD, obs = TRUE, pearson_stat = FALSE)
  
  ## simulate temperatures
  simulTemp_fromObs <- simulTemp(GSOD, meanTemp_GSOD, GMT, CESM_models[[closest_ind]]$meanTemp_coef, 
                                 CESM_models[[closest_ind]]$delta_hat, obsYearsForSimul, simulYears)
  
  ## simulate dew points
  simulDewp_fromObs <- simulDEWP(GSOD, quantModel_OBS, GMT, CESM_models[[closest_ind]]$quantModelCESMCoef,
                                 simulTemp_fromObs, obsYearsForSimul, simulYears)
  simulDPS <- simulTemp_fromObs$simul - simulDewp_fromObs
    
  ## calculate change in risk probabilities
  riskProbs_forSimul <- calculate_riskProbs_forSimul(GSOD, meanTemp_GSOD, quantModel_OBS, GMT,
                                                     CESM_models[[closest_ind]]$quantModelCESMCoef, 
                                                     CESM_models[[closest_ind]]$meanTemp_coef,
                                                     simulTemp_fromObs)
    
    # plotting for Minneapolis example
  if(grepl("MINNEAPOLIS", GSOD_metadata$NAME[s])){
    plot_combinedTemp(GSOD, simulTemp_fromObs,  CESM_models[[closest_ind]]$meanTemp_coef, 
                      CESM_models[[closest_ind]]$meanTemp_bootCoef,
                      CESM_models[[closest_ind]]$delta_hat, CESM_models[[closest_ind]]$V_deltaSmooth, 
                      obsYearsForSimul, simulYears)
    
    #re-estimate CESM quantile model to have everything needed
    meanTemp_CESM <- estimate_meanTemp(MSP_data, boot = FALSE)
    MSP_data$TempAnom <- resid(meanTemp_CESM$model)
    quantModel_CESM <- estimate_quantModel(MSP_data, obs = FALSE, pearson_stat = FALSE)
    plot_DewpVsTemp(GSOD, MSP_data, quantModel_OBS, quantModel_CESM, 
                    meanTemp_GSOD, meanTemp_CESM,
                    simulTemp_fromObs, simulDewp_fromObs, obsYearsForSimul, simulYears)
  }
  

  list(STNID = GSOD_metadata$STNID[s],
       closest_lat = lat_subs[L], 
       closest_lon = lon_convert_subs[l],
       closest_ind = closest_ind,
       TempSimul = simulTemp_fromObs,
       DewpSimul = simulDewp_fromObs,
       riskProbs_forSimul = riskProbs_forSimul)
}
Sys.time() - tic

save(file = "/research/Temp_Dewpt/Results/data_output/simulate_fromGSOD.RData", 
     list = c("GSOD_simul", "GSOD_metadata"))
q(save = "no")
