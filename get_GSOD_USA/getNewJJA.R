#nohup R CMD BATCH getNewJJA.R getNewJJA.Rout &

library(GSODR)
library(dplyr)
library(lubridate)
USA <- get_GSOD(2018, country = "United States")
station_list <- unique(USA$STNID[as.numeric(substr(USA$BEGIN, 1, 4)) <= 1973 & 
                                   as.numeric(substr(USA$END, 1, 4)) >= 2018])

date_df = data.frame(YEARMODA = seq.Date(from = ymd("1973-01-01"), to = ymd("2018-12-31"), by = 1))
date_df_JJA <- date_df %>% filter(month(YEARMODA) >= 6 & month(YEARMODA) <= 8)


keep <- logical(length(station_list))
tic <- Sys.time()
for(s in 1:length(station_list)){
  
  station_data <- get_GSOD(1973:2018, station = station_list[s])
  station_JJA <- station_data %>% filter(MONTH >= 6 & MONTH <= 8 & 
                                           as.numeric(TEMP_ATTRIBUTES) >= 4 & as.numeric(DEWP_ATTRIBUTES) >= 4)
  station_JJA <- full_join(station_JJA, date_df_JJA, by = "YEARMODA")
  station_JJA <- station_JJA %>% select(-c("YEAR", "MONTH", "DAY", "YDAY"))
  pctMissing <- station_JJA %>% group_by(year(YEARMODA)) %>% summarise(missingTemp = mean(is.na(TEMP)), 
                                                             missingDewp = mean(is.na(DEWP)))
  keep[s] <- (max(pctMissing$missingTemp[1:3]) < 1) & (max(pctMissing$missingTemp[44:46]) < 1) &
          (mean(pctMissing$missingTemp >= 0.2) <= 0.2) &
          (max(pctMissing$missingDewp[1:3]) < 1) & (max(pctMissing$missingDewp[44:46]) < 1) &
          (mean(pctMissing$missingDewp >= 0.2) <= 0.2)
 
  if(keep[s]){
    f <- paste0("/research/Temp_Dewpt/GSOD_USA/newJJA/", station_list[s], ".csv")
    write.csv(station_JJA, file = f)
  }
  print(s)
  flush.console()
}

Sys.time() - tic


f <- list.files("/research/Temp_Dewpt/GSOD_USA/newJJA")
f <- f[f != "metadata_new.csv"]
lats <- numeric(length(f))
lons <- numeric(length(f))
stnids <-  character(length(f))
names <-  character(length(f))
states <-  character(length(f))
n_missing <- numeric(length(f))
for(s in 1:length(f)){
  GSOD_loc <- read.csv(paste0("/research/Temp_Dewpt/GSOD_USA/newJJA/",f[s]), header=TRUE)
  lats[s] <- GSOD_loc$LATITUDE[1]
  lons[s] <- GSOD_loc$LONGITUDE[1]
  stnids[s] <- as.character(GSOD_loc$STNID[1])
  states[s] <- as.character(GSOD_loc$STATE[1])
  names[s] <- as.character(GSOD_loc$NAME[1])
  n_missing[s] <- sum(is.na(GSOD_loc$TEMP) | is.na(GSOD_loc$DEWP))
}

metadata <- data.frame(LAT = lats, LON = lons, STNID = stnids, STATE = states, NAME = names, n_missing = n_missing)
write.csv(metadata, file = "/research/Temp_Dewpt/GSOD_USA/newJJA/metadata_new.csv", row.names = FALSE)