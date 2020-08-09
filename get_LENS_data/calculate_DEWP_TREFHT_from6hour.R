## This file calculates daily average TREFHT and DEWPREFHT using 6-hourly data from CESM-LE. DEWPREFHT is calculated from TREFHT, PS, and QREFHT as described in the paper.

rm(list=ls())
library(ncdf4)
library(lubridate)
library(stringr)
source('/glade/u/home/apoppick/TempDewp_code/conversion_formulae.R')

runs <- c(1:35, 101:105)
dates <- c(seq.Date(mdy("01-01-1990"), mdy("12-31-2005"), by = 1),
           seq.Date(mdy("01-01-2026"), mdy("12-31-2035"), by = 1),
           seq.Date(mdy("01-01-2071"), mdy("12-31-2080"), by = 1))
dates <- dates[!(month(dates) == 2 & day(dates) == 29)]
dates_6hour <- rep(dates, rep(4, length(dates)))

dir <- "/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/hourly6/"

DEWPREFHT <- array(dim = c(48, 28, length(dates), length(runs)))
TREFHT <- array(dim = c(48, 28, length(dates), length(runs)))
for(r in 1:length(runs)){
  if(runs[r] %in% c(1:35, 104)){
    PS_6hour_f1 <- paste0(dir, "/PS/b.e11.B20TRC5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                          ".cam.h2.PS.1990010100Z-2005123118Z.nc")
  } else {
    PS_6hour_f1 <- paste0(dir, "/PS/b.e11.B20TRC5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                          ".cam.h2.PS.1920010100Z-2005123118Z.nc")
  }
  
  PS_6hour_f2 <- paste0(dir, "/PS/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                        ".cam.h2.PS.2026010100Z-2035123118Z.nc")
  PS_6hour_f3 <- paste0(dir, "/PS/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                        ".cam.h2.PS.2071010100Z-2080123118Z.nc")
  
  QREFHT_6hour_f1 <- paste0(dir, "/QREFHT/b.e11.B20TRC5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.QREFHT.1990010100Z-2005123118Z.nc")
  QREFHT_6hour_f2 <- paste0(dir, "/QREFHT/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.QREFHT.2026010100Z-2035123118Z.nc")
  QREFHT_6hour_f3 <- paste0(dir, "/QREFHT/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.QREFHT.2071010100Z-2080123118Z.nc")
  
  TREFHT_6hour_f1 <- paste0(dir, "/TREFHT/b.e11.B20TRC5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.TREFHT.1990010100Z-2005123118Z.nc")
  TREFHT_6hour_f2 <- paste0(dir, "/TREFHT/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.TREFHT.2026010100Z-2035123118Z.nc")
  TREFHT_6hour_f3 <- paste0(dir, "/TREFHT/b.e11.BRCP85C5CNBDRD.f09_g16.",str_pad(runs[r], 3, side = "left", pad = "0"),
                            ".cam.h2.TREFHT.2071010100Z-2080123118Z.nc")
  
  nc_PS_6hour_f1 <- nc_open(PS_6hour_f1)
  nc_PS_6hour_f2 <- nc_open(PS_6hour_f2)
  nc_PS_6hour_f3 <- nc_open(PS_6hour_f3)
  nc_QREFHT_6hour_f1 <- nc_open(QREFHT_6hour_f1)
  nc_QREFHT_6hour_f2 <- nc_open(QREFHT_6hour_f2)
  nc_QREFHT_6hour_f3 <- nc_open(QREFHT_6hour_f3)
  nc_TREFHT_6hour_f1 <- nc_open(TREFHT_6hour_f1)
  nc_TREFHT_6hour_f2 <- nc_open(TREFHT_6hour_f2)
  nc_TREFHT_6hour_f3 <- nc_open(TREFHT_6hour_f3)
  
  if(r == 1){
    lon <- ncvar_get(nc_PS_6hour_f1, "lon")
    lat <- ncvar_get(nc_PS_6hour_f1, "lat")
    lon_convert <- ifelse(lon <= 180, lon, lon - 360)
    min_lon_ind <- which.min(abs(lon_convert - -125))
    max_lon_ind <- which.min(abs(lon_convert - -66))
    min_lat_ind <- which.min(abs(lat - 25))
    max_lat_ind <- which.min(abs(lat - 50))
    
    lon_convert_subs <- lon_convert[min_lon_ind:max_lon_ind]
    lat_subs <- lat[min_lat_ind:max_lat_ind]
  }
  
  PS_6hour <- array(dim = c(48, 28, length(dates)*4))
  QREFHT_6hour <- array(dim = c(48, 28, length(dates)*4))
  TREFHT_6hour <- array(dim = c(48, 28, length(dates)*4))
  
  PS_6hour[,,year(dates_6hour) <= 2005] <- ncvar_get(nc_PS_6hour_f1, "PS",
                                                     start = c(min_lon_ind, min_lat_ind, 1),
                                                     count = c(max_lon_ind - min_lon_ind + 1, 
                                                               max_lat_ind - min_lat_ind + 1, 
                                                               sum(year(dates_6hour) <= 2005)))
  PS_6hour[,,year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026] <- ncvar_get(nc_PS_6hour_f2, "PS",
                                                                                 start = c(min_lon_ind, min_lat_ind, 1),
                                                                                 count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                           max_lat_ind - min_lat_ind + 1, 
                                                                  sum(year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026)))
  PS_6hour[,,year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071] <- ncvar_get(nc_PS_6hour_f3, "PS",
                                                                                 start = c(min_lon_ind, min_lat_ind, 1),
                                                                                 count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                           max_lat_ind - min_lat_ind + 1, 
                                                                  sum(year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071)))
  QREFHT_6hour[,,year(dates_6hour) <= 2005] <- ncvar_get(nc_QREFHT_6hour_f1, "QREFHT",
                                                         start = c(min_lon_ind, min_lat_ind, 1),
                                                         count = c(max_lon_ind - min_lon_ind + 1, 
                                                                   max_lat_ind - min_lat_ind + 1, 
                                                                  sum(year(dates_6hour) <= 2005)))
  QREFHT_6hour[,,year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026] <- ncvar_get(nc_QREFHT_6hour_f2, "QREFHT",
                                                                                     start = c(min_lon_ind, min_lat_ind, 1),
                                                                                     count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                               max_lat_ind - min_lat_ind + 1, 
                                                                    sum(year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026)))
  QREFHT_6hour[,,year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071] <- ncvar_get(nc_QREFHT_6hour_f3, "QREFHT",
                                                                                     start = c(min_lon_ind, min_lat_ind, 1),
                                                                                     count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                               max_lat_ind - min_lat_ind + 1, 
                                                                    sum(year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071)))
  
  TREFHT_6hour[,,year(dates_6hour) <= 2005] <- ncvar_get(nc_TREFHT_6hour_f1, "TREFHT",
                                                         start = c(min_lon_ind, min_lat_ind, 1),
                                                         count = c(max_lon_ind - min_lon_ind + 1, 
                                                                   max_lat_ind - min_lat_ind + 1, 
                                                                   sum(year(dates_6hour) <= 2005)))
  TREFHT_6hour[,,year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026] <- ncvar_get(nc_TREFHT_6hour_f2, "TREFHT",
                                                                                     start = c(min_lon_ind, min_lat_ind, 1),
                                                                                     count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                               max_lat_ind - min_lat_ind + 1, 
                                                                  sum(year(dates_6hour) <= 2035 & year(dates_6hour) >= 2026)))
  TREFHT_6hour[,,year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071] <- ncvar_get(nc_TREFHT_6hour_f3, "TREFHT",
                                                                                     start = c(min_lon_ind, min_lat_ind, 1),
                                                                                     count = c(max_lon_ind - min_lon_ind + 1, 
                                                                                               max_lat_ind - min_lat_ind + 1, 
                                                                  sum(year(dates_6hour) <= 2080 & year(dates_6hour) >= 2071)))
  # }
  
  DEWP_6hour <- Td(QREFHT_6hour, PS_6hour) + 273.15
  DEWP_mat <- array(data = DEWP_6hour, dim = c(48, 28, 4, length(dates)))
  DEWPREFHT[,,,r] <- (DEWP_mat[,,1,]+DEWP_mat[,,2,]+DEWP_mat[,,3,]+DEWP_mat[,,4,])/4
  
  TREFHT_mat <- array(data = TREFHT_6hour, dim = c(48, 28, 4, length(dates)))
  TREFHT[,,,r] <- (TREFHT_mat[,,1,]+TREFHT_mat[,,2,]+TREFHT_mat[,,3,]+TREFHT_mat[,,4,])/4
    
  print(r)
  flush.console()
}

DEWPREFHT[DEWPREFHT > TREFHT] <- TREFHT[DEWPREFHT > TREFHT]

save(list = c("lon_convert_subs", "lat_subs", "DEWPREFHT", "TREFHT"),
     file = "/glade/work/apoppick/LENS/DEWP_TREFHT_calculated6Hour_LENS_CONUS.RData")

q(save = "no")