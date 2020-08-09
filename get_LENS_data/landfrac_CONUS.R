rm(list = ls())
load('/research/Temp_Dewpt/LENS/landfrac.RData')
lon_convert <- ifelse(lon <= 180, lon, lon - 360)

min_lon_ind <- which.min(abs(lon_convert - -125))
max_lon_ind <- which.min(abs(lon_convert - -66))
min_lat_ind <- which.min(abs(lat - 25))
max_lat_ind <- which.min(abs(lat - 50))

landfrac_CONUS <- landfrac[min_lon_ind:max_lon_ind, min_lat_ind:max_lat_ind]
lon_convert_subs <- lon_convert[min_lon_ind:max_lon_ind]
lat_subs <- lat[min_lat_ind:max_lat_ind]

fields::image.plot(lon_convert_subs, lat_subs, landfrac_CONUS)

save(file = "/research/Temp_Dewpt/LENS/landfrac_CONUS.RData", list = c("landfrac_CONUS", "lon_convert_subs", "lat_subs"))