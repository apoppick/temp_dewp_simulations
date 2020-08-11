# This code produces the summary map figures shown in the paper and supplement. Note: I've hard-coded the colorbars to make the breaks correspond to nice numbers. This is dangerous and will need to be changed if the analysis changes.

#nohup R CMD BATCH map_figures.R map_figures.Rout &

rm(list=ls())
library(fields)
library(lubridate)
load('/research/Temp_Dewpt/Results/data_output/simulate_fromGSOD.RData')
load('/research/Temp_Dewpt/Results/data_output/estimate_CESM_models.RData')

load('/research/Temp_Dewpt/LENS/landfrac_CONUS.RData')

#### Extract and reformat results output for plots ####

nLon <- length(lon_convert_subs)
nLat <- length(lat_subs)

yearDates <- seq.Date(mdy("01-01-1920"), mdy("12-31-1920"), by = 1)
yearDates <- yearDates[!(month(yearDates) == 2 & day(yearDates) == 29)]
season_ind <- (month(yearDates) >= 6 &  month(yearDates) <= 8)

num_harmonics <- 2
X_c <- (cos((2*pi/365) * outer(0:(365 - 1), 1:num_harmonics)))[season_ind, ]
X_s <- (sin((2*pi/365) * outer(0:(365 - 1), 1:num_harmonics)))[season_ind, ]
extract_meanChangeCoef <- matrix(c(0, 1, rep(0, 4), X_c[1,], X_s[1,],
                                   0, 1, rep(0, 4), X_c[1+30,], X_s[1+30,],
                                   0, 1, rep(0, 4), X_c[1+30+31,], X_s[1+30+31,]),
                                 nrow = 10, ncol = 3)

# for temperature models
deltasToPlot <- array(dim = c(nLon, nLat, 47))
delta_SD <- array(dim = c(nLon, nLat,  47))
meanChangeCoefsToPlot <- array(dim = c(nLon, nLat,  3))
dev_stat <- array(dim = c(nLon, nLat))
R2 <- array(dim = c(nLon, nLat))
I1 <- array(dim = c(nLon, nLat))
I0 <- array(dim = c(nLon, nLat))

# for quantile regression models:
GMT_coef <- array(dim = c(nLon, nLat, 11))
onedegChange <- array(dim = c(nLon, nLat, 3*3, 11))
pearson_stat <- array(dim = c(nLon, nLat))

# for comparison of GSOD and CESM simulations
nStations <- length(GSOD_simul)
riskProbs_GSOD <- array(dim = c(nStations, 3, 4))
riskProbs_CESM <- array(dim = c(nLon, nLat, 3, 4))
riskProbs_CESM_closest <- array(dim = c(nStations, 3, 4))
sLat <- numeric(nStations)
sLon <- numeric(nStations)

for(l in 1:nLon){
  for(L in 1:nLat){
    ind <- (L-1)*nLon + l
    
    deltasToPlot[l,L,] <- CESM_models[[ind]]$delta_hat
    delta_SD[l,L,] <- sqrt(CESM_models[[ind]]$V_deltaSmooth)
    meanChangeCoefsToPlot[l,L,] <- t(CESM_models[[ind]]$meanTemp_coef) %*% extract_meanChangeCoef
    dev_stat[l, L] <- CESM_models[[ind]]$deviance
    R2[l, L] <- CESM_models[[ind]]$R2
    I0[l, L] <- CESM_models[[ind]]$I0
    I1[l, L] <- CESM_models[[ind]]$I1
    
    GMT_coef[l, L, ] <-  CESM_models[[ind]]$quantModelCESMCoef[2,]
    pearson_stat[l, L] <-  CESM_models[[ind]]$quantModel_pearson_stat
    riskProbs_CESM[l, L, ,] <- as.matrix(CESM_models[[ind]]$riskProbs) #month, Temp, Dewp high, Dewp low
  }
}

for(s in 1:nStations){
  sLon[s] <- GSOD_metadata$LON[s]
  sLat[s] <- GSOD_metadata$LAT[s]
  
  closest_ind <- GSOD_simul[[s]]$closest_ind
  riskProbs_GSOD[s,,] <- as.matrix(GSOD_simul[[s]]$riskProbs) #month, Temp, Dewp high, Dewp low
  riskProbs_CESM_closest[s,,] <- as.matrix(CESM_models[[closest_ind]]$riskProbs)
}

# p values for variability change terms:
p_delta <- 2*pnorm(abs(deltasToPlot/delta_SD), lower.tail = FALSE)
n_p <- prod(dim(p_delta))
p_FDR <- sort(p_delta)[max(which(sort(p_delta)  <= (1:n_p)/n_p * 0.05))]
isSig <- (p_delta <= p_FDR)



##### Figure 4 ######
## Temperature changes
pdf("/research/Temp_Dewpt/Results/figures/combined_temp_CONUS.pdf",
    height = 8/1.5, width = 14)

par(mar = c(2,2,2,2), oma = c(0,0,0,7))
# first divide screen into the figure region (left) and legend region (right)
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))

# now subdivide up the figure region into two parts
split.screen(c(2,3), screen = 1)-> ind

screen(ind[1])
image(lon_convert_subs, lat_subs, meanChangeCoefsToPlot[,,1]*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
      breaks = seq(0.8-0.05, 2.2+0.05, by = 0.1),
      col = designer.colors(15, c(start = "blue", middle = "white",end = "red"),
                            x = c(0,
                                  (1 - 0.8-0.1)/(2.2+0.1-0.8),
                                  1)),
      #col = designer.colors(10, c(start = "blue", middle = "white",end = "red"),
      #                      x = c(0, 
      #                            (1 - min(meanChangeCoefsToPlot))/(max(meanChangeCoefsToPlot) - min(meanChangeCoefsToPlot)),
      #                            1)),
      #zlim = c(min(meanChangeCoefsToPlot), max(meanChangeCoefsToPlot)),
      xlab = "", ylab = "", main = "July 1", axes = FALSE,
      cex.main = 1.5)
map("world", add = TRUE, fill = FALSE, col = "black")
map("state", add = TRUE, fill = FALSE, col = "black")

screen(ind[2])
image(lon_convert_subs, lat_subs, meanChangeCoefsToPlot[,,2]*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
      breaks = seq(0.8-0.05, 2.2+0.05, by = 0.1),
      col = designer.colors(15, c(start = "blue", middle = "white",end = "red"),
                            x = c(0,
                                  (1 - 0.8-0.1)/(2.2+0.1-0.8),
                                  1)),
      #col = designer.colors(10, c(start = "blue", middle = "white",end = "red"),
      #                      x = c(0, 
      #                            (1 - min(meanChangeCoefsToPlot))/(max(meanChangeCoefsToPlot) - min(meanChangeCoefsToPlot)),
      #                            1)),
      #zlim = c(min(meanChangeCoefsToPlot), max(meanChangeCoefsToPlot)),
      xlab = "", ylab = "", main = "July 1", axes = FALSE,
      cex.main = 1.5)
map("world", add = TRUE, fill = FALSE, col = "black")
map("state", add = TRUE, fill = FALSE, col = "black")

screen(ind[3])
image(lon_convert_subs, lat_subs, meanChangeCoefsToPlot[,,3]*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
      breaks = seq(0.8-0.05, 2.2+0.05, by = 0.1),
      col = designer.colors(15, c(start = "blue", middle = "white",end = "red"),
                            x = c(0,
                                  (1 - 0.8-0.1)/(2.2+0.1-0.8),
                                  1)),
      #col = designer.colors(10, c(start = "blue", middle = "white",end = "red"),
      #                      x = c(0, 
      #                            (1 - min(meanChangeCoefsToPlot))/(max(meanChangeCoefsToPlot) - min(meanChangeCoefsToPlot)),
      #                            1)),
      #zlim = c(min(meanChangeCoefsToPlot), max(meanChangeCoefsToPlot)),
      xlab = "", ylab = "", main = "August 1", axes = FALSE,
      cex.main = 1.5)
map("world", add = TRUE, fill = FALSE, col = "black")
map("state", add = TRUE, fill = FALSE, col = "black")


screen(ind[4])
image(lon_convert_subs, lat_subs, (exp(0.5*deltasToPlot[,,1])-1)*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
      breaks = seq(-0.15-0.05/2, 0.15+0.05/2, by = 0.05),
      col = two.colors(n=7, start = "blue", end = "red", middle = "white"),
      # col = two.colors(n=10, start = "blue", end = "red", middle = "white"),
      # zlim = c(-max(abs(exp(0.5*deltasToPlot)-1)), 
      #          max(abs(exp(0.5*deltasToPlot)-1))),
      xlab = "", ylab = "", main = "JJA Average", axes = FALSE,
      cex.main = 1.5)
map("world", add = TRUE, fill = FALSE, col = "black")
map("state", add = TRUE, fill = FALSE, col = "black")
whichIsNotSig_delta <- which(isSig[,,1] == 0 & landfrac_CONUS >= 0.5, arr.ind = TRUE)
text(lon_convert_subs[whichIsNotSig_delta[,1]], lat_subs[whichIsNotSig_delta[,2]], labels = '.')

screen(ind[5])
image(lon_convert_subs, lat_subs, (exp(0.5*deltasToPlot[,,2])-1)*ifelse(landfrac_CONUS >= 0.5, 1, NA),
      breaks = seq(-0.15-0.05/2, 0.15+0.05/2, by = 0.05),
      col = two.colors(n=7, start = "blue", end = "red", middle = "white"),
      # col = two.colors(n=10, start = "blue", end = "red", middle = "white"),
      # zlim = c(-max(abs(exp(0.5*deltasToPlot)-1)), 
      #          max(abs(exp(0.5*deltasToPlot)-1))),
      xlab = "", ylab = "", main = "3 Month Period", axes = FALSE,
      cex.main = 1.5)
map("world", add=T,fill=FALSE, col="black")
map("state", add = TRUE, fill = FALSE, col = "black")
whichIsNotSig_delta <- which(isSig[,,2] == 0 & landfrac_CONUS >= 0.5, arr.ind = TRUE)
text(lon_convert_subs[whichIsNotSig_delta[,1]], lat_subs[whichIsNotSig_delta[,2]], labels = '.')

screen(ind[6])
image(lon_convert_subs, lat_subs, (exp(0.5*deltasToPlot[,,47])-1)*ifelse(landfrac_CONUS >= 0.5, 1, NA),
      breaks = seq(-0.15-0.05/2, 0.15+0.05/2, by = 0.05),
      col = two.colors(n=7, start = "blue", end = "red", middle = "white"),
      # col = two.colors(n=10, start = "blue", end = "red", middle = "white"),
      # zlim = c(-max(abs(exp(0.5*deltasToPlot)-1)), 
      #          max(abs(exp(0.5*deltasToPlot)-1))),
      xlab = "", ylab = "", main = "2 Day Period", axes = FALSE,
      cex.main = 1.5)
map("world", add=T,fill=FALSE, col="black")
map("state", add = TRUE, fill = FALSE, col = "black")
whichIsNotSig_delta2 <- which(isSig[,,47] == 0  & landfrac_CONUS >= 0.5, arr.ind = TRUE)
text(lon_convert_subs[whichIsNotSig_delta2[,1]], lat_subs[whichIsNotSig_delta2[,2]], labels = '.')

# move to skinny region on right and draw the legend strip 
split.screen(c(2,1), screen = 2)-> ind

screen(ind[1])
image.plot(lon_convert_subs, lat_subs, meanChangeCoefsToPlot[,,3]*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
           legend.only=TRUE, 
           breaks = seq(0.8-0.05, 2.2+0.05, by = 0.1),
           col = designer.colors(15, c(start = "blue", middle = "white",end = "red"),
                                 x = c(0,
                                       (1 - 0.8-0.1)/(2.2+0.1-0.8),
                                       1)),
           # col = designer.colors(10, c(start = "blue", middle = "white",end = "red"),
           #                      x = c(0,
           #                            (1 - min(meanChangeCoefsToPlot))/(max(meanChangeCoefsToPlot) - min(meanChangeCoefsToPlot)),
           #                            1)),
           # zlim = c(min(meanChangeCoefsToPlot), max(meanChangeCoefsToPlot)),
           legend.lab = expression(paste("mean change /", degree, "C")),
           legend.line = 4.75, legend.cex = 1.4, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))

screen(ind[2])
image.plot(lon_convert_subs, lat_subs, exp(0.5*deltasToPlot[,,47])-1,
           legend.only=TRUE, 
           breaks = seq(-0.15-0.05/2, 0.15+0.05/2, by = 0.05),
           col = two.colors(n=7, start = "blue", end = "red", middle = "white"),
           # zlim = c(-max(abs(exp(0.5*deltasToPlot)-1)), 
           #          max(abs(exp(0.5*deltasToPlot)-1))),
           legend.lab = expression(paste("relative change in SD / ", degree, "C")),
           legend.line = 4.75, legend.cex = 1.4, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))
dev.off()



#### Figure S2 #####

# Temperature model validation
(worst_mean <- which.max(I1*ifelse(landfrac_CONUS == 1, 1, NA)))
(typical_mean <- which(I1*ifelse(landfrac_CONUS == 1, 1, NA) == median(I1*ifelse(landfrac_CONUS == 1, 
                                                                                 1, NA), na.rm = TRUE)))
(best_mean <- which.min(I1*ifelse(landfrac_CONUS == 1, 1, NA)))

(worst_dev <- which.max(dev_stat*ifelse(landfrac_CONUS == 1, 1, NA)))
(typical_dev <- which(dev_stat*ifelse(landfrac_CONUS == 1, 1, NA) == median(dev_stat*ifelse(landfrac_CONUS == 1, 
                                                                                 1, NA), na.rm = TRUE)))
(best_dev <- which.min(dev_stat*ifelse(landfrac_CONUS == 1, 1, NA)))

pdf('/research/Temp_Dewpt/Results/figures/validation/I1_and_deviance_maps.pdf', height = 5, width = 10)
par(mfrow = c(1,2), mar = c(2,2,2,2), oma = c(1.5,0,0,0))
image.plot(lon_convert_subs, lat_subs, I1*ifelse(landfrac_CONUS >= 0.5, 1, NA), 
           col = two.colors(n=9, start = "white", end = "black", middle = "gray"),
           breaks = seq(0.995-0.0025/2, 1.015+0.0025/2, by = 0.0025),
           legend.lab = expression(I[1]),
           legend.line = 2.5, legend.cex = 1.25, legend.width = 1, axis.args = list(cex.axis = 1.25),
           xlab = "", ylab = "", main = "Mean Model Validation", axes = FALSE,
           cex.main = 1.5, horizontal = TRUE)
map("world", add=T,fill=FALSE, col="black")
text(x=lon_convert_subs[arrayInd(worst_mean, dim(I1))[1]], 
     y=lat_subs[arrayInd(worst_mean, dim(I1))[2]], labels = "x", col = 'red', cex = 1.5)
text(x=lon_convert_subs[arrayInd(best_mean, dim(I1))[1]], 
     y=lat_subs[arrayInd(best_mean, dim(I1))[2]], labels = "o", col = 'blue', cex = 1.5)
points(x=lon_convert_subs[arrayInd(typical_mean, dim(I1))[1]], 
     y=lat_subs[arrayInd(typical_mean, dim(I1))[2]], pch = 2, col = 'green', cex = 1.5, lwd = 2)


image.plot(lon_convert_subs, lat_subs, dev_stat*ifelse(landfrac_CONUS >= 0.5, 1, NA),
           col = two.colors(n=9, start = "white", end = "black", middle = "gray"),
           breaks = seq(1400-50, 2200+50, by = 100),
           legend.lab = "Deviance",
           legend.line = 2.5, legend.cex = 1.25, legend.width = 1, axis.args = list(cex.axis = 1.25),
           xlab = "", ylab = "", main = "Variability Model Validation", axes = FALSE,
           cex.main = 1.5, horizontal = TRUE)
map("world", add=T,fill=FALSE, col="black")
text(x=lon_convert_subs[arrayInd(worst_dev, dim(I1))[1]], 
     y=lat_subs[arrayInd(worst_dev, dim(I1))[2]], labels = "x", col = 'red', cex = 1.5)
text(x=lon_convert_subs[arrayInd(best_dev, dim(I1))[1]], 
     y=lat_subs[arrayInd(best_dev, dim(I1))[2]], labels = "o", col = 'blue', cex = 1.5)
points(x=lon_convert_subs[arrayInd(typical_dev, dim(I1))[1]], 
     y=lat_subs[arrayInd(typical_dev, dim(I1))[2]], pch = 2, col = 'green', cex = 1.5, lwd  =2)
dev.off()




##### Figure 5 ######

# GMT coefficients for quantile regression model
pdf("/research/Temp_Dewpt/Results/figures/quantModel_GMTCoefs.pdf",height = 8/3, width = 14)
par(mar = c(2,2,2,2), oma = c(0,0,0,5))
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))

# now subdivide up the figure region into two parts
split.screen(c(1,3), screen = 1)-> ind

screen(ind[1])
image(x = lon_convert_subs, y = lat_subs, z= ifelse(landfrac_CONUS >= 0.5, 1, NA)*GMT_coef[,,3]/log(10),
      breaks = seq(-0.15, 0.15, by = 0.02),
      col = two.colors(n=15, start = "cyan4", end = "coral4", middle = "white"),
      # zlim = c(-max(abs(GMT_coef/log(10))), max(abs(GMT_coef/log(10)))), 
      # col = two.colors(n=10, start = "cyan4", end = "coral4", middle = "white"),
      xlab = "", ylab = "", main = "5th percentile", axes = FALSE,
      cex.main = 2)
map("world", add=T,fill=FALSE, col="black")
map("state", add=T,fill=FALSE, col="black")

screen(ind[2])
image(x = lon_convert_subs, y = lat_subs, z= ifelse(landfrac_CONUS >= 0.5, 1, NA)*GMT_coef[,,6]/log(10),
      breaks = seq(-0.15, 0.15, by = 0.02),
      col = two.colors(n=15, start = "cyan4", end = "coral4", middle = "white"),
      # zlim = c(-max(abs(GMT_coef/log(10))), max(abs(GMT_coef/log(10)))), 
      # col = two.colors(n=10, start = "cyan4", end = "coral4", middle = "white"),
      xlab = "", ylab = "", main = "50th percentile", axes = FALSE,
      cex.main = 2)
map("world", add=T,fill=FALSE, col="black")
map("state", add=T,fill=FALSE, col="black")

screen(ind[3])
image(x = lon_convert_subs, y = lat_subs, z= ifelse(landfrac_CONUS >= 0.5, 1, NA)*GMT_coef[,,9]/log(10),
      breaks = seq(-0.15, 0.15, by = 0.02),
      col = two.colors(n=15, start = "cyan4", end = "coral4", middle = "white"),
      # zlim = c(-max(abs(GMT_coef/log(10))), max(abs(GMT_coef/log(10)))), 
      # col = two.colors(n=10, start = "cyan4", end = "coral4", middle = "white"),
      xlab = "", ylab = "", main = "95th percentile", axes = FALSE, cex.main = 2)
map("world", add=T,fill=FALSE, col="black")
map("state", add=T,fill=FALSE, col="black")

# move to skinny region on right and draw the legend strip 
split.screen(c(1,1), screen = 2)-> ind

screen(ind[1])
image.plot(x = lon_convert_subs, y = lat_subs, 
           z= ifelse(landfrac_CONUS >= 0.5, 1, NA)*GMT_coef[,,3]/log(10),
           legend.only=TRUE, 
            breaks = seq(-0.15, 0.15, by = 0.02),
            col = two.colors(n=15, start = "cyan4", end = "coral4", middle = "white"),
            # zlim = c(-max(abs(GMT_coef/log(10))), max(abs(GMT_coef/log(10)))), 
            # col = two.colors(n=10, start = "cyan4", end = "coral4", middle = "white"),
            legend.lab = expression(tilde(alpha)[paste("1,",tau)]),
            legend.line = 4.75, legend.cex = 1.5, legend.width = 3,
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))
dev.off()


##### Figure S9 #####

## Pearson statistics
(worst_pearson_stat <- which.max(apply(pearson_stat,c(1,2), sum)*ifelse(landfrac_CONUS >= 0.5, 1, NA)))
(best_pearson_stat <- which.min(apply(pearson_stat,c(1,2), sum)*ifelse(landfrac_CONUS >= 0.5, 1, NA)))
(typical_pearson_stat <- 896) # the quantile regression fit is slightly random causing the location of the median Pearson statistic to change each time. Not ideal to hard code, but just an example of a typical fit.
  # which(apply(pearson_stat,c(1,2), sum)*ifelse(landfrac_CONUS >= 0.5, 1, NA) ==
  #                     median(apply(pearson_stat,c(1,2), sum)*ifelse(landfrac_CONUS >= 0.5, 1, NA),
  #                            na.rm = TRUE)))

pdf('/research/Temp_Dewpt/Results/figures/validation/quantModel_pearsonStat.pdf', height = 7, width = 11)
par(mfrow = c(1,1), mar = c(2,2,2,2), oma = c(1.5,0,0,1.5))
image.plot(x = lon_convert_subs, y = lat_subs, 
           z= ifelse(landfrac_CONUS >= 0.5, 1, NA)*pearson_stat,
           breaks = seq(0-5000/2, 35000+5000/2, by = 5000),
           col = two.colors(n=8, start = "white", end = "black", middle = "gray"),
           ylab = "", xlab = "", main = "Quantile Regression Model Validation", axes = FALSE,
           cex.main = 1.5, legend.lab = "Pearson statistic", legend.line = 3.5, 
           legend.cex = 1, legend.width = 1, axis.args = list(cex.axis = 1))
map("world", add=T,fill=FALSE, col="black")
map("state", add=T,fill=FALSE, col="black")
text(x=lon_convert_subs[arrayInd(worst_pearson_stat, dim(pearson_stat))[1]], 
     y=lat_subs[arrayInd(worst_pearson_stat, dim(pearson_stat))[2]], labels = "x", col = 'red', cex = 1.5)
text(x=lon_convert_subs[arrayInd(best_pearson_stat, dim(pearson_stat))[1]], 
     y=lat_subs[arrayInd(best_pearson_stat, dim(pearson_stat))[2]], labels = "o", col = 'blue', cex = 1.5)
points(x=lon_convert_subs[arrayInd(typical_pearson_stat, dim(pearson_stat))[1]], 
     y=lat_subs[arrayInd(typical_pearson_stat, dim(pearson_stat))[2]], 
     pch = 2, col = 'green', lwd = 2, cex = 1.5)
dev.off()


##### Figure 6 #####

## Comparison of changes in risk probabilities

logORs <- log10(riskProbs_GSOD[,,2:4]/(1- riskProbs_GSOD[,,2:4]))  -
  log10(riskProbs_CESM_closest[,,2:4]/(1- riskProbs_CESM_closest[,,2:4])) # temperature, dewp high, dewp low

# note: all the riskProbs are quantile levels, so for the "high" plot we are interested in 1-riskProb (i.e., probabiltiy of exceedance)

month_names <- c("June", "July", "August")

pdf("/research/Temp_Dewpt/Results/figures/dewp_risks_high.pdf",height = 8, width = 14)
par(mar = c(2,2,2,2), oma = c(0,0,0,4))
# first divide screen into the figure region (left) and legend region (right)
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))

# now subdivide up the figure region into two parts
split.screen(c(3,3), screen = 1)-> ind

k <- 1

for(m in 1:3){
  screen(ind[k])
  image(lon_convert_subs, lat_subs, ifelse(landfrac_CONUS >= 0.5, 1, NA)*(-log10(riskProbs_CESM[,,m,3] / (1-riskProbs_CESM[,,m,3])) - log10(0.05 / 0.95)),
        col = two.colors(n=9, start = "coral4", end = "cyan4", middle = "white"),
        breaks = seq(-4.5, 4.5, by=1),
        #zlim = c(-4, 4),
        xlab = "", ylab = "", main = "", axes = FALSE,
        cex.main = 1.5)
  map("world", add=T,fill=FALSE, col="black")
  map("state", add = TRUE, fill = FALSE, col = "black")
  mtext(month_names[m], cex = 1.5)
  if(m == 1){
    mtext("CESM", side = 2, cex = 1.5)
  }
  ind <- ind + 1
}


for(m in 1:3){
  screen(ind[k])
  usa_map <- map("usa", plot = FALSE)
  plot(0,0, type='n', xlim = c(-125, -66), ylim = c(25, 50), axes = FALSE,
       ylab = "", xlab = "", main = "")
  polygon(usa_map$x, usa_map$y, col = "gray")
  map("state", add=T,fill=FALSE, col="black")
  points(x = sLon, y = sLat, col = color.scale(-log10(riskProbs_GSOD[,m,3]/(1-riskProbs_GSOD[,m,3])) -
                                                 log10(0.05/0.95), 
                                               col = two.colors(n=9, start = "coral4", end = "cyan4", middle = "white"),
                                               zlim = c(-4, 4), eps = 0,
                                               transparent.color = "white"),
         pch = 16, cex = 1.5)
  
  k <- k+1
  if(m == 1){
    mtext("Simulation", side = 2, cex = 1.5)
  }
}

for(m in 1:3){
  screen(ind[k])
  usa_map <- map("usa", plot = FALSE)
  plot(0,0, type='n', xlim = c(-125, -66), ylim = c(25, 50), axes = FALSE,
       ylab = "", xlab = "", main = "")
  polygon(usa_map$x, usa_map$y, col = "gray")
  map("state", add=T,fill=FALSE, col="black")
  points(x = sLon, y = sLat, col = color.scale(-logORs[,m,2], 
                                               col = two.colors(n = 9, start = "orange", end = "purple", 
                                                                middle = "white"),
                                               zlim = c(-2,2), eps = 0,
                                               transparent.color = "white"),
         pch = 16, cex = 1.5)
  
  if(m == 1){
    mtext("Simulation vs. CESM", side = 2, cex = 1.5)
  }
  k <- k+1
}

split.screen(c(3,1), screen = 2)-> ind

screen(ind[1])
image.plot(lon_convert_subs, lat_subs, ifelse(landfrac_CONUS >= 0.5, 1, NA)*(-log10(riskProbs_CESM[,,m,3] / (1-riskProbs_CESM[,,m,3])) - log10(0.05 / 0.95)),
           legend.only=TRUE, 
           col = two.colors(n=9, start = "coral4", end = "cyan4", middle = "white"),
           breaks = seq(-4.5, 4.5, by=1),
           #zlim = c(-4, 4),
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))

screen(ind[2])
image.plot(lon_convert_subs, lat_subs, ifelse(landfrac_CONUS >= 0.5, 1, NA)*(-log10(riskProbs_CESM[,,m,3] / (1-riskProbs_CESM[,,m,3])) - log10(0.05 / 0.95)),
           legend.only=TRUE, 
           col = two.colors(n=9, start = "coral4", end = "cyan4", middle = "white"),
           zlim = c(-4, 4), eps = 0,
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))

screen(ind[3])
image.plot(lon_convert_subs, lat_subs, ifelse(landfrac_CONUS >= 0.5, 1, NA)*(-log10(riskProbs_CESM[,,m,3] / (1-riskProbs_CESM[,,m,3])) - log10(0.05 / 0.95)),
           legend.only=TRUE, 
           col = two.colors(n = 9, start = "orange", end = "purple", 
                            middle = "white"),
           zlim = c(-2,2), eps = 0,
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))
close.screen(all = TRUE)
dev.off()



##### Figure 7 #####

pdf("/research/Temp_Dewpt/Results/figures/dewp_risks_low.pdf",height = 8, width = 14)
par(mar = c(2,2,2,2), oma = c(0,0,0,4))
# first divide screen into the figure region (left) and legend region (right)
split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))

# now subdivide up the figure region into two parts
split.screen(c(3,3), screen = 1)-> ind

k <- 1

for(m in 1:3){
  screen(ind[k])
  image(lon_convert_subs, lat_subs, ifelse(landfrac_CONUS >= 0.5, 1, NA)*(log10(riskProbs_CESM[,,m,4] / (1-riskProbs_CESM[,,m,4])) - log10(0.05 / 0.95)),
        col = two.colors(n=9, start = "cyan4", end = "coral4", middle = "white"),
        zlim = c(-2, 2),
        xlab = "", ylab = "", main = "", axes = FALSE,
        cex.main = 1.5)
  map("world", add=T,fill=FALSE, col="black")
  map("state", add = TRUE, fill = FALSE, col = "black")
  
  mtext(month_names[m], cex = 1.5)
  if(m == 1){
    mtext("CESM", side = 2, cex = 1.5)
  }
  ind <- ind + 1
}


for(m in 1:3){
  screen(ind[k])
  usa_map <- map("usa", plot = FALSE)
  plot(0,0, type='n', xlim = c(-125, -66), ylim = c(25, 50), axes = FALSE,
       ylab = "", xlab = "", main = "")
  polygon(usa_map$x, usa_map$y, col = "gray")
  map("state", add=T,fill=FALSE, col="black")
  points(x = sLon, y = sLat, col = color.scale(log10(riskProbs_GSOD[,m,4] / (1 - riskProbs_GSOD[,m,4])) - 
                                                 log10(0.05 / 0.95), 
                                               col = two.colors(n = 9, start = "cyan4", end = "coral4", 
                                                                middle = "white"),
                                               zlim = c(-2,2),
                                               transparent.color = "white", eps = 1e-8),
         pch = 16, cex = 1.5)
  
  k <- k+1
  
  if(m == 1){
    mtext("Simulation", side = 2, cex = 1.5)
  }
  
}

for(m in 1:3){
  screen(ind[k])
  usa_map <- map("usa", plot = FALSE)
  plot(0,0, type='n', xlim = c(-125, -66), ylim = c(25, 50), axes = FALSE,
       ylab = "", xlab = "", main = "")
  polygon(usa_map$x, usa_map$y, col = "gray")
  map("state", add=T,fill=FALSE, col="black")
  points(x = sLon, y = sLat, col = color.scale(logORs[,m,3], 
                                               col = two.colors(n = 9, start = "orange", end = "purple", 
                                                                middle = "white"),
                                               zlim = c(-2,2),
                                               transparent.color = "white", eps = 1e-8),
         pch = 16, cex = 1.5)
  
  if(m == 1){
    mtext("Simulation vs. CESM", side = 2, cex = 1.5)
  }
  k <- k+1
}

split.screen(c(3,1), screen = 2)-> ind

screen(ind[1])
image.plot(legend.only=TRUE, 
           zlim = c(-2,2),
           col = two.colors(n = 9, start = "cyan4", end = "coral4", 
                            middle = "white"),
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))

screen(ind[2])
image.plot(legend.only=TRUE, 
           zlim = c(-2,2),
           col = two.colors(n = 9, start = "cyan4", end = "coral4", 
                            middle = "white"),
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))

screen(ind[3])
image.plot(legend.only=TRUE, 
           zlim = c(-2,2),
           col = two.colors(n = 9, start = "orange", end = "purple", 
                            middle = "white"),
           legend.lab = "log odds ratio",
           legend.line = 4, legend.cex = 1.5, legend.width = 3, 
           axis.args = list(cex.axis = 1.5),
           smallplot=c(.1,.2, .1,.9))
close.screen(all = TRUE)
dev.off()
