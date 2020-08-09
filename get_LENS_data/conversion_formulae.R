## See Willett et al. (2014), Table 1

e_fun <- function(q, P){
  # q is specific humidity in kg/kg
  # P is pressure in Pa
  
  P <- P/100 #convert to hPa
  q <- q*1000 #convert to g/kg
  
  q*P/(1000*(0.622 + q/1000*(1-0.622)))
}

fw_fun <- function(P){
  # P is pressure in Pa
  P <- P/100 #convert to hPa
  1 + 7e-4 + 3.46e-6 * P
}

Td <- function(q, P){
  # q is specific humidity in kg/kg
  # P is pressure in Pa
  
  e <- e_fun(q, P)
  fw <- fw_fun(P)
  a <- -1/227.3
  b <- (18.729 - log(e/(6.1121*fw)))
  c <- -257.87 * log(e/(6.1121*fw))
  
  (-b + sqrt(b^2 - 4*a*c)) / (2*a)
}

PS <- function(PSL, Z, Temp){
  # Z is height in meters
  # Temp is in Kelvin
  # PSL is in Pa
  
  PSL <- PSL/100 #convert to hPa
  
  PSL * (Temp / (Temp + 0.0065*Z))^5.625
  
}

es_fun <- function(P, Temp){
  # P is pressure in Pa
  # Temp is temperature in K

  Temp <- Temp - 273.15 # convert to C

  fw <- fw_fun(P)
  6.1121 * fw * exp((18.729 - Temp / 227.3)*Temp /
                      (257 + Temp))
}

RH_fun <- function(q, P, Temp){
  # P is pressure in Pa
  # Temp is temperature in K
  # q is specific humidity in kg/kg

  e <- e_fun(q, P)
  es <- es_fun(P, Temp)
  100 * e / es
}

Td2 <- function(RH, Temp){
  # RH is a percentage (0 to 100)
  # Temp is in K
  
  b <- 17.27
  c <- 237.7
  
  Temp <- Temp - 273.15 #convert to C
  g <- log(RH/100) + b*Temp/(c + Temp)
  c * g / (b - g)
}