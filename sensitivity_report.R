########################### DAY 2 ###########################
library(deSolve)
library(fields)

#1. I0, average surface light
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(350/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  

  
  for (i in 1:(length(x))){ 
    param$I0 <- x[i] # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
    
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

I0 <- c(87.5,175,350,700,1400)
run <- sensitivity(I0)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2.6,30,legend=c("I0 = 87.5 uE", "I0 = 175 uE", "I0 = 350 uE","I0 = 700 uE","I0 = 1400 uE"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#2. kw, light attenuation of sea water
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kc <- 0.05 #light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$kw <- x[i] # light absorbed by water (1/m)
    
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

kw <- c(0.0375/4,0.0375/2,0.0375,0.0375*2,0.0375*4)
run <- sensitivity(kw)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(1.3,20,legend=c("kw = 0.0094 m^2/mmol N", "kw = 0.0188 m^2/mmol N", "kw = 0.0375 m^2/mmol N","kw = 0.0750 m^2/mmol N","kw = 0.1500 m^2/mmol N"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#3. kc, light attenuation of phytoplankton
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  

  
  for (i in 1:(length(x))){ 
    param$kc <- x[i] #light absorbed by phyto plankton (m^2/mmol N)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

kc <- c(0.0125,0.025,0.05,0.1,0.2)
run <- sensitivity(kc)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(1.3,30,legend=c("kc = 0.0125 m^2/mmol N", "kc = 0.025 m^2/mmol N", "kc = 0.05 m^2/mmol N","kc = 0.1 m^2/mmol N","kc = 0.2 m^2/mmol N"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#4. gmax, maximum growth rate 
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$gmax <- x[i] # maximum growth (1/d)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

gmax <- c(0.375,0.75,1.5,3,6)
run <- sensitivity(gmax)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2,30,legend=c("gmax = 0.375 1/d", "gmax = 0.75 1/d", "gmax = 1.5 1/d","gmax = 3 1/d","gmax = 6 1/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#5. alpha, initial slope of the P-I curve
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$alpha <- x[i] # Slope of the PI-curve (1/(uEinstein*d)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

alpha <- c(0.025/(350/200),0.05/(350/200),0.1/(350/200),0.2/(350/200),0.4/(350/200))
run <- sensitivity(alpha)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(1.75,30,legend=c("alpha = 0.014 1/uE/d", "alpha = 0.029 1/uE/d", "alpha = 0.057 1/uE/d","alpha = 0.114 1/uE/d","alpha = 0.229 1/uE/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#6. HN, half saturation for nutrients
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$HN <- x[i] # Half saturation constant for nutrients (mmol N/m^3)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

HN <- c(0.075,0.15,0.3,0.6,1.2)
run <- sensitivity(HN)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(1.3,30,legend=c("HN = 0.075 mmol N/m^3", "HN = 0.15 mmol N/m^3", "HN = 0.3 mmol N/m^3","HN = 0.6 mmol N/m^3","HN = 1.2 mmol N/m^3"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#7. m, background mortality
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$m <- x[i] # specific loss rate (1/d)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

m <- c(0.0075,0.015,0.03,0.06,0.12)
run <- sensitivity(m)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2.3,30,legend=c("m = 0.0075 1/d", "m = 0.015 1/d", "m = 0.03 1/d","m = 0.06 1/d","m = 1.2 1/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#8. gamma, grazing parameter
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$gamma <- x[i] # Grazing parameter (m^3/(mmol N * d))
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

gamma <- c(0.025,0.05,0.1,0.2,0.4)
run <- sensitivity(gamma)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(0.55,12,legend=c("gamma = 0.025 m^3/mmol N/d", "gamma = 0.05 m^3/mmol N/d", "gamma = 0.1 m^3/mmol N/d","gamma = 0.2 m^3/mmol N/d","gamma = 0.4 m^3/mmol N/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#9. epsilon, remineralization
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$eps <- x[i] # remineralization (1/d)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

eps <- c(0.025,0.05,0.1,0.2,0.4)
run <- sensitivity(eps)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2.3,30,legend=c("eps = 0.025 1/d", "eps = 0.05 1/d", "eps = 0.1 1/d","eps = 0.2 1/d","eps = 0.4 1/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#10. Dv, diffusivity
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$Dv <- x[i] # Diffusivity (m^2/d)  
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

Dv <- c(1.25/(10**(5))*3600*24,2.5/(10**(5))*3600*24,5/(10**(5))*3600*24,10/(10**(5))*3600*24,20/(10**(5))*3600*24)
run <- sensitivity(Dv)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2.1,30,legend=c("Dv = 1.08 m^2/d", "Dv = 2.16 m^2/d", "Dv = 4.32 m^2/d","Dv = 8.64 m^2/d","Dv = 17.28 m^2/d"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#11. uD, settling velocity of detritus
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$up <- 0.5 # Settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$ud <- x[i] # detritus settling velocity (m/d)
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

ud <- c(3.75,7.5,15,30,60)
run <- sensitivity(ud)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2,42,legend=c("uD = 3.75 m/day", "uD = 7.5 m/day", "uD = 15 m/day","uD = 30 m/day","uD = 60 m/day"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#12. up, settling velocity of phytoplankton
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$up <- x[i] # phytoplankton settling velocity (m/d)
    
 #   param$ud <- d_sink[i] # Settling velocity (m/d)
    
    
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

up <- c(0.125,0.25,0.5,1,2)
run <- sensitivity(up)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(2,42,legend=c("up = 0.125 m/day", "up = 0.25 m/day", "up = 0.5 m/day","up = 1.0 m/day","up = 2.0 m/day"),
        lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.85)
#####

#13. Nb, nutrient bottom concentration
#####
sensitivity <- function(x){
  # Create an empty list to add the parameters
  
  # Create an empty list to add the parameters
  param <- list()
  
  # Define parameters 
  param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
  param$kw <- 0.0375 # light absorbed by water (1/m)
  param$kc <- 0.05 # light absorbed by phyto plankton (m^2/mmol N)
  param$gmax <- 1.5 # maximum growth (1/d)
  param$alpha <- 0.1/(param$I0/200) # Slope of the PI-curve (1/(uEinstein*d)
  param$HN <- 0.3 # Half saturation constant for nutrients (mmol N/m^3)
  param$m <- 0.03 # specific loss rate (1/d)
  param$gamma <- 0.1 # Grazing parameter (m^3/(mmol N * d))
  param$eps <- 0.1 # remineralization (1/d)
  param$Dv <- 5/(10**(5))*3600*24 # Diffusivity (m^2/d)  
  param$ud <- 15 # detritus settling velocity (m/d)
  param$up <- 0.5 # Settling velocity (m/d)
  param$dz <- 2 # grid spacing (m)
  param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
  param$n <- length(param$z) # number of grid cells
  
  
  P_res <- list()
  N_res <- list()
  D_res <- list()
  
  
  
  for (i in 1:(length(x))){ 
    param$Nb <- x[i] # nutrient content at the bottom (mmol N/m^3) 
    
    #   param$ud <- d_sink[i] # Settling velocity (m/d)
    
    
    
    # Define initial conditions - concentrations of phytoplankton and nutrients
    PND <- c(rep(0.01,param$n),rep(param$Nb/10, param$n),rep(0,param$n))
    
    # Create a function
    CalLight <- function(t,P,D, param) {
      season <- 1+cos(pi+2*pi/365*t+4*pi/73)
      
      Q <- param$kc * param$dz * ((cumsum(P) - P/2) +  (cumsum(D) - D/2))
      
      I <- param$I0 * exp(-param$kw * param$z - Q)*season
      
      return(I)
    }
    
    
    
    
    # Create function that creates the differential equation at each step
    FuncPND <- function(t, PND, param) {
      P <- PND[1:param$n]
      N <- PND[(1+param$n):(2*param$n)]
      D <- PND[(2*param$n+1):(3*param$n)]
      
      # Phytoplankton flux
      Jap <- rep(0,param$n+1)
      Jdp <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jap[i] <- param$up * P[i-1]
        Jdp[i] <- -param$Dv * (P[i] - P[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jap[1] = 0
      Jap[param$n+1] = 0
      
      # Diffusive flux boundary
      Jdp[1] = 0
      Jdp[param$n+1] = 0
      
      Jp = Jap + Jdp
      
      # Nutrient flux
      Jdn <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jdn[i] <- -param$Dv * (N[i] - N[i-1]) / param$dz
      }
      
      
      # Diffusive flux boundary
      Jdn[1] = 0
      Jdn[param$n+1] = -param$Dv*(param$Nb-N[param$n])/param$dz
      
      Jn = Jdn
      
      # Detritus flux
      Jad <- rep(0,param$n+1)
      Jdd <- rep(0,param$n+1)
      
      for (i in 2:param$n){
        Jad[i] <- param$ud * D[i-1]
        Jdd[i] <- -param$Dv * (D[i] - D[i-1]) / param$dz
      }
      
      # Advective flux boundary
      Jad[1] = 0
      Jad[param$n+1] = param$ud*D[param$n]
      
      # Diffusive flux boundary
      Jdd[1] = 0
      Jdd[param$n+1] = 0
      
      Jd = Jad + Jdd
      
      
      I <- CalLight(t,P,D,param)
      
      g <- param$gmax * pmin(param$alpha*I / sqrt((param$alpha*I)**2 + param$gmax**2),N/(N+param$HN)) 
      
      dPdT <- -(Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g- param$m-param$gamma*P) * P
      
      dNdT <- -g*P +param$eps*D-(Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
      
      dDdT <- -(Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + (param$m+param$gamma*P) * P-param$eps*D
      
      
      return(list(c(dPdT,dNdT, dDdT)))
    }
    
    
    
    
    # Define time step
    time <- seq(1,3*365, by=1)
    
    # Solve our differential equation
    res <- ode(PND, time, FuncPND, param)
    
    P_res[[i]] <- res[,2:(param$n+1)]
    N_res[[i]] <- res[,(param$n+2):(param$n*2+1)]
    D_res[[i]] <- res[,(param$n*2+2):(param$n*3+1)]
    
    
  }
  return(list(P_res,N_res,D_res))
}

Nb <- c(7.5,15,30,60,120)
run <- sensitivity(Nb)  

P <- run[[1]]
N <- run[[2]]
D <- run[[3]]


mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
par(mfrow=c(3,4))
plot_nr <- c(1,4,7,10,2,5,8,11,3,6,9,12)
month <- c("Dec","Jan","Feb", "Mar","Apr", "May","Jun","Jul","Aug","Sep","Oct", "Nov")
#paste(as.character(month[index]))
title <- c("Winter","Spring","Summer","Fall","","","","","","","","")
xtitle <- c("","","","","","","","","mmol N/m^3","mmol N/m^3","mmol N/m^3","mmol N/m^3")
xax <- c("n","n","n","n","n","n","n","n","s","s","s","s")
ytitle <- c("Depth [m]","","","","Depth [m]","","","","Depth [m]","","","")
yax <- c("s","n","n","n","s","n","n","n","s","n","n","n")
l_marg <- c(4,1.5,1.5,1.5,4,1.5,1.5,1.5,4,1.5,1.5,1.5)
r_marg <- c(0.5,2,2,2,0.5,2,2,2,0.5,2,2,2)
b_marg <- c(1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,4,4,4,4)
t_marg <- c(2,2,2,2,1,1,1,1,0.1,0.1,0.1,0.1)
depth <- (seq(2/2,100,by = 2))

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]
  
  
  plot(P[[1]][mid_of_month[index],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3, 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 1)
  lines(P[[2]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
  lines(P[[3]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
  lines(P[[4]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
  lines(P[[5]][mid_of_month[index],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
  text(2.5,5,month[index])
}
legend(1.6,55,legend=c("Nb = 7.5 mmol N/m^3", "Nb = 15 mmol N/m^3", "Nb = 30 mmol N/m^3","Nb = 60 mmol N/m^3","Nb = 120 mmol N/m^3"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 0.8)
#####


# Plotting sensitivity test in mid march

par(mfrow=c(1,2))
mid_of_month <- 2*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
depth <- (seq(2/2,100,by = 2))

par(mar=c(4,4,1,0))

plot(P[[1]][mid_of_month[4],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3,
     lwd = 1,xlab = "mmol N/m^3", ylab = "Depth  [m]")
lines(P[[2]][mid_of_month[4],], depth, ylim = rev(range(depth)),lty = 3,lwd = 1.5,type = "l")
lines(P[[3]][mid_of_month[4],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2,type = "l")
lines(P[[4]][mid_of_month[4],], depth, ylim = rev(range(depth)),lty = 3,lwd = 2.5,type = "l")
lines(P[[5]][mid_of_month[4],], depth, ylim = rev(range(depth)),lty = 3,lwd = 3,type = "l")
text(2.5,5,"Nb")

par(mar=c(4,2.3,1,1))

plot(P[[1]][mid_of_month[4],], depth, xlim = c(0,5),ylim = rev(range(depth)), type = "l",lty = 3,
     lwd = 1,xlab = "mmol N/m^3",yaxt="n",col = "white")
legend(0.5,18,legend=c("Table value/4", "Table value/2", "Table value","Table value x 2","Table value x 4"),
       lty = 3,lwd = c(1,1.5,2,2.5,3),cex = 1)



