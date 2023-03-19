########################### DAY 3 ###########################
library(deSolve)
library(fields)

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
param$ud <- 15 # Settling velocity (m/d)
param$up <- 0.5 # Settling velocity (m/d)
param$Nb <- 30 # nutrient content at the bottom (mmol N/m^3) 
param$dz <- 1 # grid spacing (m)
param$z <- seq(param$dz/2,100,by = param$dz) # depth (m)
param$n <- length(param$z) # number of grid cells



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
time <- seq(1,365*4, by=1)

# Solve our differential equation
start.time <- Sys.time()
res <- ode(PND, time, FuncPND, param)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Runtime =",as.character(round(time.taken,1)),"s"))


P <- res[,2:(param$n+1)]
N <- res[,(param$n+2):(param$n*2+1)]
D <- res[,(param$n*2+2):(param$n*3+1)]

first_of_month <- 3*365+c(1,32,60,91,121,152,182,213,244,274,305,335)-31+15 # Day number at the 15th of each month in the 4th year of simulation. Starting from december in year 3
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

for (i in 1:12){  
  par(mar=c(b_marg[i],l_marg[i],t_marg[i],r_marg[i]))
  index <- plot_nr[i]

  
  plot(P[first_of_month[index],], param$z, xlim = c(0,3),ylim = rev(range(param$z)), type = "l",pch = 19, col = "darkseagreen", 
       main = paste(as.character(title[i])), ylab = paste(as.character(ytitle[i])),yaxt=paste(as.character(yax[i])), 
       xaxt=paste(as.character(xax[i])),xlab = paste(as.character(xtitle[i])),lwd = 3)
  lines(N[first_of_month[index],]/10, param$z, ylim = rev(range(param$z)),col = "blue3",lwd = 3,type = "l")
  lines(D[first_of_month[index],], param$z, ylim = rev(range(param$z)),col = "brown4",lwd = 3,type = "l")
  text(1.5,5,month[index])
}
legend(1.2,42,legend=c("Phytoplankton", expression("Nutrients x " ~ 10^{-1} ), "Detritus"),col=c("darkseagreen", "blue3", "brown4"), pch = 19,cex = 1)

