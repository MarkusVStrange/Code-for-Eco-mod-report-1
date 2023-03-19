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
time <- seq(1,365*5, by=1)

# Solve our differential equation
start.time <- Sys.time()
res <- ode(PND, time, FuncPND, param)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Runtime =",as.character(round(time.taken,1)),"s"))
P <- res[,2:(param$n+1)]
N <- res[,(param$n+2):(param$n*2+1)]
D <- res[,(param$n*2+2):(param$n*3+1)]

cm <- rep(0,365)
dm <- rep(0,365)
count <- 1
for (i in (1+4*365):(5*365)){
  cm[count] <- max(P[i,])
  dm[count] <- match(cm[count],P[i,])
  count <- count+1
}
par(mfrow=c(2,1))
par(mar=c(2,4,3,1))
plot(seq(1,365),cm, pch = 19,cex = 0.5, col = "darkseagreen",xaxt="n",xlab = "",
      ylab = substitute(paste(bold("mmol N/m^3"))), main = "Size of chlorophyll max")
par(mar=c(4,4,1,1))
plot(seq(1,365),dm,ylim = rev(range(c(0,100))), pch = 19,cex = 0.5, col = "darkseagreen",
     xlab = "Day of the year", ylab = substitute(paste(bold("depth [m]"))), main = "Depth of chlorophyll max")











# Surface plot of phytoplankton
par(mfrow=c(1,3))
par(mar=c(4,4,2,1))
image.plot((time)/365, param$z, P, ylim = rev(range(param$z)),col = topo.colors(100)[20:100] ,
           xlab = substitute(paste(bold("Year"))), ylab = substitute(paste(bold("Depth [m]"))),
           main = "Phytoplankton [mmol N/m^3]")
box()


# Surface plot of Nutrients
par(mar=c(4,1,2,4))
image.plot((time)/365, param$z, N, ylim = rev(range(param$z)), col = rev(rainbow(300,start = 0.5,end = 0.6)),
           xlab = substitute(paste(bold("Year"))), ylab = "",
           main = "Nutrient [mmol N/m^3]",yaxt="n")
box()

# Surface plot of Nutrients
par(mar=c(4,1,2,5))

image.plot((time)/365, param$z, D, ylim = rev(range(param$z)), col = terrain.colors(100)[45:100],
           xlab = substitute(paste(bold("Year"))), ylab = "",
           main = "Detritus [mmol N/m^3]",yaxt="n")
box()

