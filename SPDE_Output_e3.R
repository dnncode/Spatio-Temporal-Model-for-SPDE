#******************************************************
#******************************************************
#******************************************************
# Paper: Statistical Modeling for Spatio-Temporal Data from Physical Convection-Diffusion Processes

# ------------------------------------------------------------------------------------------
# Description (See "readme.docx" for more details): 
# Visualization and process the output of SPDE_Analysis_e3.R
# Visualization 1: plot the first actual image, actual wind, and 10 simulated figures
# Visualization 2: plot the first actual image, actual wind, estimated wind, kernel location
# Visualization 3: Kalman Filter and Predictions
# Visualization 4: visualize the filtering results
# Visualization 5: Visualize growth
# Visualization 6: visualize the prediction results
# ------------------------------------------------------------------------------------------
# Main code: 
# SPDE_Output_e3.R
# ------------------------------------------------------------------------------------------
# Subroutines:
# 1). "SPDE_Simulation_Sub.R"
# 2). "SPDE_Analysis_Sub_G.R"
# 3). "Kalman.R"
# ------------------------------------------------------------------------------------------
# Data:
# 1). data_0313_example3.RData (simulated data)
# 2). data_0313_example3_process.RData (processed data by SPDE_Process_e3.R)
# 3). J2_example3_b.RData (output data using 2 kernels for the velocity field)
# 4). J4_example3_c.RData (output data using 4 kernels for the velocity field)
# ------------------------------------------------------------------------------------------
# Authors: XXX
# Last revision: 06/01/2019
#******************************************************
#******************************************************
#******************************************************

rm(list=ls())
wd = "C:/Users/xl027/Desktop/25.SPDE_SS/code/20180814/example3/Github"
setwd(wd)


library(sp)
library(RColorBrewer)
library(fields)
library(MASS)
library(mvtnorm)
library(dlm)
library(expm)
library(spectral)
# load original input data:
load("data_0313_example3.RData")
load("data_0313_example3_process.RData")
#load("processedData/like_opt_onenode_100iter.RData")
source("SPDE_Simulation_Sub.R")
source("SPDE_Analysis_Sub_G.R")
source("Kalman.R")

if (TRUE){
  J=2
  k.center.2 = cbind( c(0.75, 0.25),
                      c(0.5, 0.5))
  # kernel functions; bi-variate normal; sigma = 0.5 (half with)
  kernel.list = list()
  FF = cbind(rep(1,length(grd.original)),coordinates(grd.original)) 
  X.list=list()
  kernel.sigma = 1/J
  for (j in 1:J){
    kernel.list[[j]] = dmvnorm(coordinates(grd.original), mean=c(k.center.2[j,1],k.center.2[j,2]), 
                               sigma=matrix( c(kernel.sigma,0,0,kernel.sigma),ncol=2), log=FALSE)
    X.list[[j]] = diag(kernel.list[[j]] ) %*% FF
  }
  X.2 = do.call(cbind,X.list)
}
if (TRUE){
  J=4
  k.center = cbind( c(0.225, 0.725, 0.225, 0.725),
                    c(0.225, 0.725, 0.725, 0.225))
  # kernel functions; bi-variate normal; sigma = 0.5 (half with)
  kernel.list = list()
  FF = cbind(rep(1,length(grd.original)),coordinates(grd.original)) 
  X.list=list()
  kernel.sigma = 1/J
  for (j in 1:J){
    kernel.list[[j]] = dmvnorm(coordinates(grd.original), mean=c(k.center[j,1],k.center[j,2]), 
                               sigma=matrix( c(kernel.sigma,0,0,kernel.sigma),ncol=2), log=FALSE)
    X.list[[j]] = diag(kernel.list[[j]] ) %*% FF
  }
  X = do.call(cbind,X.list)
}


# load output from SPDE_Analysis
v.max = 0.19
#mle.para = round(mle.para,2)
load("J2_example3_b.RData")
mle.para.J2 =like.opt$par
load("J4_example3_c.RData")
mle.para =like.opt$par
# mle.para needs to be obtained from SPDE_Analysis


colPalette <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'YlOrRd'))(100), .85)
colPalette2 <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Blues'))(100), 0.99)
colPalette1 <- rev(rainbow(34, start=5/6, end=3/6, alpha=0.8))


# ------------------------------------------------
# ------------------------------------------------
# Visualization 1: plot the first actual image, actual wind, and 10 simulated figures
# ------------------------------------------------
# ------------------------------------------------
vortex.center = matrix(c(0.75,0.35),nrow=1)

png(filename = "figures/data.png",
    width = 6, height = 4, unit="in", pointsize = 12,
    bg = "white", res = 600)

layout.matrix = (matrix(c(1,1,1,1,1,2,3,4,5,6,7,8,9,10,11), nrow = 5, ncol = 3))
layout(mat = layout.matrix,
       heights = c(2,2,2,2,2), # Heights of the two rows
       widths = c(5.3,2,2)) # Widths of the two columns

# actual velocity
par(mar=c(2,2,2,2))
ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,1])))
plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
image(ksi.spdf, add=TRUE, col=colPalette1, zlim=c(0,34.1))
plot(grd, col="red",xlim=c(0,1),ylim=c(0,1),add=TRUE)
arrows(coordinates(grd)[,1],coordinates(grd)[,2],
       coordinates(grd)[,1]+mu[,1]/1.2,
       coordinates(grd)[,2]+mu[,2]/1.2,length=0.05,col="blue",lwd=1)
text(x=0.1, y = -0.05, labels = "spatially-varying velocity field with", 
     pos = 4,  
     cex = 1.2)
text(x=0.1, y = -0.10, labels = "a vortex in the southeast region (circled)", 
     pos = 4,  
     cex = 1.2)
vortex.center.spdf = SpatialPoints(vortex.center)
plot(vortex.center.spdf, pch=1, col="green", cex=25, add=TRUE,lwd=2)
title("initial condition (t=0) and velocity field")

for (i in 2:11){
  par(mar=c(1,0.5,1,0.5))
  ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,i])))
  plot(grd, axes=FALSE,col="grey",ylim=c(0,1),xlim=c(0,1),xlab="",ylab="")
  image(ksi.spdf, add=TRUE, col=colPalette1, zlim=c(0,34.1))
  title(paste("t=",i-1,sep=""))
}
dev.off()





# ------------------------------------------------
# ------------------------------------------------
# Visualization 2: plot the first actual image, actual wind, estimated wind, kernel location
# ------------------------------------------------
# ------------------------------------------------
vortex.center = matrix(c(0.75,0.35),nrow=1)
png(filename = "figures/velocity.png",
    width = 8, height = 4, unit="in", pointsize = 12,
    bg = "white", res = 600)
par(mfrow = c(3, 2))
par(las=TRUE)
par(mar=c(4,2.5,2,1))

par(mfrow=c(1,2))
# actual velocity
ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,1])))
plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
image(ksi.spdf, add=TRUE, col=colPalette1)
plot(grd, col="red",xlim=c(0,1),ylim=c(0,1),add=TRUE)
arrows(coordinates(grd)[,1],coordinates(grd)[,2],
       coordinates(grd)[,1]+mu[,1],
       coordinates(grd)[,2]+mu[,2],length=0.05,col="blue",lwd=1)
title("actual velocity field")
mtext("(a)", side=1,line = 2)

# estimated velocity
velocity = convection_fun_1(X.input=X, gamma_x=mle.para[1:(3*J)],gamma_y=mle.para[(3*J+1):(6*J)])
velocity = v.max * tanh(velocity)
ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,1])))
plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
image(ksi.spdf, add=TRUE, col=colPalette1)
plot(grd, col="red",xlim=c(0,1),ylim=c(0,1),add=TRUE)
arrows(coordinates(grd.original)[,1],coordinates(grd.original)[,2],
       coordinates(grd.original)[,1]+velocity[,1]/1,
       coordinates(grd.original)[,2]+velocity[,2]/1,length=0.05,col="blue",lwd=1)
k.center.spdf = SpatialPoints(k.center)
plot(k.center.spdf, pch=15, col="yellow", cex=3, add=TRUE,lwd=1)
vortex.center.spdf = SpatialPoints(vortex.center)
plot(vortex.center.spdf, pch=1, col="green", cex=15, add=TRUE,lwd=2)
for (i in 1:J){
  text(k.center[i,1], k.center[i,2], i, col="darkgreen",cex=1.4)
}
title("estimated velocity field")
mtext("(b)", side=1,line = 2)
#save(velocity, file="data/velocity.RData")
dev.off()


wind.list = list()
png(filename = "figures/velocityonly.png",
    width = 8, height = 4, unit="in", pointsize = 10,
    bg = "white", res = 600)
par(las=TRUE)
par(mar=c(4,4,4,1))
par(mfrow=c(1,2))
# actual velocity
#ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,1])))
plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
#image(ksi.spdf, add=TRUE, col=colPalette1)
#plot(grd, col="grey",xlim=c(0,1),ylim=c(0,1),add=TRUE)
arrows(coordinates(grd)[,1],coordinates(grd)[,2],
       coordinates(grd)[,1]+mu[,1]/1,
       coordinates(grd)[,2]+mu[,2]/1,length=0.05,col="darkgreen",lwd=0.5)
#title("actual velocity field")
# estimated velocity
J = 2
velocity = convection_fun_1(X.input=X.2, gamma_x=mle.para.J2[1:(3*J)],
                            gamma_y=mle.para.J2[(3*J+1):(6*J)])
velocity = v.max * tanh(velocity)
velocity = convection_fun_1(X.input=X.2, gamma_x=mle.para.J2[1:(3*J)],
                            gamma_y=mle.para.J2[(3*J+1):(6*J)])
velocity = v.max * tanh(velocity)

#ksi.spdf = SpatialPixelsDataFrame(grd, data.frame(as.matrix(KSI.matrix[,1])))
#plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
#image(ksi.spdf, add=TRUE, col=colPalette1)
#plot(grd, col="red",xlim=c(0,1),ylim=c(0,1),add=TRUE)
arrows(coordinates(grd.original)[,1],coordinates(grd.original)[,2],
       coordinates(grd.original)[,1]+velocity[,1]/1,
       coordinates(grd.original)[,2]+velocity[,2]/1,length=0.05,code=2,col="red",lwd=1.5)
k.center.spdf = SpatialPoints(k.center.2)
#plot(k.center.spdf, pch=10, col="blue", cex=2, add=TRUE)
for (i in 1:J){
  text(k.center.2[i,1], k.center.2[i,2], i, col="blue")
}
title("actual and estimated velocity field (J=2)")
mtext("(a)", side=1,line = 2)

#mse.x.2 = mean( (velocity-mu)[,1]^2 )
#mse.y.2 =mean( (velocity-mu)[,2]^2 )
wind.list[[1]] = mu
wind.list[[2]] = velocity

J = 4
plot(grd, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
arrows(coordinates(grd)[,1],coordinates(grd)[,2],
       coordinates(grd)[,1]+mu[,1]/1,
       coordinates(grd)[,2]+mu[,2]/1,length=0.05,col="darkgreen",lwd=0.5)
velocity = convection_fun_1(X.input=X, gamma_x=mle.para[1:(3*J)],gamma_y=mle.para[(3*J+1):(6*J)])
velocity = v.max * tanh(velocity)
arrows(coordinates(grd.original)[,1],coordinates(grd.original)[,2],
       coordinates(grd.original)[,1]+velocity[,1]/1,
       coordinates(grd.original)[,2]+velocity[,2]/1,length=0.05,code=2,col="red",lwd=1.5)
#k.center.spdf = SpatialPoints(k.center)
#plot(k.center.spdf, pch=10, col="blue", cex=2, add=TRUE)
for (i in 1:J){
  text(k.center[i,1], k.center[i,2], i, col="blue")
}
title("actual and estimated velocity field (J=4)")
mtext("(b)", side=1,line = 2)

#mse.x.4 = mean( (velocity-mu)[,1]^2 )
#mse.y.4 =mean( (velocity-mu)[,2]^2 )
wind.list[[3]] = velocity

dev.off()

tmp1 = c( wind.list[[2]][,1] - wind.list[[1]][,1], wind.list[[3]][,1] - wind.list[[1]][,1] )
tmp2 = c( wind.list[[2]][,2] - wind.list[[1]][,2], wind.list[[3]][,2] - wind.list[[1]][,2] )
tmp = data.frame( matrix(c(tmp1, tmp2),ncol=1) )
colnames(tmp) = "error"
tmp$number_of_kernals = rep( c(rep(2,400),rep(4,400)), 2)
tmp$directions = c(rep("horizontal",800),rep("vertical",800))
tmp$number_of_kernals = as.factor(tmp$number_of_kernals)
tmp$directions = as.factor(tmp$directions)

library(ggplot2)
png(filename = "figures/box.png",
    width = 6, height = 3, unit="in", pointsize = 10,
    bg = "white", res = 600)
ggplot(tmp, aes(x=number_of_kernals, y=error, color=directions)) +
  geom_boxplot()
dev.off()
# Change the position



# ------------------------------------------------
# ------------------------------------------------
# Visualization 3: Kalman Filter and Predictions
# ------------------------------------------------
# ------------------------------------------------
source("SPDE_Analysis_Sub_G.R") # reload the G function; J4.RData was calculated based on the different G.function
source("Kalman.R") # reload the G function; J4.RData was calculated based on the different G.function


para = mle.para
mu = convection_fun_1(X.input=X, gamma_x=para[1:(3*J)],gamma_y=para[(3*J+1):(6*J)])
r.diffusion = exp(para[6*J+1])

sigma = lapply(c(1:(n.x*n.y)), diffusion_fun, r=r.diffusion,
               x.Coords=coordinates(grd.original)[,1], y.Coords=coordinates(grd.original)[,2], 
               Mu = mu) # a n^2 list, for each list, 2 by 2 matrix
for (i in 1:(n.x*n.y)){
  sigma[[i]][1,1] = mu[i,1] *sigma[[i]][1,1]
  sigma[[i]][2,2] = mu[i,2] *sigma[[i]][2,2]
}
div.sigma = list()
tmp = solve(matrix( c(kernel.sigma,0,0,kernel.sigma),ncol=2))
for (i in 1:(n.x*n.y)){
  x = coordinates(grd.original)[i,1]
  y = coordinates(grd.original)[i,2]
  tmp2 = tmp3 = list()
  for (j in 1:J){
    term1 = -kernel.list[[j]][i] * tmp[1,1] * (x-k.center[j,1])
    term2 = (term1*x + kernel.list[[j]][i] )
    term3 = (term1*y)
    tmp2[[j]] = matrix( c(term1, term2, term3),nrow=1 )
    
    term4 = -kernel.list[[j]][i] * tmp[2,2] * (x-k.center[j,2])
    term5 = (term4*x)
    term6 = (term4*y + kernel.list[[j]][i] )
    tmp3[[j]] = matrix( c(term4, term5, term6),nrow=1 )
  }
  tmp2 = do.call(cbind, tmp2)
  tmp3 = do.call(cbind, tmp3)
  tmp2 = tmp2%*%matrix(Gamma_x,ncol=1)
  tmp3 = tmp3%*%matrix(Gamma_y,ncol=1)
  
  div.sigma[[i]] = matrix(c(tmp2*r.diffusion,tmp3*r.diffusion),nrow=1)
}

J.l = J
zeta.G = 1/(1+exp(-(para[6*J.l+2])))
G.1 = G.function(N.alpha=n.alpha, N.k1=n.k1, N.k2=n.k2,
                 Omega.1.k1=omega.1.k1, Omega.2.k1=omega.2.k1, 
                 Omega.1.k2=omega.1.k2, Omega.2.k2=omega.2.k2,
                 Phi.c.mat.k1=phi.c.mat.k1, Phi.s.mat.k1=phi.s.mat.k1, 
                 Phi.c.mat.k2=phi.c.mat.k2, Phi.s.mat.k2=phi.s.mat.k2,
                 N.x=n.x, N.y=n.y, Delta.x=delta.x, Delta.y=delta.y,
                 MU=mu, SIGMA=sigma, DIV.SIGMA=div.sigma, ZETA=zeta.G,
                 Par=FALSE)  # Par is set to FALSE because optimParallel is used for MLE.



Alpha.l.fun=Alpha
d.t=delta.t
precision=int.precision
V.max=v.max
X.l=X
J.l=J
n.x.l=n.x
n.y.l=n.y
grd.original.l=grd.original
kernel.sigma.l=kernel.sigma
kernel.list.l=kernel.list
k.center.l=k.center
n.alpha.l=n.alpha
n.k1.l= n.k1
n.k2.l=n.k2
omega.1.k1.l=omega.1.k1
omega.2.k1.l=omega.2.k1
omega.1.k2.l=omega.1.k2
omega.2.k2.l=omega.2.k2
phi.c.mat.k1.l= phi.c.mat.k1
phi.s.mat.k1.l= phi.s.mat.k1
phi.c.mat.k2.l= phi.c.mat.k2
phi.s.mat.k2.l= phi.s.mat.k2
delta.x.l= delta.x
delta.y.l= delta.y
Growth=TRUE


d.t = delta.t
precision = 10
# compute GG
# Initial Condition
m0 = Alpha[,1]
C0 = 0.01^2*diag(n.alpha)
GG = as.matrix( expm(G.1*d.t) )
if (Growth){
  RHO=exp(para[6*J.l+5]) * d.t 
  m0 = c( Alpha.l.fun[,1], array(0, dim=c(1,length(Alpha.l.fun[,1]))) )
  C0 = 0.01^2*diag(2 * n.alpha.l)
  GG.tmp1 = as.matrix( expm(G.1*d.t) )
  GG.tmp2 = diag(nrow(G.1))*d.t
  GG.tmp3 = array(0,dim=c(nrow(G.1),ncol(G.1)))
  GG.tmp4 = RHO * diag(nrow(G.1))
  tmp1 = rbind(GG.tmp1,GG.tmp3)
  tmp2 = rbind(GG.tmp2,GG.tmp4)
  GG = cbind(tmp1, tmp2)
}

# calculate the covariance matrix for W.
H = diag( exp(para[6*J+4]),nrow(G.1))
interval.t = seq(d.t/precision/2, d.t, d.t/precision)
if (FALSE){  # parallel foreach
  dummy = foreach(i=1:length(interval.t)) %dopar% {
    tau = interval.t[i]
    aa = as.matrix( expm(G.1*(d.t-tau)) )  %*% H %*% (as.matrix( expm(t(G.1)*(d.t-tau)) ) )
    return(aa)
  }
  S_W = Reduce("+", dummy) * d.t/precision
}
S_W = array(0, dim=dim(G.1))
for (i in 1:length(interval.t)){
  tau = interval.t[i]
  S_W = S_W +  expm(G.1*(d.t-tau))   %*% H %*% (expm(t(G.1)*(d.t-tau)))  * d.t/precision
}
if (FALSE){
  S_W = round(S_W,8)
  for (i in (nrow(S_W)-1)){
    for (j in (i+1):nrow(S_W)){
      S_W[j,i] = S_W[i,j]
    }
  } #S_W[1,2] == S_W[2,1]
}

if (Growth){
  tmp1 = array(0,dim=c(nrow(S_W),ncol(S_W)))
  tmp2 = exp(para[6*J.l+6]) * diag(nrow(S_W))
  S_W = cbind( rbind(S_W,tmp1), rbind(tmp1, tmp2) )
}


FF = diag(nrow(Alpha))
if (Growth){
  FF = cbind(FF, array(0,dim=c(nrow(FF),ncol(FF))))
}

if (TRUE){ # run the KL filter using my own function
  Kalman.output = Kalman(data.alpha=Alpha[,1:10], m0, C0, FF, GG=GG, s_W=S_W, s_V=exp(para[6*J.l+3]))
  # evaluate MLE. 
  f = Kalman.output[[3]]
  Q = Kalman.output[[4]]
}
if (FALSE){ # run the KL filter using the function from dlm
  DSTM = dlm(m0=m0, C0=C0, FF=FF, GG=GG, W=S_W, V=exp(para[6*J+3])^2*diag(nrow(Alpha)))
  Filt = dlmFilter(Alpha[,1:10], DSTM)
  f = Filt$f
  Q = list()
  for (i in 1:ncol(f)){
    #R = Filt$U.R[[i]] %*% diag(Filt$D.R[i,]^2) %*% t(Filt$U.R[[i]])
    R = dlmSvd2var(Filt$U.R[[i]], Filt$D.R[i,])
    Q[[i]] = FF%*%R%*%t(FF) + exp(para[6*J+3])^2*diag(n.alpha)
  }
}



#******************************************************
#******************************************************
# Visualization 4: visualize the filtering results
#******************************************************
#******************************************************
a = Kalman.output[[5]]
#a = t( Filt$a )
n.a = nrow(a)
a = a[1:(n.a/2),]
for (i.t in 1:(ncol(a)-1)){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  #alpha.c.1 = matrix(Filt$m[i.t+1,1:4] , ncol=1)
  #alpha.c.2 = matrix(Filt$m[i.t+1,5: ((ncol(Filt$m)-4)/2+4)], ncol=1)
  #alpha.s.2 = matrix( Filt$m[i.t+1, ((ncol(Filt$m)-4)/2+5):ncol(Filt$m)], ncol=1)
  
  alpha.c.1 = matrix(a[1:4,i.t] , ncol=1)
  alpha.c.2 = matrix(a[5: ((nrow(a)-4)/2+4), i.t], ncol=1)
  alpha.s.2 = matrix(a[((nrow(a)-4)/2+5):nrow(a), i.t], ncol=1)
  
  # alpha.c.1 = matrix(Alpha[1:4,i.t] , ncol=1)
  # alpha.c.2 = matrix(Alpha[5: ((ncol(Filt$m)-4)/2+4), i.t], ncol=1)
  # alpha.s.2 = matrix( Alpha[((ncol(Filt$m)-4)/2+5):ncol(Filt$m), i.t], ncol=1)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  #alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    #alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  # alpha.c.2 = array(0/0,dim=c(n.k2,1))
  #alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    #alpha.c.2[i,1]= Re( my.s[case1,case2] )
    #alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  #par(mfrow = c(1, 2))
  #rasterImage2(x = grd.x,
  #             y = grd.y,
  #            z = z,palette=colPalette1)
  #rasterImage2(x = grd.x,
  #            y = grd.y,
  #            z = (t(matrix((z.recover.3),nrow=n.y))),palette=colPalette1)
  
  Min = 0
  Max = max(max(z), max(z.recover.3))
  
  # output figures
  png(filename = paste("figures/filter_", i.t, ".png", sep=""),
      width = 12, height =6, unit="in", pointsize = 12,
      bg = "white", res = 600)
  par(mfrow = c(1, 2))
  par(las=TRUE)
  par(mar=c(4,4,2,1))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z,xlab="",ylab="",main=paste("actual data (","time:",i.t,") ",sep=""),
               zlim=c( Min,Max),palette=colPalette1)
  rasterImage2(x = grd.x,
               y = grd.y,
               z = (t(matrix((z.recover.3),nrow=n.y))),
               xlab="",ylab="",main=paste("filtered data (","time:",i.t,") ",sep=""),zlim=c( Min,Max),palette=colPalette1)
  dev.off()
  
}

# only visualize 2, 5, 10
png(filename = paste("figures/filter_big", ".png", sep=""),
    width = 12, height =4.1, unit="in", pointsize = 12,
    bg = "white", res = 600)
par(mfcol = c(2, 5))
par(las=TRUE)
par(mar=c(2,2,2,1))
for (i.t in seq(2,10,2)){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  #alpha.c.1 = matrix(Filt$m[i.t+1,1:4] , ncol=1)
  #alpha.c.2 = matrix(Filt$m[i.t+1,5: ((ncol(Filt$m)-4)/2+4)], ncol=1)
  #alpha.s.2 = matrix( Filt$m[i.t+1, ((ncol(Filt$m)-4)/2+5):ncol(Filt$m)], ncol=1)
  
  alpha.c.1 = matrix(a[1:4,i.t] , ncol=1)
  alpha.c.2 = matrix(a[5: ((nrow(a)-4)/2+4), i.t], ncol=1)
  alpha.s.2 = matrix(a[((nrow(a)-4)/2+5):nrow(a), i.t], ncol=1)
  
  # alpha.c.1 = matrix(Alpha[1:4,i.t] , ncol=1)
  # alpha.c.2 = matrix(Alpha[5: ((ncol(Filt$m)-4)/2+4), i.t], ncol=1)
  # alpha.s.2 = matrix( Alpha[((ncol(Filt$m)-4)/2+5):ncol(Filt$m), i.t], ncol=1)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  #alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    #alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  # alpha.c.2 = array(0/0,dim=c(n.k2,1))
  #alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    #alpha.c.2[i,1]= Re( my.s[case1,case2] )
    #alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  
  Min = 0
  Max = max(max(z), max(z.recover.3))
  
  # output figures
  
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z,xlab="",ylab="",main=paste("actual image (","time:",i.t,") ",sep=""),
               zlim=c( Min,Max),palette=colPalette1, cex.main=1.5)
  rasterImage2(x = grd.x,
               y = grd.y,
               z = (t(matrix((z.recover.3),nrow=n.y))),cex.main=1.5,
               xlab="",ylab="",main=paste("filtered image (","time:",i.t,") ",sep=""),zlim=c( Min,Max),palette=colPalette1)
  
  
}
dev.off()


#******************************************************
#******************************************************
# Visualization 5: Visualize growth
#******************************************************
#******************************************************
a = Kalman.output[[5]]
#a = t( Filt$a )
n.a = nrow(a)
a = a[(n.a/2):n.a,]
#******************************************************

for (i.t in 1:(ncol(a)-1)){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  #alpha.c.1 = matrix(Filt$m[i.t+1,1:4] , ncol=1)
  #alpha.c.2 = matrix(Filt$m[i.t+1,5: ((ncol(Filt$m)-4)/2+4)], ncol=1)
  #alpha.s.2 = matrix( Filt$m[i.t+1, ((ncol(Filt$m)-4)/2+5):ncol(Filt$m)], ncol=1)
  
  alpha.c.1 = matrix(a[1:4,i.t] , ncol=1)
  alpha.c.2 = matrix(a[5: ((nrow(a)-4)/2+4), i.t], ncol=1)
  alpha.s.2 = matrix(a[((nrow(a)-4)/2+5):nrow(a), i.t], ncol=1)
  
  # alpha.c.1 = matrix(Alpha[1:4,i.t] , ncol=1)
  # alpha.c.2 = matrix(Alpha[5: ((ncol(Filt$m)-4)/2+4), i.t], ncol=1)
  # alpha.s.2 = matrix( Alpha[((ncol(Filt$m)-4)/2+5):ncol(Filt$m), i.t], ncol=1)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  #alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    #alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  # alpha.c.2 = array(0/0,dim=c(n.k2,1))
  #alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    #alpha.c.2[i,1]= Re( my.s[case1,case2] )
    #alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  #par(mfrow = c(1, 2))
  #rasterImage2(x = grd.x,
  #             y = grd.y,
  #            z = z,palette=colPalette1)
  #rasterImage2(x = grd.x,
  #            y = grd.y,
  #            z = (t(matrix((z.recover.3),nrow=n.y))),palette=colPalette1)
  
  Min = 0
  Max = max(max(z), max(z.recover.3))
  
  # output figures
  png(filename = paste("figures/source_", i.t, ".png", sep=""),
      width = 12, height =6, unit="in", pointsize = 12,
      bg = "white", res = 600)
  par(mfrow = c(1, 2))
  par(las=TRUE)
  par(mar=c(4,4,2,1))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z,xlab="",ylab="",main=paste("actual data (","time:",i.t,") ",sep=""),
               zlim=c( Min,Max),palette=colPalette1)
  rasterImage2(x = grd.x,
               y = grd.y,
               z = (t(matrix((z.recover.3),nrow=n.y))),
               xlab="",ylab="",main=paste("filtered data (","time:",i.t,") ",sep=""),palette=colPalette1)
  dev.off()
  
}


#******************************************************
#******************************************************
# Visualization 6: visualize the prediction results
#******************************************************
#******************************************************

a = Kalman.output[[5]]
aa = matrix( a[,ncol(a)], ncol=1)
jj = 1
n.aa = nrow(aa)
for (i.t in (ncol(a):(ncol(a)+8))){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  #alpha.c.1 = matrix(Filt$m[i.t+1,1:4] , ncol=1)
  #alpha.c.2 = matrix(Filt$m[i.t+1,5: ((ncol(Filt$m)-4)/2+4)], ncol=1)
  #alpha.s.2 = matrix( Filt$m[i.t+1, ((ncol(Filt$m)-4)/2+5):ncol(Filt$m)], ncol=1)
  aa = GG %*% aa
  
  aaa = aa[1:(n.aa/2),1]
  alpha.c.1 = matrix(aaa[1:4] , ncol=1)
  alpha.c.2 = matrix(aaa[5: ((length(aaa)-4)/2+4)], ncol=1)
  alpha.s.2 = matrix(aaa[((length(aaa)-4)/2+5):length(aaa)], ncol=1)
  
  # alpha.c.1 = matrix(Alpha[1:4,i.t] , ncol=1)
  # alpha.c.2 = matrix(Alpha[5: ((ncol(Filt$m)-4)/2+4), i.t], ncol=1)
  # alpha.s.2 = matrix( Alpha[((ncol(Filt$m)-4)/2+5):ncol(Filt$m), i.t], ncol=1)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  #alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    #alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  # alpha.c.2 = array(0/0,dim=c(n.k2,1))
  #alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    #alpha.c.2[i,1]= Re( my.s[case1,case2] )
    #alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  z.pred = (t(matrix((z.recover.3),nrow=n.y)))
  Min = 0
  Max = max(max(z), max(z.recover.3))
  
  
  # output figures
  png(filename = paste("figures/pred_", i.t, ".png", sep=""),
      width = 12, height =6, unit="in", pointsize = 12,
      bg = "white", res = 600)
  par(mfrow = c(1, 2))
  par(las=TRUE)
  par(mar=c(4,4,2,1))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z,xlab="",ylab="",main=paste("actual data (","time: +",jj,") ",sep=""),
               zlim=c( Min,Max))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = (t(matrix((z.recover.3),nrow=n.y))),
               xlab="",ylab="",main=paste("predicted data (","time: +",jj,") ",sep=""),zlim=c( Min,Max))
  dev.off()
  jj = jj+1
}


# Visualization 2: only visualize 2, 5, 10
png(filename = paste("figures/pred_big", ".png", sep=""),
    width = 8, height =5, unit="in", pointsize = 12,
    bg = "white", res = 600)
par(mfcol = c(2, 3))
par(las=TRUE)
par(mar=c(4,4,2,1))

a = Kalman.output[[5]]
aa = matrix( a[,ncol(a)], ncol=1)
jj = 1
n.aa = nrow(aa)
for (i.t in (ncol(a):(ncol(a)+10))){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  #alpha.c.1 = matrix(Filt$m[i.t+1,1:4] , ncol=1)
  #alpha.c.2 = matrix(Filt$m[i.t+1,5: ((ncol(Filt$m)-4)/2+4)], ncol=1)
  #alpha.s.2 = matrix( Filt$m[i.t+1, ((ncol(Filt$m)-4)/2+5):ncol(Filt$m)], ncol=1)
  aa = GG %*% aa
  
  aaa = aa[1:(n.aa/2),1]
  alpha.c.1 = matrix(aaa[1:4] , ncol=1)
  alpha.c.2 = matrix(aaa[5: ((length(aaa)-4)/2+4)], ncol=1)
  alpha.s.2 = matrix(aaa[((length(aaa)-4)/2+5):length(aaa)], ncol=1)
  
  # alpha.c.1 = matrix(Alpha[1:4,i.t] , ncol=1)
  # alpha.c.2 = matrix(Alpha[5: ((ncol(Filt$m)-4)/2+4), i.t], ncol=1)
  # alpha.s.2 = matrix( Alpha[((ncol(Filt$m)-4)/2+5):ncol(Filt$m), i.t], ncol=1)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  #alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    #alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  # alpha.c.2 = array(0/0,dim=c(n.k2,1))
  #alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    #alpha.c.2[i,1]= Re( my.s[case1,case2] )
    #alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  
  Min = 0
  Max = max(max(z), max(z.recover.3))
  
  if (i.t %in% c(11,14,20)){
    rasterImage2(x = grd.x,
                 y = grd.y,
                 z = z,xlab="",ylab="",main=paste("actual image (","time: +",i.t-10,") ",sep=""),
                 zlim=c( Min,Max))
    rasterImage2(x = grd.x,
                 y = grd.y,
                 z = (t(matrix((z.recover.3),nrow=n.y))),
                 xlab="",ylab="",main=paste("predicted image (","time: +",i.t-10,") ",sep=""),zlim=c( Min,Max))
  }
  
  jj = jj+1
}
dev.off()























