#******************************************************
#******************************************************
#******************************************************
# Paper: Statistical Modeling for Spatio-Temporal Data from Physical Convection-Diffusion Processes

# ------------------------------------------------------------------------------------------
# Description (See "readme.docx" for more details): 
# Subroutines for SPDE_Analylsis_e3.R
# ------------------------------------------------------------------------------------------
# Authors: XXX
# Last revision: 06/01/2019
#******************************************************
#******************************************************
#******************************************************

library(sp)
#library(rgdal)
#library(rgeos)
#library(gstat)
library(RColorBrewer)
#library(spacetime)
#library(selextR)
library(fields)
library(MASS)
library(mvtnorm)
library(matrixcalc)
library(spectral)
library(dlm)
library(expm)
library(optimParallel)
#library(Matrix)



# clusterEvalQ(cl, source("SPDE_Analysis_Cluster_Loading.R"))
# lines 58 to 318 are loaded by the above line 
# lines 58 to 371 are loaded again below in main program

#******************************************************
#******************************************************
# Load data
#******************************************************
#******************************************************
load("data/data_1030_example2_process.RData")
#load("data/data_1030_example2.RData")
load("data/data_1030_example2_ksi.RData")
load("data/data_1030_example2_grd_original.RData")
load("data/data_1030_example2_G.RData")
load("data/data_1030_example2_Noise.RData")
load("data/data_1030_example2_para.RData")
delta.t=para.setting[3]
n.alpha = nrow(Alpha)

my.x = seq(-n.x/2+1, n.x/2, 1)
my.y = seq(-n.y/2+1, n.y/2, 1)
grd.x = unique(coordinates(grd.original)[,1])
grd.y = ( unique(coordinates(grd.original)[,2]) ) 
grd.xx = c(0:(length(grd.x)-1))
grd.yy = c(0:(length(grd.y)-1))

omega.1.k1 = 2*pi*k1.set[,1]
omega.2.k1 =2*pi*k1.set[,2]
omega.1.k2 = 2*pi*k2.set[,1]
omega.2.k2 =2*pi*k2.set[,2]

delta.x = 1/n.x # spatial resolution;
delta.y = 1/n.y


tmp = lapply(c(1:n.k1), phi.c, Kx=omega.1.k1, Ky=omega.2.k1, Grd=grd.original)
phi.c.mat.k1 = do.call(cbind, tmp)  # the cosine values at each location for each basis
tmp = lapply(c(1:n.k1), phi.s, Kx=omega.1.k1, Ky=omega.2.k1, Grd=grd.original)
phi.s.mat.k1 = do.call(cbind, tmp)  # the sine values at each location for each basis
tmp = lapply(c(1:n.k2), phi.c, Kx=omega.1.k2, Ky=omega.2.k2, Grd=grd.original)
phi.c.mat.k2 = do.call(cbind, tmp)  # the cosine values at each location for each basis
tmp = lapply(c(1:n.k2), phi.s, Kx=omega.1.k2, Ky=omega.2.k2, Grd=grd.original)
phi.s.mat.k2 = do.call(cbind, tmp)  # the sine values at each location for each basis
phi.c.mat = cbind(phi.c.mat.k1,phi.c.mat.k2)
phi.s.mat = cbind(phi.s.mat.k1,phi.s.mat.k2)

v.max = max(grd.x)/5  # the maximum velocity allowed

#******************************************************
#******************************************************
# MLE
#******************************************************
#******************************************************

# parameter
# this is the assumed convection_fun and diffusion_fun; kernal approximation

if (FALSE){
  J=4
  k.center = cbind( c(0.225, 0.725, 0.225, 0.725),
                    c(0.225, 0.725, 0.725, 0.225))
}
if (TRUE){
  J=3
  k.center = cbind( c(0.742, 0.475, 0.208),
                    c(0.292, 0.78, 0.292))
}



# kernel functions; bi-variate normal; sigma = 0.5 (half with)
kernel.list = list()
FF = cbind(rep(1,length(grd.original)),coordinates(grd.original)) 
X.list=list()
kernel.sigma = 0.25
for (j in 1:J){
  kernel.list[[j]] = dmvnorm(coordinates(grd.original), mean=c(k.center[j,1],k.center[j,2]), 
                             sigma=matrix( c(kernel.sigma,0,0,kernel.sigma),ncol=2), log=FALSE)
  X.list[[j]] = diag(kernel.list[[j]] ) %*% FF
}
X = do.call(cbind,X.list)

# parameter
# this is the assumed convection_fun and diffusion_fun
if (FALSE){ # set to FALSE, because convection_fun_1 has been loaded from SPDE_Simulation_Sub.R
  convection_fun_1 = function(X.input, gamma_x, gamma_y){
    mu_x = X %*% gamma_x
    mu_y = X %*% gamma_y
    output = cbind(mu_x, mu_y)
    output
  }
}
if (FALSE){
  velocity = convection_fun_1(X.input=X, gamma_x=rnorm(J*3,0,0.1),gamma_y=rnorm(J*3,0,0.1))
  velocity = v.max * tanh(velocity )  # this is used in the analysis part; need to bound the velocity
  plot(grd.original)
  arrows(coordinates(grd.original)[,1],coordinates(grd.original)[,2],
         coordinates(grd.original)[,1]+velocity[,1]/2,
         coordinates(grd.original)[,2]+velocity[,2]/2,length=0.1,col="darkgreen",lwd=1)
}
if (FALSE){ # set to FALSE, because convection_fun_1 has been loaded from SPDE_Simulation_Sub.R
  diffusion_fun = function(i, x.Coords, y.Coords, Mu, r){
    if (TRUE){ # diffusion formulation I
      output =  matrix(c(r,0,0,r),ncol=2)
    }
    output
  }
}

# A list of parameters include:
# 1). Gamma_x, Gamma_y --> velocity field
# 2). r.diffusion --> assuming diffusion formulation I
# 3)a. zeta --> assuming constant decay
# 3)b. Gamma_z --> assuming spatially varying decay
# 4). tau_tv --> tau tilde(v)
# 5). rho and tau_beta --> beta
# 6). spectral density of tau epsilon

#set.seed(1)
#Gamma_x = rnorm(J*3,0,0.5)/1.5
#Gamma_y = rnorm(J*3,0,0.5)/1.5
#para=c(Gamma_x,Gamma_y,0.01,0.98,0.01,0.01)
#para=c(Gamma_x,Gamma_y,log(0.1),0.9,log(0.1),log(0.01))


int.precision = 10 # number of intervals for calculate the covariance matrix of W


l.fun = function(para, Alpha.l.fun, d.t, precision, V.max,
                 X.l, J.l, n.x.l, n.y.l, grd.original.l,
                 kernel.sigma.l,kernel.list.l,
                 k.center.l, n.alpha.l,
                 n.k1.l, n.k2.l,
                 omega.1.k1.l,
                 omega.2.k1.l,
                 omega.1.k2.l ,
                 omega.2.k2.l,
                 phi.c.mat.k1.l,
                 phi.s.mat.k1.l,
                 phi.c.mat.k2.l,
                 phi.s.mat.k2.l,
                 delta.x.l,
                 delta.y.l,
                 Growth){  # Growth==TRUE, consider growth and decay
  
  # para:
  # para[1:3J]:Gamma_x
  # para[(3J+1):6J]:Gamma_y
  # para[6J+1]:r in the diffusion matrix
  # para[6J+2]:zeta
  # para[6J+3]:tau_tv --> tau tilde(v)
  # para[6J+4]:spectral density of tau epsilon
  # para[6J+5]:rho used when Q is TRUE
  
  # Alpha.l.fun = Alpha
  
  mu = convection_fun_1(X.input=X.l, gamma_x=para[1:(3*J.l)],gamma_y=para[(3*J.l+1):(6*J.l)])
  mu = V.max * tanh(mu) # this is used in the analysis part; need to bound the velocity
  r.diffusion = exp(para[6*J.l+1])
  
  sigma = lapply(c(1:(n.x.l*n.y.l)), diffusion_fun, r=r.diffusion,
                 x.Coords=coordinates(grd.original.l)[,1], y.Coords=coordinates(grd.original.l)[,2], 
                 Mu = mu) # a n^2 list, for each list, 2 by 2 matrix
  for (i in 1:(n.x.l*n.y.l)){
    sigma[[i]][1,1] = mu[i,1] *sigma[[i]][1,1]
    sigma[[i]][2,2] = mu[i,2] *sigma[[i]][2,2]
  }
  div.sigma = list()
  tmp = solve(matrix( c(kernel.sigma.l,0,0,kernel.sigma.l),ncol=2))
  for (i in 1:(n.x.l*n.y.l)){
    x = coordinates(grd.original.l)[i,1]
    y = coordinates(grd.original.l)[i,2]
    tmp2 = tmp3 = list()
    for (j in 1:J.l){
      term1 = -kernel.list.l[[j]][i] * tmp[1,1] * (x-k.center.l[j,1])
      term2 = (term1*x + kernel.list.l[[j]][i] )
      term3 = (term1*y)
      tmp2[[j]] = matrix( c(term1, term2, term3),nrow=1 )
      
      term4 = -kernel.list.l[[j]][i] * tmp[2,2] * (x-k.center.l[j,2])
      term5 = (term4*x)
      term6 = (term4*y + kernel.list.l[[j]][i] )
      tmp3[[j]] = matrix( c(term4, term5, term6),nrow=1 )
    }
    tmp2 = do.call(cbind, tmp2)
    tmp3 = do.call(cbind, tmp3)
    tmp2 = tmp2%*%matrix(para[1:(3*J.l)],ncol=1)  
    tmp3 = tmp3%*%matrix(para[(3*J.l+1):(6*J.l)],ncol=1)
    
    div.sigma[[i]] = matrix(c(tmp2*r.diffusion,tmp3*r.diffusion),nrow=1)
  }
  
  G = G.function(N.alpha=n.alpha.l, N.k1=n.k1.l, N.k2=n.k2.l,
                 Omega.1.k1=omega.1.k1.l, Omega.2.k1=omega.2.k1.l, 
                 Omega.1.k2=omega.1.k2.l, Omega.2.k2=omega.2.k2.l,
                 Phi.c.mat.k1=phi.c.mat.k1.l, Phi.s.mat.k1=phi.s.mat.k1.l, 
                 Phi.c.mat.k2=phi.c.mat.k2.l, Phi.s.mat.k2=phi.s.mat.k2.l,
                 N.x=n.x.l, N.y=n.y.l, Delta.x=delta.x.l, Delta.y=delta.y.l,
                 MU=mu, SIGMA=sigma, DIV.SIGMA=div.sigma, ZETA=para[6*J.l+2], RHO=para[6*J+5], 
                 Q=Growth, Par=FALSE)
  
  # compute GG
  # Initial Condition
  m0 = Alpha.l.fun[,1]
  C0 = 0.01^2*diag(n.alpha.l)
  GG = as.matrix( expm(G*d.t) )
  if (Growth){
    m0 = c( Alpha.l.fun[,1], array(0, dim=c(1,length(Alpha.l.fun[,1]))) )
    C0 = 0.01^2*diag(2 * n.alpha.l)
    GG = as.matrix( expm(G*d.t) )
  }
  
  # calculate the covariance matrix for W.
  H = diag( exp(para[6*J.l+4]),nrow(GG))
  interval.t = seq(d.t/precision/2, d.t, d.t/precision)
  if (FALSE){  # parallel foreach
    dummy = foreach(i=1:length(interval.t)) %dopar% {
      tau = interval.t[i]
      aa = as.matrix( expm(G*(d.t-tau)) )  %*% H %*% (as.matrix( expm(t(G)*(d.t-tau)) ) )
      return(aa)
    }
    S_W = Reduce("+", dummy) * d.t/precision
  }
  S_W = array(0, dim=dim(G))
  for (i in 1:length(interval.t)){
    tau = interval.t[i]
    S_W = S_W +  expm(G*(d.t-tau))   %*% H %*% (expm(t(G)*(d.t-tau)))  * d.t/precision
  }
  #S_W = round(S_W,3)
  for (i in (nrow(S_W)-1)){
    for (j in (i+1):nrow(S_W)){
      S_W[j,i] = S_W[i,j]
    }
  }
  if (Growth){
    tmp1 = array(0,dim=c(nrow(S_W),ncol(S_W)))
    tmp2 = exp(para[6*J+6]) * diag(nrow(S_W))
    S_W = cbind( rbind(S_W,tmp1), rbind(tmp1, tmp2) )
  }
  
  
  #DSTM = dlm(m0=m0, C0=C0, FF=diag(n.alpha), GG=GG, W=S_W, V=exp(para[6*J+3])^2*diag(n.alpha))
  #Filt = dlmFilter(Alpha.l.fun[,1:10], DSTM)
  FF = diag(nrow(Alpha.l.fun[,1:10]))
  if (Growth){
    FF = cbind(FF, array(0,dim=c(nrow(FF),ncol(FF))))
  }
  Kalman.output = Kalman(data.alpha=Alpha.l.fun[,1:10], m0, C0, FF, GG=G.list[[2]], s_W=S_W, s_V=exp(para[6*J.l+3]))
  
  # evaluate MLE. 
  f = Kalman.output[[3]]
  Q = Kalman.output[[4]]
  if (FALSE){
    l = foreach(i=1:ncol(f)) %dopar% {
      q = Q[[i]]
      #output = -log(det(q))/2 - t(Alpha[,i]-f[,i])  %*% solve(q) %*% (Alpha[,i]-f[,i])/2
      output = dmvnorm(Alpha.l.fun[,i], 
                       mean = f[,i], sigma = q, log = TRUE)
      return(output)
    }
    L = do.call(sum, l)
  }
  L = 0
  for (i in 1:ncol(f)){
    q = Q[[i]]
    #output = -log(det(q))/2 - t(Alpha[,i]-f[,i])  %*% solve(q) %*% (Alpha[,i]-f[,i])/2
    output = dmvnorm(Alpha.l.fun[,i], 
                     mean = f[,i], sigma = q, log = TRUE)
    L = L + output
  }
  L = -L #nagative loglikelihood
  return(L)
}



