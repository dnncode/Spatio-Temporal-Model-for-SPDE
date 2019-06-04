#******************************************************
#******************************************************
#******************************************************
# Paper: Statistical Modeling for Spatio-Temporal Data from Physical Convection-Diffusion Processes

# ------------------------------------------------------------------------------------------
# Description (See "readme.docx" for more details): 
# This code is used to process the data generated from SPDE_simulation_e3.R or other sources
# The processed data will be used by SPDE_Analysis_e3.R
# ------------------------------------------------------------------------------------------
# Main code: 
# SPDE_Process_e3.R
# ------------------------------------------------------------------------------------------
# Subroutines:
# 1). N/A
# ------------------------------------------------------------------------------------------
# Data:
# 1). data_0313_example3_ksi.RData (generated from SPDE_Simulation_e3.R)
# 2). data_0313_example3_grd_original.RData (generated from SPDE_Simulation_e3.R)
# ------------------------------------------------------------------------------------------
# Authors: XXX
# Last revision: 06/01/2019
#******************************************************
#******************************************************
#******************************************************


rm(list=ls())
wd = "C:/Users/xl027/Desktop/25.SPDE_SS/code/20180814"
setwd(wd)

#******************************************************
#******************************************************
# Packages needed
#******************************************************
#******************************************************
library(sp)
library(rgdal)
library(rgeos)
library(gstat)
library(RColorBrewer)
library(spacetime)
library(selextR)
library(timeSeries)
library(fields)
library(MASS)
library(snow)
library(mvtnorm)
library(matrixcalc)
library(spectral)

colPalette <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'YlOrRd'))(100), .85)
colPalette2 <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Blues'))(100), 0.99)
colPalette1 <- rev(rainbow(34, start=5/6, end=3/6, alpha=0.8))


#******************************************************
#******************************************************
# Load data to be processed (need the grid setting and the data)
#******************************************************
#******************************************************
if (TRUE){  #for example 3
  load("data/data_0313_example3_ksi.RData")
  load("data/data_0313_example3_grd_original.RData")
}
#******************************************************
#******************************************************
# Grid information
#******************************************************
#******************************************************
grd.x = unique(coordinates(grd.original)[,1])
grd.y = ( unique(coordinates(grd.original)[,2]) ) 
grd.xx = c(0:(length(grd.x)-1))
grd.yy = c(0:(length(grd.y)-1))
KSI.matrix = data
n.t = ncol(KSI.matrix)
n.x = length(grd.x)
n.y = length(grd.y)
n.alpha = nrow(KSI.matrix)

#******************************************************
#******************************************************
# FFT for each time step, and get the alpha in my setting (the default setting from mvfft is not the same as mine)
#******************************************************
#******************************************************
Alpha = array(0/0,dim=c(n.alpha, n.t))
for (i.t in 1:n.t){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  # FFT
  tmp = t(mvfft(t(z)))
  FT = mvfft(tmp)/n.x/n.y
  mvfft.x = seq(0,n.x-1,1)
  mvfft.y = seq(0,n.y-1,1)
  mvfft.z =  FT
  
  # convert FFT results to my settings
  my.x = seq(-n.x/2+1, n.x/2, 1)
  my.y = seq(-n.y/2+1, n.y/2, 1)
  my.s = array(0/0, dim=c(n.x,n.y))
  for (i in 1:n.x){
    for (j in 1:n.y){
      target.x = my.x[i]
      target.y = my.y[j]
      
      case1 = which(mvfft.x == target.x)
      case2 = which(mvfft.y == target.y) 
      
      if (length(case1)+length(case2)==2){
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case1,case2]), 
                             imaginary = Im(FT[case1,case2]))
      }
      
      if (length(case1)-length(case2)==1){
        case3 = which(mvfft.x == target.x)
        case4 = which(mvfft.y == target.y+n.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)-length(case1)==1){
        case3 = which(mvfft.x == target.x+n.x)
        case4 = which(mvfft.y == target.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)+length(case1)==0){
        case3 = which(mvfft.x == -target.x)
        case4 = which(mvfft.y == -target.y) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = -Im(FT[case3,case4]))
      }
    }
  }
  
  k1.set = cbind(c(0,0,n.x/2,n.x/2),c(0,n.y/2,0,n.y/2))
  k2.set.part1 = cbind( rep(seq(0,n.x/2,1),each=(n.y/2+1)), rep(seq(0,n.y/2,1),(n.x/2+1)))
  tmp1 = paste(k1.set[,1],k1.set[,2],"-")
  tmp2 = paste(k2.set.part1[,1],k2.set.part1[,2],"-")
  case = which( is.na(match(tmp2, tmp1)) )
  k2.set.part1 = k2.set.part1[case,]
  k2.set.part2 = cbind( rep(seq(1,n.x/2-1,1),each=(n.y/2-1)), rep(seq(-1,-(n.y/2-1),-1),(n.x/2-1)))
  k2.set = rbind(k2.set.part1,k2.set.part2)
  n.k1 = nrow(k1.set)
  n.k2 = nrow(k2.set)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  alpha.c.2 = array(0/0,dim=c(n.k2,1))
  alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    alpha.c.2[i,1]= Re( my.s[case1,case2] )
    alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  par(mfrow = c(1, 2))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z)
  rasterImage2(x = grd.x,
               y = grd.y,
               z = (t(matrix((z.recover.3),nrow=n.y))))
  
  alpha = c(alpha.c.1, alpha.c.2, alpha.s.2)
  Alpha[,i.t] = alpha
  
  # output figures
  png(filename = paste("../../figures/numerical3/figure_", i.t, ".png", sep=""),
      width = 8, height =3, unit="in", pointsize = 12,
      bg = "white", res = 600)
  par(mfrow = c(1, 3))
  par(las=TRUE)
  par(mar=c(4,4,2,1))
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z)
  rasterImage2(x = my.x,
               y = my.y,
               z = round((matrix(Re(my.s),nrow=n.x,byrow = FALSE)),2),
               xlab=expression('waver number, k'[1]),
               ylab=expression('waver number, k'[2]),
               main="real part")
  rasterImage2(x = my.x,
               y = my.y,
               z = round((matrix(Im(my.s),nrow=n.x,byrow = FALSE))),
               xlab=expression('waver number, k'[1]),
               ylab=expression('waver number, k'[2]),
               main="imaginary part")
  dev.off()
  
}


# visualization 2: put the first 6 radar images on the same figure
Alpha = array(0/0,dim=c(n.alpha, n.t))
png(filename = paste("../../figures/numerical3/figure_big", ".png", sep=""),
       width = 8, height = 5.5, unit="in", pointsize = 12,
        bg = "white", res = 600)
par(mfcol = c(2, 3))
par(las=TRUE)
par(mar=c(4,4,2,1))

for (i.t in 1:6){
  # data
  z = t( matrix(KSI.matrix[,i.t],nrow=n.y) )
  
  # FFT
  tmp = t(mvfft(t(z)))
  FT = mvfft(tmp)/n.x/n.y
  mvfft.x = seq(0,n.x-1,1)
  mvfft.y = seq(0,n.y-1,1)
  mvfft.z =  FT
  
  # convert FFT results to my settings
  my.x = seq(-n.x/2+1, n.x/2, 1)
  my.y = seq(-n.y/2+1, n.y/2, 1)
  my.s = array(0/0, dim=c(n.x,n.y))
  for (i in 1:n.x){
    for (j in 1:n.y){
      target.x = my.x[i]
      target.y = my.y[j]
      
      case1 = which(mvfft.x == target.x)
      case2 = which(mvfft.y == target.y) 
      
      if (length(case1)+length(case2)==2){
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case1,case2]), 
                             imaginary = Im(FT[case1,case2]))
      }
      
      if (length(case1)-length(case2)==1){
        case3 = which(mvfft.x == target.x)
        case4 = which(mvfft.y == target.y+n.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)-length(case1)==1){
        case3 = which(mvfft.x == target.x+n.x)
        case4 = which(mvfft.y == target.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)+length(case1)==0){
        case3 = which(mvfft.x == -target.x)
        case4 = which(mvfft.y == -target.y) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = -Im(FT[case3,case4]))
      }
    }
  }
  
  k1.set = cbind(c(0,0,n.x/2,n.x/2),c(0,n.y/2,0,n.y/2))
  k2.set.part1 = cbind( rep(seq(0,n.x/2,1),each=(n.y/2+1)), rep(seq(0,n.y/2,1),(n.x/2+1)))
  tmp1 = paste(k1.set[,1],k1.set[,2],"-")
  tmp2 = paste(k2.set.part1[,1],k2.set.part1[,2],"-")
  case = which( is.na(match(tmp2, tmp1)) )
  k2.set.part1 = k2.set.part1[case,]
  k2.set.part2 = cbind( rep(seq(1,n.x/2-1,1),each=(n.y/2-1)), rep(seq(-1,-(n.y/2-1),-1),(n.x/2-1)))
  k2.set = rbind(k2.set.part1,k2.set.part2)
  n.k1 = nrow(k1.set)
  n.k2 = nrow(k2.set)
  
  z.recover.set.1 = array(0, dim=c(n.x*n.y,1))
  alpha.c.1 = array(0/0,dim=c(n.k1,1))
  for (i in 1:n.k1){
    case1 = which( abs(my.x-k1.set[i,1])<0.001 )
    case2 = which( abs(my.y-k1.set[i,2])<0.001 )
    alpha.c.1[i,1]= Re( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.1 = z.recover.set.1+ alpha.c.1[i,1] * f.c
  }
  
  z.recover.set.2 = array(0, dim=c(n.x*n.y,1))
  alpha.c.2 = array(0/0,dim=c(n.k2,1))
  alpha.s.2 = array(0/0,dim=c(n.k2,1))
  for (i in 1:n.k2){
    case1 = which( abs(my.x-k2.set[i,1])<0.001 )
    case2 = which( abs(my.y-k2.set[i,2])<0.001 )
    alpha.c.2[i,1]= Re( my.s[case1,case2] )
    alpha.s.2[i,1]= -Im( my.s[case1,case2] )
    f.c = cos(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    f.s = sin(my.x[case1]*2*pi*rep(grd.xx,each=n.y)/n.x+my.y[case2]*2*pi*rep(grd.yy,n.x)/n.y)
    z.recover.set.2 = z.recover.set.2 + 2*alpha.c.2[i,1] * f.c + 2*alpha.s.2[i,1] * f.s
  }
  
  z.recover.3 = (z.recover.set.1 + z.recover.set.2)
  #par(mfrow = c(1, 2))
  #rasterImage2(x = grd.x,
  #             y = grd.y,
  #             z = z)
  #rasterImage2(x = grd.x,
  #             y = grd.y,
  #             z = (t(matrix((z.recover.3),nrow=n.y))))
  
  alpha = c(alpha.c.1, alpha.c.2, alpha.s.2)
  Alpha[,i.t] = alpha
  
  rasterImage2(x = grd.x,
               y = grd.y,
               z = z, main=paste("time:", i.t, sep=""),xlab="",ylab="",palette=colPalette1)
  #rasterImage2(x = my.x,
              # y = my.y,
              # z = round((matrix(Re(my.s),nrow=n.x,byrow = FALSE)),2),
              # xlab=expression('waver number, k'[1]),
              # ylab=expression('waver number, k'[2]),
              # main="real part")
 # rasterImage2(x = my.x,
             #  y = my.y,
             #  z = round((matrix(Im(my.s),nrow=n.x,byrow = FALSE))),
             #  xlab=expression('waver number, k'[1]),
             #  ylab=expression('waver number, k'[2]),
             #  main="imaginary part")
  
}
dev.off()


#******************************************************
#******************************************************
# Save the output, which is to be used by SPDE_Analysis_e3
#******************************************************
#******************************************************
if (FALSE){
  save.image("data/data_0313_example3_process.RData")
}







