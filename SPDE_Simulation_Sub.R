#******************************************************
#******************************************************
#******************************************************
# Paper: Statistical Modeling for Spatio-Temporal Data from Physical Convection-Diffusion Processes

# ------------------------------------------------------------------------------------------
# Description (See "readme.docx" for more details): 
# Subroutines for SPDE_Simulation_e3.R
# ------------------------------------------------------------------------------------------
# Authors: XXX
# Last revision: 06/01/2019
#******************************************************
#******************************************************
#******************************************************

# Create Grid
grid.create = function(n.x,n.y,PLOT){
  
  # create grid
  #seed.x <- seq(-(n/2-1),n/2,1)
  seed.x <- seq(0,1-1/n.x,1/n.x)
  #seed.y <- seq(-(n/2-1),n/2,1)
  seed.y <- seq(0,1-1/n.y,1/n.y)
  grd.out.x <- rep( seed.x , each=length(seed.y) )
  grd.out.y <- rep( seed.y , length(seed.x))
  
  #min.x = min( grd.out.x)
  #max.x = max( grd.out.x)
  #range.x = max.x - min.x
  #min.y = min( grd.out.y)
  #max.y = max( grd.out.y)
  #range.y = max.y - min.y

  #grd.out.x = (grd.out.x-min.x)/range.x
  #grd.out.y = (grd.out.y-min.y)/range.y
  
  grd.out <- SpatialPoints( cbind(grd.out.x,grd.out.y) )
  #plot(grd.out,axes=TRUE)
  
  if (PLOT){
    plot(grd.out,axes=TRUE,col="blue")
  }
  
  return(grd.out)
}


# Calculation of the cos and sin basis
phi.c = function(i, Kx, Ky, Grd){
  tmp = matrix( coordinates(Grd), ncol=2 ) %*% matrix( c(Kx[i], Ky[i]), nrow=2) 
  output = round( apply(tmp, 1, cos), 4)
  return(output)
}

phi.s = function(i, Kx, Ky, Grd){
  tmp = matrix( coordinates(Grd), ncol=2 ) %*% matrix( c(Kx[i], Ky[i]), nrow=2) 
  output = round( apply(tmp, 1, sin), 4)
  return(output)
}

# generate of mu and sigma
mu_fun = function(i, x.Coords, y.Coords){
  
  if (FALSE){
    tmp = 1 * exp( -2* sqrt( x.Coords[i]^2 + y.Coords[i]^2 ) )
    output = c(tmp, tmp)
  }
  
  if (FALSE){
    mu_x = 1*exp(  -x.Coords[i] )
    mu_y = 1*exp(  -(1-y.Coords[i])*2 )
    output = c(mu_x, mu_y)
  }
  
  if (FALSE){
    mu_x = 0.1    
    mu_y = 0.1
    output = c(mu_x, mu_y)
  }
  
  
  if (TRUE){  # for the numerical example 1
    mu_x = 0.1 - 0.1*x.Coords[i]
    mu_y = 0 + 0.1*y.Coords[i]
    output = c(mu_x, mu_y)
  }
  output
  
}

s_fun = function(i, x.Coords, y.Coords, Mu){
  
  if (TRUE){
    output =  matrix(c(0.1,0,0,0.1),ncol=2)
  }
  output
  
}


# sub-functions required for computing the transition matrix:
My.fun1 = function(j, l, Mu, Kx, Ky, base1, base2, N=(n.x*n.y), Ds.x = delta.x, Ds.y=delta.y){
  k = matrix(c(Kx[j],Ky[j]),ncol=1)
  term1 = Mu %*% k
  term2 = base2[,j] * base1[,l]
  output = sum( term1 * term2 * Ds.x * Ds.y )
  return(output)
}


My.fun2 = function(j, l, Sigma, Div.Sigma, Kx, Ky, base1, base2.c, base2.s, N=(n.x*n.y), Ds.x = delta.x, Ds.y=delta.y){
  k = matrix(c(Kx[j],Ky[j]),ncol=1)
  # for a given s
  term3 = array(0/0, dim=c(sqrt(N),1))
  for (ss in 1:sqrt(N)){
    term1 = - t(k) %*% Sigma[[ss]] %*% (k) * base2.c[ss,j]
    term2 = - Div.Sigma[[ss]] %*% (k) * base2.s[ss,j]
    term3[ss,1] = (term1 + term2) * base1[ss,l] 
  }
  output = sum( term3 * Ds.x * Ds.y )
  return(output)
}

My.fun3 = function(j, l, Sigma, Div.Sigma, Kx, Ky, base1, base2.c, base2.s, N=(n.x*n.y), Ds.x = delta.x, Ds.y=delta.y){
  k = matrix(c(Kx[j],Ky[j]),ncol=1)
  # for a given s
  term3 = array(0/0, dim=c(sqrt(N),1))
  for (ss in 1:sqrt(N)){
    term1 = - t(k) %*% Sigma[[ss]] %*% (k) * base2.s[ss,j]
    term2 = - Div.Sigma[[ss]] %*% (k) * base2.c[ss,j]
    term3[ss,1] = (term1 + term2) * base1[ss,l] 
  }
  output = sum( term3 * Ds.x * Ds.y )
  return(output)
}

My.fun4 = function(j, l, Zeta, base1, base2, N=(n.x*n.y), Ds.x = delta.x, Ds.y=delta.y){
  term1 = Zeta
  term2 = base2[,j] * base1[,l]
  output = -sum( term1 * term2 * Ds.x * Ds.y )
  return(output)
}


# used for mu and sigma: approach 2
convection_fun_1 = function(X.input, gamma_x, gamma_y){
  mu_x = X.input %*% gamma_x
  mu_y = X.input %*% gamma_y
  output = cbind(mu_x, mu_y)
  output
}
if (FALSE){
  velocity = convection_fun_1(X.input=X, gamma_x=rnorm(J*3,0,0.1),gamma_y=rnorm(J*3,0,0.1))
  plot(grd.original)
  arrows(coordinates(grd.original)[,1],coordinates(grd.original)[,2],
         coordinates(grd.original)[,1]+velocity[,1]/2,
         coordinates(grd.original)[,2]+velocity[,2]/2,length=0.1,col="darkgreen",lwd=1)
}
diffusion_fun = function(i, x.Coords, y.Coords, Mu, r){
  if (TRUE){ # diffusion formulation I
    output =  matrix(c(r,0,0,r),ncol=2)
  }
  output
}

