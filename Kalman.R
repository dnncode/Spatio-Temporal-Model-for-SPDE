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
#data.alpha = Alpha[,1:6]
#DSTM = dlm(m0=m0, C0=C0, FF=diag(n.alpha), GG=GG, W=S_W, V=s_V^2*diag(n.alpha))
#Filt = dlmFilter(data.alpha, DSTM)


Kalman=function(data.alpha, m0, C0, GG, FF, s_W, s_V){
  #data.alpha = Alpha[,1:6]
  # Initial Condition
  #m0 = Alpha[,1]
  #C0 = 0.01^2*diag(n.alpha)
  #GG = as.matrix( expm(G*0.1) )
  
  nn.t = ncol(data.alpha)
  n_alpha_1 = length(m0)
  n_alpha_2 = nrow(data.alpha)
  
  state.filt = array(0/0, dim=c(n_alpha_1,nn.t+1)) 
  a = array(0/0, dim=c(n_alpha_1,nn.t)) 
  f = e = array(0/0, dim=c(n_alpha_2,nn.t)) 
  state.filt[,1] = m0
  C = C0
  W = s_W
  #W = round(W,2)
  #FF = diag(n_alpha)
  V=s_V^2*diag(n_alpha_2)
  
  R = Q = C = vector("list", length = nn.t)
  a[,1] = GG %*% state.filt[,1]
  R[[1]] = GG %*% C0 %*% t(GG) + W
  #R[[1]] = round(R[[1]],2)
  
  f[,1] = FF %*% a[,1]
  Q[[1]] = FF%*%R[[1]]%*%t(FF) + V
  
  e[,1] = data.alpha[,1] - f[,1]
  #cma <- chol(Q[[1]],tol=1e-4)
  #inv.Q = chol2inv(cma)
  #inv.Q = round(inv.Q,4)
  inv.Q = solve(Q[[1]])
  
  state.filt[,2] = a[,1] + R[[1]]%*%t(FF)%*%inv.Q%*%e[,1]
  #C[[1]] = R[[1]] - R[[1]]%*% t(FF)%*%inv.Q%*%FF %*%R[[1]]
  
  
  tmp1 = chol(R[[1]])
  tmp2 = chol(inv.Q)
  C[[1]] = t(tmp1)%*%tmp1 - t(tmp1) %*% tmp1 %*% t(FF)%*%t(tmp2)%*%tmp2 %*%FF %*% t(tmp1)  %*%tmp1

  #C[[1]] = round(C[[1]],4)
  #isSymmetric(C[[1]])
  #isSymmetric((tmp1 %*% t(FF)%*%t(tmp2)%*%tmp2 %*%FF %*% t(tmp1) ) )
  #isSymmetric(R[[1]])
  #isSymmetric(Q[[1]])
  #isSymmetric(inv.Q)
  #isSymmetric(V)
  
  
  for (i in 2:nn.t){
    a[,i] = GG %*% state.filt[,i]
    R[[i]] = GG %*% C[[i-1]] %*% t(GG) + W
    
    f[,i] = FF %*% a[,i]
    Q[[i]] = FF%*%R[[i]]%*%t(FF) + V
    
    e[,i] = data.alpha[,i] - f[,i]
    #cma <- chol(Q[[i]])
    #inv.Q = chol2inv(cma)
    inv.Q = solve(Q[[i]])
    state.filt[,i+1] = a[,i] + R[[i]]%*%t(FF)%*%inv.Q%*%e[,i]
    C[[i]]= R[[i]] - R[[i]]%*%t(FF)%*%inv.Q%*%FF%*%R[[i]]
    
    # process the digits:
    # C[[i]] = round(C[[i]],1)
    #isSymmetric(W)
    #isSymmetric(R[[i]])
    #isSymmetric(Q[[i]])
    #isSymmetric(inv.Q)
    #isSymmetric(C[[i]])
  }
  
  
  
  
  
  
  Kalman.list = vector("list", length = 6)
  Kalman.list[[1]] = a
  Kalman.list[[2]] = R
  Kalman.list[[3]] = f
  Kalman.list[[4]] = Q
  Kalman.list[[5]] = state.filt
  Kalman.list[[6]] = C
  return(Kalman.list) # a, R, f, Q, m, C  
}
