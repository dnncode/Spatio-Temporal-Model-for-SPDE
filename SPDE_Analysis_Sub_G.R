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

if (FALSE){ # values for testing the G.function
  wd = "C:/Users/xl027/Desktop/25.SPDE_SS/code/20180814"
  setwd(wd)
  source("SPDE_Analysis_Sub_G_Sub.R")
  
  N.alpha = n.alpha.l
  N.k1 = n.k1.l
  N.k2 = n.k2.l
  Omega.1.k1 = omega.1.k1.l
  Omega.2.k1 = omega.2.k1.l
  Omega.1.k2 = omega.1.k2.l
  Omega.2.k2 = omega.2.k2.l
  Phi.c.mat.k1 = phi.c.mat.k1.l
  Phi.s.mat.k1 = phi.s.mat.k1.l
  Phi.c.mat.k2 = phi.c.mat.k2.l
  Phi.s.mat.k2 = phi.s.mat.k2.l
  N.x = n.x.l
  N.y = n.y.l
  Delta.x = delta.x.l
  Delta.y = delta.y.l
  
  MU = mu
  SIGMA =sigma
  DIV.SIGMA = div.sigma
  
  ZETA=zeta.G
  Par=FALSE
}




G.function = function(N.alpha, N.k1, N.k2,
                      Omega.1.k1, Omega.2.k1, Omega.1.k2, Omega.2.k2,
                      Phi.c.mat.k1, Phi.s.mat.k1, Phi.c.mat.k2, Phi.s.mat.k2,
                      N.x, N.y, Delta.x, Delta.y,
                      MU, SIGMA, DIV.SIGMA, ZETA,
                      Par){  
                      # Par==TRUE, parallel computing using foreach

  C = 1/2 # the constant term
  
  if (Par){  # parallel computing for G
    
    #******************************************************
    # the transition matrix g related to MU (convection)
    #******************************************************
    g.mat.1 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to MU 
    #tmp.s = proc.time()
    dummy = foreach(ll=1:4) %dopar% {
      aa = sapply(1:N.k1, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k1, base2=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb= 2* sapply(1:N.k2, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.c.mat.k1, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = -2*sapply(1:N.k2, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k1, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      return(c(aa,bb,cc))
    }
    g.mat.1[1:4,] = do.call(rbind,dummy)
    dummy = foreach(ll = (N.k1+N.k2+1):N.alpha) %dopar% {
      aa = sapply(1:N.k1, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.s.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.s.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = - 2*sapply(1:N.k2, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                      base1=Phi.s.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.1[(N.k1+N.k2+1):N.alpha,] = do.call(rbind,dummy)
    dummy = foreach(ll = 5:(N.k1+N.k2)) %dopar% {
      aa = sapply(1:N.k1, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k2, base2=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = - 2*sapply(1:N.k2, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                      base1=Phi.c.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.1[5:(N.k1+N.k2),] = do.call(rbind,dummy)
    #print(proc.time()-tmp.s)
    
    
    
    #******************************************************
    # the transition matrix g related to Sigma(diffusion)
    #******************************************************
    g.mat.2 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to mu 
    
    dummy = foreach(ll=1:4) %dopar% {
      aa = sapply(1:N.k1, My.fun2,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k1, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb = 2*sapply(1:N.k2, My.fun2,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.c.mat.k1, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = 2*sapply(1:N.k2, My.fun3,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k1, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      return(c(aa,bb,cc))
    }
    g.mat.2[1:4,] = do.call(rbind,dummy)
    
    dummy = foreach(ll = (N.k1+N.k2+1):N.alpha) %dopar% {
      aa= sapply(1:N.k1, My.fun2,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                 base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)  / (2*C)
      bb= 2* sapply(1:N.k2, My.fun2,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      cc = 2*sapply(1:N.k2, My.fun3,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.2[(N.k1+N.k2+1):N.alpha,] = do.call(rbind,dummy)
    
    dummy = foreach(ll = 5:(N.k1+N.k2)) %dopar% {
      aa = sapply(1:N.k1, My.fun2,l=(ll-N.k1), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun2,l=(ll-N.k1),Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc =  2*sapply(1:N.k2, My.fun3,l=(ll-N.k1), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.2[5:(N.k1+N.k2),] = do.call(rbind,dummy)
    
    
    #******************************************************
    # the transition matrix g related to Zeta
    #******************************************************
    g.mat.3 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to mu 
    
    dummy = foreach(ll=1:4) %dopar% {
      aa = sapply(1:N.k1, My.fun4,l=ll, Zeta = ZETA,  
                  base1=Phi.c.mat.k1, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb = 2*sapply(1:N.k2, My.fun4,l=ll, Zeta = ZETA,
                    base1=Phi.c.mat.k1, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = 2*sapply(1:N.k2, My.fun4,l=ll, Zeta = ZETA,
                    base1=Phi.c.mat.k1, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      return(c(aa,bb,cc))
    }
    g.mat.3[1:4,] = do.call(rbind,dummy)
    
    dummy = foreach(ll = (N.k1+N.k2+1):N.alpha) %dopar% {
      aa = sapply(1:N.k1, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                  base1=Phi.s.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                     base1=Phi.s.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc= 2*sapply(1:N.k2, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                   base1=Phi.s.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.3[(N.k1+N.k2+1):N.alpha,] = do.call(rbind,dummy)
    
    dummy = foreach(ll = 5:(N.k1+N.k2)) %dopar% {
      aa = sapply(1:N.k1, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                  base1=Phi.c.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                     base1=Phi.c.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = 2*sapply(1:N.k2, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                    base1=Phi.c.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      return(c(aa,bb,cc))
    }
    g.mat.3[5:(N.k1+N.k2),] = do.call(rbind,dummy)
    
    
  }else{ # non-parallel for loop
    
    #******************************************************
    # the transition matrix g related to MU (convection)
    #******************************************************
    g.mat.1 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to MU 
    #tmp.s = proc.time()
    tmp = array(0/0, dim=c(4, N.alpha))
    i = 1
    for (ll in c(1:4)){
      aa = sapply(1:N.k1, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k1, base2=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb= 2* sapply(1:N.k2, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.c.mat.k1, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = -2*sapply(1:N.k2, My.fun1,l=ll, Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k1, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.1[1:4,] = tmp
    
    tmp = array(0/0, dim=c(length(c((N.k1+N.k2+1):N.alpha)), N.alpha))
    i = 1
    for (ll in c((N.k1+N.k2+1):N.alpha)){
      aa = sapply(1:N.k1, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.s.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.s.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = - 2*sapply(1:N.k2, My.fun1,l=(ll-N.k1-N.k2), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                      base1=Phi.s.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.1[(N.k1+N.k2+1):N.alpha,] = tmp
    
    tmp = array(0/0, dim=c(length(c(5:(N.k1+N.k2))), N.alpha))
    i = 1
    for (ll in c(5:(N.k1+N.k2))){
      aa = sapply(1:N.k1, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k2, base2=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = - 2*sapply(1:N.k2, My.fun1,l=(ll-N.k1), Mu=MU, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                      base1=Phi.c.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.1[5:(N.k1+N.k2),] = tmp
    #print(proc.time()-tmp.s)
    
    
    
    #******************************************************
    # the transition matrix g related to Sigma(diffusion)
    #******************************************************
    g.mat.2 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to mu 
    
    tmp = array(0/0, dim=c(4, N.alpha))
    i = 1
    for (ll in c(1:4)){
      aa = sapply(1:N.k1, My.fun2,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k1, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb = 2*sapply(1:N.k2, My.fun2,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.c.mat.k1, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = 2*sapply(1:N.k2, My.fun3,l=ll, Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k1, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.2[1:4,] = tmp
    
    tmp = array(0/0, dim=c(length(c((N.k1+N.k2+1):N.alpha)), N.alpha))
    i = 1
    for (ll in c((N.k1+N.k2+1):N.alpha)){
      aa= sapply(1:N.k1, My.fun2,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                 base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)  / (2*C)
      bb= 2* sapply(1:N.k2, My.fun2,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      cc = 2*sapply(1:N.k2, My.fun3,l=(ll-N.k1-N.k2), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                    base1=Phi.s.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.2[(N.k1+N.k2+1):N.alpha,] = tmp
    
    tmp = array(0/0, dim=c(length(c(5:(N.k1+N.k2))), N.alpha))
    i = 1
    for (ll in c(5:(N.k1+N.k2))){
      aa = sapply(1:N.k1, My.fun2,l=(ll-N.k1), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k1, Ky=Omega.2.k1, 
                  base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k1, base2.s=Phi.s.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun2,l=(ll-N.k1),Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc =  2*sapply(1:N.k2, My.fun3,l=(ll-N.k1), Sigma=SIGMA, Div.Sigma=DIV.SIGMA, Kx=Omega.1.k2, Ky=Omega.2.k2, 
                     base1=Phi.c.mat.k2, base2.c=Phi.c.mat.k2, base2.s=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.2[5:(N.k1+N.k2),] = tmp
    
    
    #******************************************************
    # the transition matrix g related to Zeta
    #******************************************************
    g.mat.3 = array(0/0, dim=c(N.alpha, N.alpha))  # the transition matrix g related to mu 
    
    tmp = array(0/0, dim=c(4, N.alpha))
    i = 1
    for (ll in c(1:4)){
      aa = sapply(1:N.k1, My.fun4,l=ll, Zeta = ZETA,  
                  base1=Phi.c.mat.k1, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      bb = 2*sapply(1:N.k2, My.fun4,l=ll, Zeta = ZETA,
                    base1=Phi.c.mat.k1, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      cc = 2*sapply(1:N.k2, My.fun4,l=ll, Zeta = ZETA,
                    base1=Phi.c.mat.k1, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / C
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.3[1:4,] = tmp
    
    tmp = array(0/0, dim=c(length(c((N.k1+N.k2+1):N.alpha)), N.alpha))
    i = 1
    for (ll in c((N.k1+N.k2+1):N.alpha)){
      aa = sapply(1:N.k1, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                  base1=Phi.s.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                     base1=Phi.s.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc= 2*sapply(1:N.k2, My.fun4,l=(ll-N.k1-N.k2), Zeta = ZETA,
                   base1=Phi.s.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.3[(N.k1+N.k2+1):N.alpha,] = tmp
    
    tmp = array(0/0, dim=c(length(c(5:(N.k1+N.k2))), N.alpha))
    i = 1
    for (ll in c(5:(N.k1+N.k2))){
      aa = sapply(1:N.k1, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                  base1=Phi.c.mat.k2, base2=Phi.c.mat.k1, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y) / (2*C)
      bb = 2* sapply(1:N.k2, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                     base1=Phi.c.mat.k2, base2=Phi.c.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      cc = 2*sapply(1:N.k2, My.fun4,l=(ll-N.k1), Zeta = ZETA,
                    base1=Phi.c.mat.k2, base2=Phi.s.mat.k2, N=(N.x*N.y), Ds.x = Delta.x, Ds.y=Delta.y)/ (2*C)
      tmp[i,]=c(aa,bb,cc)
      i = i + 1
    }
    g.mat.3[5:(N.k1+N.k2),] = tmp
    
  }
  G.1 = g.mat.1+g.mat.2+g.mat.3
  
 
  
  G=G.1
  return(G)
}
  
