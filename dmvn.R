#******************************************************
#******************************************************
#******************************************************
# Paper: Statistical Modeling for Spatio-Temporal Data from Physical Convection-Diffusion Processes

# ------------------------------------------------------------------------------------------
# Description (See "readme.docx" for more details): 
# Subroutines for SPDE_Analylsis_e3.R
# modified dmvnorm from its source code
# ------------------------------------------------------------------------------------------
# Authors: XXX
# Last revision: 06/01/2019
#******************************************************
#******************************************************
#******************************************************

dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE)
{
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if(!missing(mean)) {
    if(!is.null(dim(mean))) dim(mean) <- NULL
    if (length(mean) != p)
      stop("mean and sigma have non-conforming size")
  }
  if(!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = 1000*sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("sigma must be a symmetric matrix (modified)")
  }
  
  ## <faster code contributed by Matteo Fasiolo mf364 at bath.ac.uk
  dec <- tryCatch(chol(sigma), error=function(e)e)
  if (inherits(dec, "error")) {
    ## warning("cannot compute chol(sigma)"); return(NaN)
    ## behave the same as dnorm(): return Inf or 0
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf # and all other f(.) == 0
  } else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp ^ 2)
    logretval <- - sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if(log) logretval else exp(logretval)
}