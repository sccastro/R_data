################################################################################
# Uniform/Lognormal B = b-A parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,k=2,st0=.05)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with b, A, k, s, ter 
  p$B <- p$b-p$A
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p$b = p$A + p$B
  p
}

