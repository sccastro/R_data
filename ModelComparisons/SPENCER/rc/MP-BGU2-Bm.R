################################################################################
# Uniform/Gamma B = b-A parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,K=2,st0=.05)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with b, A, k, s, ter 
  p$B <- p$b-p$A
  p$v <- p$k*p$s
  p$sv <- sqrt(p$k)*p$s
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p$b <- p$A + p$B
  p$k <- (p$v/p$sv)^2
  p$s <- ((p$sv)^2)/p$v
  p
}

