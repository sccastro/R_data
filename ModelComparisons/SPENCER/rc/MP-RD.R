################################################################################
# Original RD
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams)
}

M2P <- function(p,Di=NULL) {
  # p is a data frame (made by start.model) with c("a","v","ter","z","sz","sv","st") 
  p$Z <- p$z/p$a                                     # z as a proprotion of a
  p$SZ <- p$sz/(2*p$a*apply(cbind(p$Z,1-p$Z),1,min))  # sz as a proprtion of min(z,a-z)
  p
}


P2Mfit <- function(p,pmap,D) {
  # for accuracte responding
  #  stimulus corresponding to first response level has negative v
  #  stimulus corresponding to second response level has positive v
  # as always pass row correspodning to first response level flip sign
  # for correct responses
  p$v[as.logical(D$D[[D$SC[1]]])] <- -p$v[as.logical(D$D[[D$SC[1]]])]
  p
}

P2Mnatural <- function(p,pmap,D) {
  # on natural scale change to pfun and dfun parameters
  p$z <- p$Z*p$a
  p$sz <- 2*p$SZ*p$a*apply(cbind(p$Z,1-p$Z),1,min)
  p
}


