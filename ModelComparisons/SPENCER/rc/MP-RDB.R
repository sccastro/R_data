################################################################################
# RD with seperate LBA like "B" parameter 
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with c("a","v","ter","z","sz","sv","st") 
  sz2 <- p$sz/2
  p$B0 <- p$z - sz2
  p$Ba <- p$a - p$z - sz2
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
  sz2 <- p$sz/2
  p$z <- p$B0+sz2
  p$a <- p$z+sz2+p$Ba
  p
}


