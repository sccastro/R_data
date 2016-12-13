################################################################################
# RD for coherance (CO) experiments, 
# assumes factor CO (with % entries) and S with stimulus (e.g. left or right)
###  drift = vi + vs * (CO), vi(first level of S) = -vi(second level of S)
################################################################################

# bound.tol=.01
start.model <- function(dati,pnams,Bounds,bound.tol=.01) {
  p <- sfun(dati,pnams)
  Di=attributes(dati)$D
  vs <- function(y,x) {lm("y~x-1",data.frame(y=y,x=x))$coefficients[1]}
  CO <- rep(as.numeric(as.character(levels(Di$D$CO))),each=2) # 2 accumulators
  db <- diff(Bounds$vs)*bound.tol
  s=pmin(pmax(tapply(p$v,Di$D[,Di$F[Di$F!="CO"]],vs,x=CO),
              Bounds$vs[1]+db),Bounds$vs[2]-db)
  p$vi <- rep(0,dim(p)[1])
  p$vs <- rep(s,each=length(CO))
  p
}

M2P <- function(p,Di=NULL) {
  # p is a data frame (made by start.model) with c("a","v","ter","z","sz","sv","st") 
  p$Z <- p$z/p$a                                     # z as a proprotion of a
  p$SZ <- p$sz/(2*p$a*apply(cbind(p$Z,1-p$Z),1,min))  # sz as a proprtion of min(z,a-z)
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # on natural scale change to pfun and dfun parameters
  p$z <- p$Z*p$a
  p$sz <- 2*p$SZ*p$a*apply(cbind(p$Z,1-p$Z),1,min)
  
  # USE COHERANCE TO GET DRIFT RATE
  is.left <- D$D$S==levels(D$D$S)[1] # e.g., "left"
  p$v[is.left] <- p$vi[is.left]+
    p$vs[is.left]*as.numeric(as.character(D$D$CO[is.left]))
  p$v[!is.left] <- -p$vi[!is.left]+
    p$vs[!is.left]*as.numeric(as.character(D$D$CO[!is.left]))
  # for accuracte responding
  #  stimulus corresponding to first response level has negative v
  #  stimulus corresponding to second response level has positive v
  # as always pass row correspodning to first response level flip sign
  # for correct responses
  p$v[as.logical(D$D[[D$SC[1]]])] <- -p$v[as.logical(D$D[[D$SC[1]]])] 
  p
}


