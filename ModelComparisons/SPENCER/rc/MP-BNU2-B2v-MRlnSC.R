################################################################################
# LOGNROMAL NON-DECISION TIME, ANGLE EFFECTS SCALE 
# SPECIFIC TO MR PROJECT ASSUMING $A ANGLE FACTOR and $AS parameter (rate meanlog)
# $st0 is rate sdlog
# Original LBA B = b-A, sv=1, parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  p <- sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,unitv=F,st0=.05)
  # gamma angle slope (AS)
  p$AS <- -1 # scale about .35s effect at 180
  p$AI <- -2 # mirror intercept for shift, about effect at .13s
  p
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with b, A, v, sv, ter 
  p$B <- p$b-p$A
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p$b <- p$A + p$B
  p$scale <- ifelse(D$D$A=="0",-Inf,-p$AS+log(as.numeric(as.character(D$D$A))/180))
  if (any(names(p)=="AI"))
    p$scale[D$D$S=="mirror"] <- p$scale[D$D$S=="mirror"] + p$AI[D$D$S=="mirror"] 
  p
}

