################################################################################
# GAMMA NON-DECISION TIME, ANGLE EFFECTS SCALE 
# SPECIFIC TO MR PROJECT ASSUMING $A ANGLE FACTOR and $AS parameter
# Original LBA B = b-A, sv=1, parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  p <- sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,unitv=F,st0=.05)
  # gamma angle slope (AS)
  p$AS <- 1   # scale = 1 at 180
  p$AI <- 0.1 # mirror intercept for shift
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
  p$scale <- as.numeric(as.character(D$D$A))*p$AS/180
  if (any(names(p)=="AI"))
    p$scale[D$D$S=="mirror"] <- p$scale[D$D$S=="mirror"] + p$AI[D$D$S=="mirror"] 
  p
}

