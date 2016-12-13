################################################################################
# LNR seperate mu for each accumulator parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.5)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with m, v, r, ter, pc 
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p
}

