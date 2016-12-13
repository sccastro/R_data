################################################################################
# Original LBA a=b/A, v(Correct)=1-v(Error), sv free parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,unitv=T,st0=.05)
}

M2P <- function(p) {
  # starts is a data frame (made by start.model) with b, A, v, sv, ter 
  p$a <- p$A/p$b
  p
}

P2Mfit <- function(p,pmap,D) {
  # flip sign for error on logit scale to get 1-v on natural 
  p$v <- c(p$v*(2*as.numeric(as.logical(D$D[[D$SC[1]]]))-1))
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p$A = p$a*p$b
  p
}
