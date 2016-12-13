################################################################################
# Uniform/Lognormal B = b-A and v/sv parameterization 
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.5)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with m, v, r, ter, pc 
  S=log(1+(p$s^2)/exp(2*log(p$m)))
  p$v <- -S/2+log(p$m)
  p$sv <- sqrt(S)
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  S <- p$sv^2
  p$m <- exp(p$v+S/2)
  p$s <- sqrt((exp(S)-1)*exp(2*p$v+S))
  p
}

## Expressions for mean and variance of lognormal, based on underlying
## normal distribution with parameters mu and sigma:
#mean = exp(mu + 0.5*sigma^2)
#var = (exp(sigma^2)-1)*exp(2mu+sigma^2)
#
## sigma appears only as sigma^2, so clean that up for ease:
#Let S=sigma^2
#mean = exp(mu + 0.5*S)
#var = (exp(S)-1)*exp(2mu+S)
#
## Re-arrange to get mu as a fucntion of S and observed mean.
#log(mean) = mu + S/2
#mu = -S/2+log(mean)
#
## Now take the expression for observed variance, re-arrange:
#var = (exp(S)-1)*exp(-S+2*log(mean)+S)
#var = exp(S)*exp(-S + 2*log(mean)+S) - exp(-S+2*log(mean)+S)
#var = exp(S)*exp(2*log(mean)) - exp(2*log(mean))
#var = exp(2*log(mean))*(exp(S)-1)
#1+var/exp(2*log(mean)) = exp(S)
#S=log(1+var/exp(2*log(mean)))
