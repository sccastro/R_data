################################################################################
# LNR seperate mu for each accumulator parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  p<- sfun(dati,pnams,pc=.01,aprior=.5,tertune=.5)
  p$a <- min(p$v)/2
  tmp <- cbind.data.frame(y=p$v-p$a,x=exp(p$m))
  p$c <- -log(max(c(coef(lm(y~x-1,tmp)),.0001)))
  p$e <- p$c-p$m
  p$u <- 1
  p
}

M2P <- function(p) {
  # called by sfun to change to parameters specified by model.rc
  # not required in this case
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  es <- p$e/p$u # mean evidence scaled for urgency
  p$m <- p$c-es
  p$v <- p$a + exp(-es)
  p
}

