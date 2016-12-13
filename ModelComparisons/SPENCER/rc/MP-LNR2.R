################################################################################
# LNR variance linear in mean parameterizaton
################################################################################

start.model <- function(dati,pnams,Bounds) {
  p<- sfun(dati,pnams,pc=.01,aprior=.5,tertune=.5)
  p$a <- min(p$v)/2                        # mean ~ variance intercept
  tmp <- cbind.data.frame(y=p$v-p$a,x=exp(p$m))
  p$b <- max(c(coef(lm(y~x-1,tmp)),.0001)) # mean ~ variance slope
  p$e <- -p$m  # log evidence rate
  p$c <- 0     # log criterion
  p$d <- 1     # ???
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
  p$m <- p$c-p$e*p$d
  p$v <- p$a + p$b*exp(p$m)
  p
}

