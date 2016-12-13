################################################################################
# ARCHITECTURE = 1 Binary ACCUMULATOR
# Rogers contamination model
#
#
# called by augment.design in mt.R to setup fast mapping of 
# parameters by fitp2modelp
# DUMMY CALL IN THIS ARCHITECTURE
get.reorder <- function(D,Cpars) {
  list(P=NULL,D=NULL,N=NULL)
}


fitp2modelp <- function(p,pmap,D) {
  # takes fitting parameterization (untransformed) and changes
  # to parameters used by dfun, pfun and rfun etc.
  
  # expand paramater vector to parameter list
  p <- mapfitp(p,pmap,trans=F)
  # parameterization specific transform on fitting scale
  p <- P2Mfit(p,pmap,D)
  # transform to natural scale
  p <- fitp2p(pmap,p)
  # parameterization specific transform on natural scale
  p <- P2Mnatural(p,pmap,D)
  # add in contamination if not in model parameterization
  # also make sure no zeros occur
  if (!any(names(p)=="pc")) 
    p[["pc"]] <- rep(1e-10,length(p[[1]])) else
    p$pc[p$pc<1e-10] <- 1e-10
  if (!any(names(p)=="pg")) p[["pg"]] <- 0.5 else 
    p$pg[p$pg<1e-10] <- 1e-10
  data.frame(p)
}


  
# p=strt
objective <- function(p,dat,pmap,pu,qmpe=F,precision=2.5) {
  D <- attr(dat,"D")
#   cat("\n")
#   print(p["sv.I"])   
  p <- fitp2modelp(p,pmap,D)
#   print(p$sv[1])
  cp <- 1/diff(pu) # contamination density
  if ( !qmpe ) { # Likelihood
#     L <- vector(length=0)
#     dat <- dat[is.finite(dat[,D$RT]),]
#     for (i in levels(D$D$lcell)) {
#       dati <- dat[dat$cell==i,]
#       isin <- D$D$lcell==i
#       if (dim(dati)[1]>0) {
#         L <- c(L, p$pc[isin]*cp*p$pg[isin] + (1-p$pc[isin])*
#           dfun(t=dati[[D$RT]],parlist=p[isin,])) # likelihood
#       }
#     }
#     L[!is.finite(L) | L<0] <- 1e-10
#     mll <- -sum(log(L))
    stop("NOT IMPLEMENTED")
  } else { # qmpe       
    qps <- numeric(dim(dat)[1])
    cells <- matrix(levels(dat$cell),ncol=2,byrow=T)
    row.names(cells) <- levels(dat$rcell)
    for (i in levels(D$D$rcell) ) { # get predicted probs (qps) and contaminate
      isr <- dat$rcell==i
      ts <- list(lower=numeric(0),upper=numeric(0))
      for (j in 1:2)
        if (any(dat$cell==cells[i,j])) ts[[j]] <- 
          dat[dat$cell==cells[i,j] & is.finite(dat[[D$RT]]),D$RT]
        # quantile probabilities under model
      qps[isr] <- unlist(pfun(parlist=p[D$D$rcell==i,],ts=ts,doCDF=F,
        precision=precision),use.names=F)
      # contaminated quantile probabilities
      if ( length(ts$lower)!=0 ) 
        cpl <- diff(c(pu[1],ts$lower,pu[2]))*p[D$D$rcell==i,"pg"] else 
        cpl <- numeric(0)
      if ( length(ts$upper)!=0 ) 
        cpu <- diff(c(pu[1],ts$upper,pu[2]))*(1-p[D$D$rcell==i,"pg"]) else 
        cpu <- numeric(0)
      qps[isr] <- 
        p$pc[D$D$rcell==i]*c(cpl,cpu)*cp + (1-p$pc[D$D$rcell==i])*qps[isr]      
    }
    mll <- -sum(dat$qn*log(qps))
  }
  mll # minus log-likelihood
}

# trial by trail qps not yet implemented
simdat <- function(best,pmap,dat,sim.mult=100,maxsamp=1e6,qps=c(.1,.3,.5,.7,.9)) {
  
  qs <- function(x){c(1:length(x))/(length(x)+1)}
  
  D <- attr(dat,"D")
  cnam <- D$SC[1]
  rtnam <- D$RT
  rcnams <- levels(D$D$rcell)
  nc <- D$nc
  tlohi=as.numeric(dimnames(nc)$minmax)
  cp <- 1/diff(tlohi) 
  p <- fitp2modelp(best,pmap,D)
#   if (exists("dfun")) { # likelihood may not be available for some models
#     # add in model and contaminant likelihood
#     dat$lm <- rep(NA,dim(dat)[1])
#     # following not changed if no comtaminant model
#     dat$lc <- dat$lm
#     # add likelihoods to dat
#     for (i in levels(dat$cell) ) {
#       isrt <- dat$cell==i & is.finite(dat[,D$RT])
#       if (any(isrt)) {
#         dat$lm[isrt] <- dfun(ts=dat[isrt,D$RT],parlist=p[i,])
#         if (any(names(p)=="pc")) { # contaminant model
#           dat$lc[isrt] <- p[i,"pc"]*cp*p[i,"pg"]
#           dat$lm[isrt] <- dat$lm[isrt]*(1-p[i,"pc"])
#         }  
#       }
#     }
#   }
  dat$lc <- NA
  dat$lm <- NA
  # get prdicted RT via simulation
  srtnam <- paste(rtnam,"sim",sep=".")
  dat[[srtnam]] <- dat[[rtnam]]
  sqpnam <- "qp.sim"
  dat[[sqpnam]] <- dat$qp
  isrt <- is.finite(dat[[rtnam]])  
  sim <- vector(mode="list",length=length(levels(D$D$rcell)))
  names(sim) <- levels(D$D$rcell)
  for (i in levels(D$D$rcell) ) {
    Di <- D$D[D$D$rcell==i,]
    pi <- p[D$D$rcell==i,][
      c(1:2)[Di[[D$R]]==levels(Di[[D$R]])[1]],,drop=F]
    dati <- dat[dat$rcell==i,]
    dati$cell <- factor(as.character(dati$cell),row.names(Di))
    nsim <- sim.mult*tapply(dati$qn,dati$cell,sum)
    nsim[is.na(nsim)] <- 0
    ns <- nsim; names(ns) <- levels(Di[[D$R]])
    tmp <- rfun(ns,pi,tlohi,maxsamp=maxsamp,sim.mult=sim.mult,roger.contaminate=TRUE)
    sim[[i]] <- tmp$sims
    srt <- data.frame(tmp$rt[!is.na(tmp$rt[,1]),])    
    srt$choice <- factor(srt$choice,levels=0:1,labels=names(tmp$n))
    for (j in 1:length(nsim)) {
      isin1 <- dat$cell==names(nsim)[j] 
      isin2 <- isin1 & is.finite(dat[[rtnam]])
      if ( sum(isin1)==0 | dim(srt)[1]==0 ) {
        dat[isin1,srtnam] <- rep(NA,sum(isin1)) 
        dat[isin1,sqpnam] <- rep(NA,sum(isin1))
      } else {
        dat[isin2,srtnam] <- quantile(srt[srt$choice==names(ns)[j],"rt"],
                                      type=1,probs=dat[isin2,"qp"])
        dat[isin2,sqpnam] <- table(cut(srt[srt$choice==names(ns)[j],"rt"],
                                       c(-Inf,dat[isin1,rtnam]) ))[-(sum(isin1))]/
          length(srt[srt$choice==names(ns)[j],"rt"])
        dat[isin1 &!isin2,sqpnam] <- tmp$n[names(ns)[j]]/sum(tmp$n)
        dat[isin1 &!isin2,sqpnam] <- tmp$n[names(ns)[j]]/sum(tmp$n)
      }
    }
  }
  list(dat=dat[,c("lm","lc",srtnam,sqpnam)],sim=sim) 
}

sim1dat <- function(best,pmap,dat) {
  
  D <- attr(dat,"D")
  cnam <- D$SC[1]
  rtnam <- D$RT
  rcnams <- levels(D$D$rcell)
  nc <- D$nc
  tlohi=as.numeric(dimnames(nc)$minmax)
  cp <- 1/diff(tlohi) 
  p <- fitp2modelp(best,pmap,D)
  # get prdicted RT via simulation
  isrt <- is.finite(dat[[rtnam]])
  sim <- vector(mode="list",length=length(levels(D$D$rcell)))
  names(sim) <- levels(D$D$rcell)
  for (i in levels(D$D$rcell) ) {
    Di <- D$D[D$D$rcell==i,]
    pi <- p[D$D$rcell==i,][
      c(1:2)[Di[[D$R]]==levels(Di[[D$R]])[1]],,drop=F]
    dati <- dat[dat$rcell==i,]
    dati$cell <- factor(as.character(dati$cell),row.names(Di))
    nsim <- tapply(dati$qn,dati$cell,sum)
    nsim[is.na(nsim)] <- 0
    ns <- nsim; names(ns) <- levels(Di[[D$R]])
    sim[[i]] <- 
      rfun(ns,pi,tlohi,maxsamp=sum(ns),onen=sum(ns),sim.mult=1,
           warn=FALSE,roger.contaminate=TRUE)$sims
  }
  sim 
}

