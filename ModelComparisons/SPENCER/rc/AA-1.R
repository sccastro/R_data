################################################################################
# ARCHITECTURE = 1 Binary ACCUMULATOR
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


  
# objective <- function(p,dat,pmap,pu,qmpe=F) {
#   D <- attr(dat,"D")
# #   cat("\n")
# #   print(p["sv.I"])   
#   p <- fitp2modelp(p,pmap,D)
# #   print(p$sv[1])
#   cp <- 1/diff(pu) # contamination density
#   if ( !qmpe ) { # Likelihood
# #     L <- vector(length=0)
# #     dat <- dat[is.finite(dat[,D$RT]),]
# #     for (i in levels(D$D$lcell)) {
# #       dati <- dat[dat$cell==i,]
# #       isin <- D$D$lcell==i
# #       if (dim(dati)[1]>0) {
# #         L <- c(L, p$pc[isin]*cp*p$pg[isin] + (1-p$pc[isin])*
# #           dfun(t=dati[[D$RT]],parlist=p[isin,])) # likelihood
# #       }
# #     }
# #     L[!is.finite(L) | L<0] <- 1e-10
# #     mll <- -sum(log(L))
#     stop("NOT IMPLEMENTED")
#   } else { # qmpe       
#     qps <- numeric(dim(dat)[1])
#     cells <- matrix(levels(dat$cell),ncol=2,byrow=T)
#     row.names(cells) <- levels(dat$rcell)
#     for (i in levels(D$D$rcell) ) { # get predicted probs (qps) and contaminate
#       isr <- dat$rcell==i
#       ts <- list(lower=numeric(0),upper=numeric(0))
#       for (j in 1:2)
#         if (any(dat$cell==cells[i,j])) ts[[j]] <- 
#           dat[dat$cell==cells[i,j] & is.finite(dat[[D$RT]]),D$RT]
#       # quantile probabilities under model
#       qps[isr] <- unlist(pfun(parlist=p[D$D$rcell==i,],ts=ts,doCDF=F),use.names=F)
#       # contaminated quantile probabilities
#       if ( length(ts$lower)!=0 ) 
#         cpl <- diff(c(pu[1],ts$lower,pu[2]))*p[D$D$rcell==i,"pg"] else 
#         cpl <- numeric(0)
#       if ( length(ts$upper)!=0 ) 
#         cpu <- diff(c(pu[1],ts$upper,pu[2]))*(1-p[D$D$rcell==i,"pg"]) else 
#         cpu <- numeric(0)
#       qps[isr] <- 
#         p$pc[D$D$rcell==i]*c(cpl,cpu)*cp + (1-p$pc[D$D$rcell==i])*qps[isr]      
#     }
#     mll <- -sum(dat$qn*log(qps))
#   }
#   mll # minus log-likelihood
# }

objective <- function(p,dat,pmap,pu,qmpe=F,precision=2.5) {
  D <- attr(dat,"D")
  p <- fitp2modelp(p,pmap,D)
  cp <- 1/diff(pu) # contamination density
  if ( !qmpe ) { # Likelihood
    L <- vector(length=0)
    cells <- matrix(levels(dat$cell),ncol=2,byrow=T)
    row.names(cells) <- levels(dat$rcell)
    for (i in levels(D$D$rcell) ) { 
      isr <- dat$rcell==i & is.finite(dat[[D$RT]])
      isin <- D$D$rcell==i
      ts <- list(lower=numeric(0),upper=numeric(0))
      for (j in 1:2)
        if (any(dat$cell==cells[i,j])) ts[[j]] <- 
          dat[dat$cell==cells[i,j] & is.finite(dat[[D$RT]]),D$RT]
      L <- c(L,p$pc[isin][1]*cp*p$pg[isin][1] + (1-p$pc[isin][1])*
        unlist(dfun(parlist=p[isin,],ts=ts,precision=precision),use.names=FALSE))
    }
    L[!is.finite(L) | L<0] <- 1e-10
    mll <- -sum(log(L))
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
      qps[isr] <- unlist(pfun(parlist=p[D$D$rcell==i,],
        ts=ts,precision=precision,doCDF=F),use.names=F)
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
# qps=c(.1,.3,.5,.7,.9)
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
  if ( exists("dfun") ) { # likelihood may not be available for some models
    # add in model and contaminant likelihood
    dat$lm <- rep(NA,dim(dat)[1])
    # following not changed if no comtaminant model
    dat$lc <- dat$lm
    # add likelihoods to dat
    cells <- matrix(levels(dat$cell),ncol=2,byrow=T)
    row.names(cells) <- levels(dat$rcell)
    for (i in levels(D$D$rcell) ) { 
      isr <- dat$rcell==i & is.finite(dat[[D$RT]])
      isin <- D$D$rcell==i
      ts <- list(lower=numeric(0),upper=numeric(0))
      for (j in 1:2)
        if (any(dat$cell==cells[i,j])) ts[[j]] <- 
          dat[dat$cell==cells[i,j] & is.finite(dat[[D$RT]]),D$RT]
      L <- dfun(parlist=p[isin,],ts=ts)
      for (j in 1:2) {
        isin1 <- dat$cell==cells[i,j] & is.finite(dat[[D$RT]])
        dat[isin1,"lm"] <- L[[j]]
        dat[isin1,"lc"] <- p[isin,"pc"][1]*cp*p[isin1,"pg"][1] 
      }
    }
  }
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
    tmp <- rfun(ns,pi,tlohi,maxsamp=maxsamp,sim.mult=sim.mult)
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
      rfun(ns,pi,tlohi,maxsamp=sum(ns),onen=sum(ns),sim.mult=1,warn=FALSE)$sims
  }
  sim 
}

# Used by hydra.sfits
makedat <- function(sim,D,qp,ss,ssi) {
  # make empty data frame, no s column
  datj <- cbind(cell=numeric(0),rcell=numeric(0),
                s=character(0),D$D[0,-c(1:2)],D$D[0,D$R],
                numeric(0),numeric(0),numeric(0))   
  for ( rcell in levels(D$D$rcell) ) {
    Drc <- D$D[D$D$rcell==rcell,]
    isr1 <- sim[[rcell]]$r==Drc[[D$R]][1]
    nr <- length(isr1)
    nr1 <- sum(isr1)
    nr2 <- nr-nr1    
    if (nr1>0) { # reponse 1
      rt <- sim[[rcell]][isr1,"rt"]
      if ( !is.null(qp) ) {
        eqpij <- c(1:nr1)/(nr1+1)
        qpij <- qp[qp>=min(eqpij) & qp<=max(eqpij)] 
        nq <- length(qpij)
        rt <- cbind(c(quantile(rt,qpij,type=6),Inf),
                    qp=c(qpij,nr1/nr),qn=diff(c(0,qpij,1))*nr1)
      } else
        rt <- cbind(c(sort(rt),Inf),qp=c(c(1:nr1)/(nr1+1),nr1/nr),
                    qn=rep(nr1/(nr1+1),nr1+1))          
      R <- as.character(Drc[1,D$R])
      facs <- Drc[1,D$F]
      C <- Drc[1,D$SC[1]]
      datj <- rbind(datj,cbind(cell=rep(row.names(Drc)[1],dim(rt)[1]),
                               rcell=rep(rcell,dim(rt)[1]),s=rep(ssi,dim(rt)[1]),
                               facs[rep(1,dim(rt)[1]),],
                               rep(C,dim(rt)[1]),rep(R,dim(rt)[1]),rt)) 
    }                                
    if (nr2>0) { # reponse 2
      rt <- sim[[rcell]][!isr1,"rt"]
      if ( !is.null(qp) ) {
        eqpij <- c(1:nr2)/(nr2+1)
        qpij <- qp[qp>=min(eqpij) & qp<=max(eqpij)] 
        nq <- length(qpij)
        rt <- cbind(c(quantile(rt,qpij,type=6),Inf),
                    qp=c(qpij,nr2/nr),qn=diff(c(0,qpij,1))*nr2)
      } else
        rt <- cbind(c(sort(rt),Inf),qp=c(c(1:nr2)/(nr2+1),nr2/nr),
                    qn=rep(nr2/(nr2+1),nr2+1))  
      R <- as.character(Drc[2,D$R])
      facs <- Drc[2,D$F]
      C <- Drc[2,D$SC]
      datj <- rbind(datj,cbind(cell=rep(row.names(Drc)[2],dim(rt)[1]),
                               rcell=rep(rcell,dim(rt)[1]),s=rep(ssi,dim(rt)[1]),
                               facs[rep(1,dim(rt)[1]),],
                               rep(C,dim(rt)[1]),rep(R,dim(rt)[1]),rt))         
    }
  }
  names(datj) <- c("cell","rcell",D$S,D$F,D$SC,D$R,D$RT,"qp","qn")
  datj[[D$S]] <- factor(datj[[D$S]],levels=ss)
  datj[[D$R]] <- factor(datj[[D$R]],levels=levels(D$D[[D$R]]))
  datj[[D$SC[1]]] <- as.logical(as.character(datj[[D$SC[1]]]))
  datj[["cell"]] <- factor(datj[["cell"]],levels=1:dim(D$D)[1])
  datj  
}
