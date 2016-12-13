################################################################################
# ARCHITECTURE = 1 Binary ACCUMULATOR 
# TASK SWITCHING VERSION, ASSUMES SEVERAL NAMED VARIABLE PRESENT 
# NUMERIC TRIAL COVARIATE (cv) CALLED CI = CUE-TARGET INTERVAL UNIQUE IN EACH rcell
# t0 & tc: LOWER BOUND OF NON-DECISION AND CUE PROCESSING TIME
# stc: WIDTH OF UNIFORM CUE PROCESSING TIME (CAN BE ABSENT)
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
objective <- function(p,dat,pmap,pu,qmpe=F) {
  D <- attr(dat,"D")
#   cat("\n")
#   print(p["sv.I"])   
  p <- fitp2modelp(p,pmap,D)
#   print(p$sv[1])
  cp <- 1/diff(pu) # contamination density
  if (!qmpe) { # Likelihood
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
      
      ci <- dat[isr,"CI"][1]
      P <- p[D$D$rcell==i,]
      if (!any(names(P)=="stc")) {
        P$t0 <- P$t0 + pmax(0,P$tc-ci)
        qps[isr] <- unlist(pfun(parlist=P,ts=ts,doCDF=F),use.names=F)        
      } else {
        CI <- ci - P$tc[1]
        if (CI<=0) { # Target before lower bound of cue processing time
          P$t0 <- P$t0 - CI
          P$st <- P$stc
          qps[isr] <- unlist(pfun(parlist=P,ts=ts,doCDF=F),use.names=F)                
        } else {
          PFP <- P # parameters for fully processed cue
          PFP$t0 <- 1e-6
          PFP$st <- 1e-6
          if (CI >= P$stc[1]) { # cue fully processed
            qps[isr] <- unlist(pfun(parlist=PFP,ts=ts,doCDF=F),use.names=F)                
          } else {
            pcfp <- CI/p$stc[1] # probability cue fully processed
            P$t0 <- 1e-6
            P$st <- P$stc-CI
            qps[isr] <- pcfp*unlist(pfun(parlist=PFP,ts=ts,doCDF=F),use.names=F) +               
              (1-pcfp)*unlist(pfun(parlist=P,ts=ts,doCDF=F),use.names=F)                
          }
        }
      }
      # contaminated quantile probabilities
      if (length(ts$lower)!=0) 
        cpl <- diff(c(pu[1],ts$lower,pu[2]))*p[D$D$rcell==i,"pg"] else 
        cpl <- numeric(0)
      if (length(ts$upper)!=0) 
        cpu <- diff(c(pu[1],ts$upper,pu[2]))*(1-p[D$D$rcell==i,"pg"]) else 
        cpu <- numeric(0)
      qps[isr] <- 
        p$pc[D$D$rcell==i]*c(cpl,cpu)*cp + (1-p$pc[D$D$rcell==i])*qps[isr]      
    }
    mll <- -sum(dat$qn*log(qps))
  }
  mll # minus log-likelihood
}



simdat <- function(best,pmap,dat,sim.mult=100,maxsamp=1e6) {
  
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
#         dat$lm[isrt] <- dfun(t=dat[isrt,D$RT],parlist=p[i,])
#         if (any(names(p)=="pc")) { # contaminant model
#           dat$lc[isrt] <- p[i,"pc"]*cp*p[i,"pg"]
#           dat$lm[isrt] <- dat$lm[isrt]*(1-p[i,"pc"])
#         }  
#       }
#     }
#   }
  dat$lm <- rep(NA,dim(dat)[1])
  dat$lc <- dat$lm
  # get prdicted RT via simulation
  srtnam <- paste(rtnam,"sim",sep=".")
  dat[[srtnam]] <- dat[[rtnam]]
  sqpnam <- "qp.sim"
  dat[[sqpnam]] <- dat$qp
  isrt <- is.finite(dat[[rtnam]])
  for (i in levels(D$D$rcell) ) {
    Di <- D$D[D$D$rcell==i,]
    pi <- p[D$D$rcell==i,][
      c(1:2)[Di[[D$R]]==levels(Di[[D$R]])[1]],,drop=F]
    dati <- dat[dat$rcell==i,]
    dati$cell <- factor(as.character(dati$cell),row.names(Di))
    nsim <- sim.mult*tapply(dati$qn,dati$cell,sum)
    nsim[is.na(nsim)] <- 0
    ci <- dati[,"CI"][1]
    P <- pi
    if (!any(names(P)=="stc")) {
      P$t0 <- P$t0 + pmax(0,P$tc-ci)
      P$st <- 1e-6
      tmp <- rfun(nsim,P,tlohi,maxsamp=maxsamp)        
    } else {
      CI <- ci - P$tc[1]
      if (CI<=0) { # Target before lower bound of cue processing time
        P$t0 <- P$t0 - CI
        P$st <- P$stc
        tmp <- rfun(nsim,P,tlohi,maxsamp=maxsamp)        
      } else {
        PFP <- P # parameters for fully processed cue
        PFP$t0 <- 1e-6
        PFP$st <- 1e-6
        if (CI >= P$stc[1]) { # cue fully processed
          tmp <- rfun(nsim,PFP,tlohi,maxsamp=maxsamp)        
        } else {
          pcfp <- CI/p$stc[1] # probability cue fully processed
          P$t0 <- 1e-6
          P$st <- P$stc-CI
          tmp <- rfun(ceiling(nsim*pcfp),PFP,tlohi,maxsamp=maxsamp)
          tmp1 <- rfun(ceiling(nsim*(1-pcfp)),P,tlohi,maxsamp=maxsamp)
          tmp$n <- tmp$n + tmp1$n
          tmp$rt <- rbind(tmp$rt,tmp1$rt)
        }
      }
    }
    srt <- data.frame(tmp$rt[!is.na(tmp$rt[,1]),])    
    srt$choice <- factor(srt$choice,levels=0:1,labels=names(tmp$n))
    for (j in names(nsim)) {
      isin1 <- dat$cell==j 
      isin2 <- isin1 & is.finite(dat[[rtnam]])
      if ( sum(isin1)==0 | dim(srt)[1]==0 ) {
        dat[isin1,srtnam] <- rep(NA,sum(isin1)) 
        dat[isin1,sqpnam] <- rep(NA,sum(isin1))
      } else {
        dat[isin2,srtnam] <- quantile(srt[srt$choice==j,"rt"],
          type=1,probs=dat[isin2,"qp"])
        dat[isin2,sqpnam] <- table(cut(srt[srt$choice==j,"rt"],
          c(-Inf,dat[isin1,rtnam]) ))[-(sum(isin1))]/
          length(srt[srt$choice==j,"rt"])
        dat[isin1 &!isin2,sqpnam] <- tmp$n[as.character(j)]/sum(tmp$n)
      }
    }
  }
  dat[,c("lm","lc",srtnam,sqpnam)] 
}

