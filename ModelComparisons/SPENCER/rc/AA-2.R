################################################################################
# ARCHITECTURE = 2 RACING ACCUMULATORS
#
# ONE ACCUMULATOR FOR *EACH* BINARY CORRECT/INCORRECT RESPONSE
# called by augment.design in mt.R to setup fast mapping of 
# parameters by fitp2modelp
get.reorder <- function(D,Cpars) {
  # D arguement is one subjects entry in fl$hmodelname$M
  # $D output reorders whole D matrix to put accumulators for response
  #    made (i.e., actually given) first in each lcell 
  # $P output used to map v parameter, assumes same (redundant) v for each entry 
  #    in lcell, correct v for correct lcell, error v for error lcell
  #    removes redundant and orders v to fit with D !! BEFORE $D reorder !!
  podr <- 1:dim(D$D)[1]
  dodr <- podr
  rcell <- D$D$rcell
  # !! assumes first score variable indicates correct !!
  C <- as.logical(as.character(D$D[[D$SC[1]]])) 
  R <- toupper(D$D[[D$R]])
  # !! assumes first latent is response !!
  lR <- toupper(D$D[[D$L[1]]]) 
  # matrix rcell x C position of first v
  pos <- tapply(podr,list(rcell,C),function(x){x[1]})
  # Correct and error response for rcell
  CR <- matrix(R[C],nrow=2)[1,][as.numeric(rcell)]
  # row within each lcell that should be put first
  isone <- C
  isone[C] <- lR[C]==CR[C]    # For correct lcells lR == CR
  isone[!C] <- lR[!C]!=CR[!C] # For error lcells lR != CR
  dn <- dimnames(pos)
  for (i in dn[[1]]) for (j in dn[[2]]) { # rcell and correct
    isin <- rcell==i & C==j 
    is1 <- isin & isone
    isnot1 <- isin & !isone 
    podr[is1] <- pos[i,j]
    podr[isnot1] <- pos[i,dn[[2]]!=j]          	
    dodr[isin] = c(dodr[is1],dodr[isnot1])
  }
  list(P=podr,D=dodr,N=Cpars)
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
  if (!any(names(p)=="pc")) 
    p[["pc"]] <- rep(1e-10,length(p[[1]])) else
    p$pc[p$pc<1e-10] <- 1e-10
  # Sort Cpars within lcell 
  if (!any(is.na(attr(D,"reorder")$N))) for (i in attr(D,"reorder")$N) 
    p[[i]] <- p[[i]][attr(D,"reorder")$P]
  # sorts in node 1 first format
  data.frame(p)[attr(D,"reorder")$D,] 
}

  
# p=strt
objective <- function(p,dat,pmap,pu,qmpe=F) {
  D <- attr(dat,"D")
  p <- fitp2modelp(p,pmap,D)
  nlcell <- length(levels(D$D$lcell))/length(levels(D$D$rcell))
  cp <- 1/diff(pu) # uniform contaminant density
  cp <- cp/nlcell  # unbiased guess
  if (!qmpe) { # Likelihood
    L <- vector(length=0)
    dat <- dat[is.finite(dat[,D$RT]),]
    for ( i in levels(D$D$lcell) ) {
      dati <- dat[dat$cell==i,]
      isin <- D$D$lcell==i
      if (dim(dati)[1]>0) {
        pc <- p$pc[isin][1] 
        L <- c(L, pc*cp + (1-pc)*
          dfun(t=dati[[D$RT]],parlist=p[isin,])) # likelihood
      }
    }
    L[!is.finite(L) | L<0] <- 1e-10
    mll <- -sum(log(L))
  } else { # qmpe        
    ns <- numeric(dim(dat)[2])
    qps <- ns
    cps <- ns
    for (i in levels(D$D$lcell)) { # get observed n (ns) and predicted probs (qps)
      isind <- dat$cell==i
      qs <- dat[isind,D$RT]
      ns[isind] <- dat[isind,"qn"]
      qps[isind] <- pfun(t=qs,parlist=p[D$D$lcell==i,],doCDF=F)
      cps[isind] <- diff(c(pu[1],qs[-length(qs)],pu[2]))*cp
    }
    # fix problems and contaminate
    bads <- qps<=1e-10
    for (i in levels(D$D$rcell)) { 
      isin <- dat$rcell==i
      if (any(bads[isin])) qps[isin]/sum(qps[isin])
      pc <- p$pc[D$D$lcell==i][1]
      qps[isin] <- pc*cps[isin] + (1-pc)*qps[isin]
    }
    mll <- -sum(ns*log(qps))

  }
  mll # minus log-likelihood
}


# maxsamp=1e6
simdat <- function(best,pmap,dat,sim.mult=100,maxsamp=1e6) {
  
  qs <- function(x){c(1:length(x))/(length(x)+1)}
  
  D <- attr(dat,"D")
  cnam <- D$SC[1]
  rtnam <- D$RT
  rcnams <- levels(D$D$rcell)
  nc <- D$nc
  tlohi=as.numeric(dimnames(nc)$minmax)
  nlcell <- length(levels(D$D$lcell))/length(levels(D$D$rcell))
  cp <- 1/(diff(tlohi)*nlcell) # equiprobable guess 
  p <- fitp2modelp(best,pmap,D)
  if (exists("dfun")) { # likelihood may not be available for some models
    # add in model and contaminant likelihood
    dat$lm <- rep(NA,dim(dat)[1])
    # following not changed if no comtaminant model
    dat$lc <- dat$lm
    # add likelihoods to dat
    for (i in levels(D$D$lcell)) {
      isrt <- dat$cell==i & is.finite(dat[,D$RT])
      isin <- D$D$lcell==i
      if (any(isrt)) {
        dat$lm[isrt] <- dfun(t=dat[isrt,D$RT],parlist=p[isin,])
        if (any(names(p)=="pc")) {           # contaminant model
          pc <- p$pc[isin][1]/nlcell
          dat$lc[isrt] <- pc*cp 
          dat$lm [isrt] <- dat$lm[isrt]*(1-pc)
        }  
      }
    }
  }
  # get prdicted RT via simulation
  D <- D$D[attr(D,"reorder")$D,]
  srtnam <- paste(rtnam,"sim",sep=".")
  dat[[srtnam]] <- dat[[rtnam]]
  sqpnam <- "qp.sim"
  dat[[sqpnam]] <- dat$qp
  isrt <- is.finite(dat[[rtnam]])
  #  attr(dat,"D") <- NULL
  for (i in 1:length(rcnams)) {
    # correct so first v biggest
    isin <- as.character(D$rcell) == rcnams[i] & 
      as.logical(as.character(D[[cnam]]))
    # simuluate sim.mult many times as observed
    crct <- dat[dat$rcell==rcnams[i] & isrt,cnam]
    
    # specific to binary true/false case
    nsim <- sim.mult*tapply(dat$qn[dat$rcell==rcnams[i]],
                            dat[[cnam]][dat$rcell==rcnams[i]],sum)    
    if (length(nsim)==1) if (names(nsim)=="TRUE") 
    {nsim <- c(0,nsim); names(nsim)[1]="FALSE"} else
    {nsim <- c(nsim,0); names(nsim)[2]="TRUE"}    
    nsim[is.na(nsim)] <- 0
    nsim <- nsim[2:1] # Correct first    
    tmp <- rfun(nsim,p[isin,],tlohi,maxsamp=maxsamp)
    #     
    #     print(cbind(D[D$rcell==i & as.logical(D$C),],p[isin,]))
    #     print(c(tmp$n,tmp$n[1]/sum(tmp$n)))
    #     
    srt <- data.frame(tmp$rt[!is.na(tmp$rt[,1]),])    
    srt$choice <- factor(srt$choice,levels=1:2,labels=names(tmp$n))
    for (j in names(nsim)) {
      isin1 <- as.character(dat$rcell)==rcnams[i] & dat[[cnam]]==j 
      isin2 <- isin1 & isrt
      if (sum(isin2)==0 | dim(srt)[1]==0) {
        dat[isin1,srtnam] <- rep(NA,sum(isin1)) 
        dat[isin1,sqpnam] <- rep(NA,sum(isin1))
      } else {
        dat[isin2,srtnam] <- quantile(srt[as.logical(srt$choice)==j,"rt"],
                                      type=1,probs=dat[isin2,"qp"])
        dat[isin2,sqpnam] <- table(cut(srt[as.logical(srt$choice)==j,"rt"],
                                       c(-Inf,dat[isin1,rtnam]) ))[-(sum(isin1))]/
                                         length(srt[srt$choice==j,"rt"])
        dat[isin1 &!isin2,sqpnam] <- tmp$n[as.character(j)]/sum(tmp$n)
      }
    }
  }
  dat[,c("lm","lc",srtnam,sqpnam)] 
}


