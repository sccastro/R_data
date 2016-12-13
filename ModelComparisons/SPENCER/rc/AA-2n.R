################################################################################
# ARCHITECTURE = n RACING ACCUMULATORS
#
# ONE ACCUMULATOR FOR *EACH* RESPONSE
# called by augment.design in mt.R to setup fast mapping of 
# parameters by fitp2modelp
get.reorder <- function(D,Cpars,ok,SR=NULL) {
  # D arguement is one subjects entry in fl$hmodelname$M
  # $D output reorders whole D matrix to put accumulators for response
  #    made (i.e., actually given) first in each lcell 
  # $P output used to map v parameter, assumes same (redundant) v for each entry 
  #    in lcell, correct v for correct lcell, error v for each error lcell
  #    removes redundant and orders v to fit with D !! BEFORE $D reorder !!
  podr <- 1:dim(D$D)[1]
  dodr <- podr
  rcell <- D$D$rcell
  lcell <- D$D$lcell
  # !! assumes first score variable indicates correct !!
  C <- as.logical(as.character(D$D[[D$SC[1]]])) 
  R <- toupper(D$D[[D$R]])
  # !! assumes first latent is response !!
  lR <- toupper(D$D[[D$L[1]]]) 
  # matrix rcell x C position of first v
  pos <- tapply(podr,list(rcell,C),function(x){x[1]})
  CR <- rep(R[C],each=length(levels(D$D[[D$R]])))[ok]
  CR <- CR == lR # correct response before dodr sort
  isone <- lR==R # latent = observed response
  for ( i in levels(rcell) ) {
    podr[rcell==i & CR] <- pos[i,"TRUE"]
    podr[rcell==i & !CR] <- pos[i,"FALSE"]
  }
  for ( i in levels(lcell) )
    dodr[lcell==i] <- c(dodr[lcell==i & isone],dodr[lcell==i & !isone])
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
objective <- function(p,dat,pmap,pu,qmpe=FALSE,precision=2.5) {
  D <- attr(dat,"D")  
# badp <<- p
  p <- fitp2modelp(p,pmap,D)  
  nlcell <- length(levels(D$D$lcell))/length(levels(D$D$rcell))
  cp <- 1/diff(pu) # uniform contaminant density
  cp <- cp/nlcell  # unbiased guess
  if ( !qmpe ) { # Likelihood
    L <- vector(length=0)
    dat <- dat[is.finite(dat[,D$RT]),]
    for ( i in levels(D$D$lcell) ) {    
#  print(i)     
      dati <- dat[dat$cell==i,]
      isin <- D$D$lcell==i
      if ( dim(dati)[1]>0 ) {
        pc <- p$pc[isin][1]         
        L <- c(L, pc*cp + (1-pc)*
          dfun(t=dati[[D$RT]],parlist=p[isin,],cv=dati[,D$cv,drop=F]),
          precision=precision) # likelihood
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
      qps[isind] <- pfun(t=qs,parlist=p[D$D$lcell==i,],doCDF=F,
        precision=precision)
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


# maxsamp=1e6; qps <- c(.1,.3,.5,.7,.9)
simdat <- function(best,pmap,dat,sim.mult=100,maxsamp=1e6,qps=c(.1,.3,.5,.7,.9)) {
  
  qs <- function(x){c(1:length(x))/(length(x)+1)}
  
  D <- attr(dat,"D")
  cnam <- D$SC[1]
  rtnam <- D$RT
  rnam <- D$R
  lrnam <- D$L[1]
  rcnams <- levels(D$D$rcell)
  nc <- D$nc
  tlohi=as.numeric(dimnames(nc)$minmax)
  nlcell <- length(levels(D$D$lcell))/length(levels(D$D$rcell))
  cp <- 1/(diff(tlohi)*nlcell) # equiprobable guess 
  p <- fitp2modelp(best,pmap,D)
  if ( exists("dfun") ) { # likelihood may not be available for some models
    # add in model and contaminant likelihood
    dat$lm <- rep(NA,dim(dat)[1])
    # following not changed if no comtaminant model
    dat$lc <- dat$lm
    # add likelihoods to dat
    for (i in levels(D$D$lcell)) {
      isrt <- dat$cell==i & is.finite(dat[,D$RT])
      isin <- D$D$lcell==i
      if ( any(isrt) ) {
        dat$lm[isrt] <- dfun(t=dat[isrt,D$RT],parlist=p[isin,],cv=dat[isrt,D$cv,drop=F])
        if (any(names(p)=="pc")) {           # contaminant model
          pc <- p$pc[isin][1]/nlcell
          dat$lc[isrt] <- pc*cp 
          dat$lm [isrt] <- dat$lm[isrt]*(1-pc)
        }  
      }
    }
  }
  if (is.null(D$cv)) cv="" else cv <- D$cv
  RT <- D$RT
  R <- D$R
  # get prdicted RT via simulation
  D <- D$D[attr(D,"reorder")$D,]
  srtnam <- paste(rtnam,"sim",sep=".")
  dat[[srtnam]] <- dat[[rtnam]]
  sqpnam <- "qp.sim"
  dat[[sqpnam]] <- dat$qp
  if ( cv!="" ) {
    Csimnam <- "C.sim"
    sqsnam <- paste("q",qps,sep="-")
    for (i in sqsnam) dat[[i]] <- dat$qp
  } else {
    Csimnam <- NULL
    sqsnam <- NULL
  }
  isrt <- is.finite(dat[[rtnam]])
  #  attr(dat,"D") <- NULL
  sim <- vector(mode="list",length=length(rcnams))
  names(sim) <- rcnams
  for (i in 1:length(rcnams)) {
    
#               print(i)
    
    # correct so first v biggest
    isin <- as.character(D$rcell) == rcnams[i] & 
      as.logical(as.character(D[[cnam]]))
    # simuluate sim.mult many times as observed
#     crct <- dat[dat$rcell==rcnams[i] & isrt,cnam]
    nsim <- sim.mult*tapply(dat$qn[dat$rcell==rcnams[i]],
      dat[[rnam]][dat$rcell==rcnams[i]],sum)    
    nsim[is.na(nsim)] <- 0
    names(nsim) <- tolower(names(nsim))
    tmp=D[D$rcell==rcnams[i],]
    nsim <- nsim[tolower(as.character(tmp[as.logical(tmp[[cnam]]),lrnam]))]
    if ( cv!="" ) {
      dati <- dat[dat$rcell==rcnams[i] & is.finite(dat[[RT]]),]
      pi=parlistt(p[isin,],dati[,cv,drop=F])
      for ( j in 1:dim(dati)[1] ) {
        
#                 print(j)
        
        pit <- data.frame(lapply(pi,function(x){x[j,]}))
        ns  <- nsim
        ns[names(ns)==as.character(dati[j,R])] <- sim.mult
        ns[names(ns)!=as.character(dati[j,R])] <- 0        
        tmpi <- rfun(ns=ns,pi=pit,tlohi=tlohi,maxsamp=maxsamp,sim.mult=sim.mult)
        if (j==1) tmp <- tmpi else {
          tmp$n <- rbind(tmp$n,tmpi$n)
          tmp$rt <- rbind(tmp$rt,tmpi$rt)
          tmp$sims <- c(tmp$sims,tmpi$sims)
        }
      }
      tmp$n <- data.frame(tmp$n)
    } else tmp <- rfun(ns=nsim,pi=p[isin,],tlohi=tlohi,maxsamp=maxsamp,sim.mult=sim.mult)
    sim[[i]] <- tmp$sims
    srt <- data.frame(tmp$rt[!is.na(tmp$rt[,1]),])    
    srt$choice <- factor(srt$choice,levels=1:length(names(tmp$n)),labels=names(tmp$n))
    for ( j in names(nsim) ) {
      isin1 <- as.character(dat$rcell)==rcnams[i] & tolower(as.character(dat[[rnam]]))==j 
      isin2 <- isin1 & isrt
      if ( sum(isin2)==0 | dim(srt)[1]==0 ) {
        dat[isin1,srtnam] <- rep(NA,sum(isin1)) 
        dat[isin1,sqpnam] <- rep(NA,sum(isin1))
      } else {
        dat[isin2,srtnam] <- quantile(srt[as.character(srt$choice)==j,"rt"],
            type=1,probs=dat[isin2,"qp"])
        dat[isin2,sqpnam] <- table(cut(srt[as.character(srt$choice)==j,"rt"],
          c(-Inf,dat[isin1,rtnam])) )[-(sum(isin1))]/
          length(srt[as.character(srt$choice)==j,"rt"])
        if ( cv!="" ) 
          dat[isin1 &!isin2,sqpnam] <- sum(tmp$n[,as.character(j)])/sum(tmp$n) else 
          dat[isin1 &!isin2,sqpnam] <- tmp$n[as.character(j)]/sum(tmp$n)
      }
    }
    if ( cv!="" ) {
      isin3 <- dat$rcell==rcnams[i] & is.finite(dat[[RT]])
      temp <- table(dat[isin3,c(R,cnam)])
      # get name of correct response
      if (any(dimnames(temp)[[2]]=="TRUE"))
        crnam <- dimnames(temp)[[1]][temp[,"TRUE"]!=0] else
        crnam <- dimnames(temp)[[1]][temp[,"FALSE"]==0]
      dat[isin3,Csimnam] <- tmp$n[,crnam]/apply(tmp$n,1,sum)
      dat[isin3,Csimnam][is.nan(dat[isin3,Csimnam])] <- 0
      for (k in 1:sum(isin3)) {
        indx <- (sim.mult*(k-1)+1):(sim.mult*k)
        tmprt <- tmp$rt[indx,3]
        tmprt <- tmprt[!is.na(tmprt)]
        if ( length(qps)>length(tmprt) )
          warning("More quantiles than simulated data points!")
        dat[isin3,sqsnam][k,] <- quantile(tmprt,qps)
      }
      dat[isin1 &!isin2,sqpnam] <- sum(tmp$n[,as.character(j)])/sum(tmp$n) 
    }
  }
  list(dat=dat[,c("lm","lc",srtnam,sqpnam,Csimnam,sqsnam)],sim=sim)
}




sim1dat <- function(best,pmap,dat) {
  
  D <- attr(dat,"D")
  cnam <- D$SC[1]
  rtnam <- D$RT
  rnam <- D$R
  lrnam <- D$L[1]
  rcnams <- levels(D$D$rcell)
  nc <- D$nc
  tlohi=as.numeric(dimnames(nc)$minmax)
  nlcell <- length(levels(D$D$lcell))/length(levels(D$D$rcell))
  cp <- 1/(diff(tlohi)*nlcell) # equiprobable guess 
  p <- fitp2modelp(best,pmap,D)
  D <- D$D[attr(D,"reorder")$D,]
  srtnam <- paste(rtnam,"sim",sep=".")
  dat[[srtnam]] <- dat[[rtnam]]
  sqpnam <- "qp.sim"
  dat[[sqpnam]] <- dat$qp
  isrt <- is.finite(dat[[rtnam]])
  #  attr(dat,"D") <- NULL
  sim <- vector(mode="list",length=length(rcnams))
  names(sim) <- rcnams
  for (i in 1:length(rcnams) ) {
    isin <- as.character(D$rcell) == rcnams[i] & 
      as.logical(as.character(D[[cnam]]))
    nsim <- tapply(dat$qn[dat$rcell==rcnams[i]],
                            dat[[rnam]][dat$rcell==rcnams[i]],sum)    
    nsim[is.na(nsim)] <- 0
    names(nsim) <- tolower(names(nsim))
    tmp=D[D$rcell==rcnams[i],]
    nsim <- nsim[tolower(as.character(tmp[as.logical(tmp[[cnam]]),lrnam]))]
    sim[[i]] <- 
      rfun(nsim,p[isin,],tlohi,maxsamp=sum(nsim),onen=sum(nsim),sim.mult=1,warn=FALSE)$sims[[1]]
  }
  sim 
}


# Used by hydra.sfits
makedat <- function(sim,D,qp,ss,ssi) {
  # make empty data frame, no s column
  datj <- cbind(cell=numeric(0),rcell=numeric(0),
                s=character(0),D$D[0,-c(1:2)],D$D[0,D$R],
                numeric(0),numeric(0),numeric(0))   
  rlevs <- levels(D$D[,D$R])
  for ( rcell in levels(D$D$rcell) ) {
    Drc <- D$D[D$D$rcell==rcell,]
    for (j in rlevs) {
      isr1 <- sim[[rcell]]$r==j
      nr1 <- sum(isr1)
      nr <- length(isr1)
      rt <- sim[[rcell]][isr1,"rt"]
      if ( any(isr1) ) {
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
        C <- Drc[Drc[[D$R]]==j,D$SC[1]][1]
        datj <- rbind(datj,cbind(cell=rep(Drc[Drc[[D$R]]==j,]$lcell[1],dim(rt)[1]),
                                 rcell=rep(rcell,dim(rt)[1]),s=rep(ssi,dim(rt)[1]),
                                 facs[rep(1,dim(rt)[1]),],
                                 rep(C,dim(rt)[1]),rep(j,dim(rt)[1]),rt)) 
      }
    }
  }
  names(datj) <- c("cell","rcell",D$S,D$F,D$SC[1],D$R,D$RT,"qp","qn")
  datj[[D$S]] <- factor(datj[[D$S]],levels=ss)
  datj[[D$R]] <- factor(datj[[D$R]],levels=levels(D$D[[D$R]]))
  datj[[D$SC[1]]] <- as.logical(as.character(datj[[D$SC[1]]]))
  datj  
}


