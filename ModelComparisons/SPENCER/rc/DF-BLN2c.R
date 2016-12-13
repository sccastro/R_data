# Start, distribution and random funcitons for ballistic lognromal drift 
# and distance to criterion 
# generic function names start=sfun PDF=dfun, CDF=pfun, random sample=rfun
# sfun called by mt.R, dfun and pfun by fitting.R, and rfun by simdat in latter

# tertune=.5; aprior=.5; pc=.01
sfun <-function(dati,pnams,pc=.01,aprior=.5,tertune=.5,st0=.05) {
# Generic LNR start routine for binary (correct/error) choice 
# Outputs standard m, s, r, ter parameterization by adding columns to
# design attribute of dati (data from one "rlevls" cell with both correct 
# and error data) assuming:
#  Pr(error)=(num_err+aprior)/(num_resp+2*aprior)
#  Equal correct and error variance, no correlation

  Di <- attr(dati,"D")
  minter <- as.numeric(dimnames(Di$nc)$minmax[1])*1.01
  # Define output
  ncols=3
  cnams <- c("m","v","ter")
  if (any(pnams=="st0")) {
    cnams <- c(cnams,"st0")
    ncols <- ncols+1
  }
  if (any(pnams=="pc")) {
    cnams <- c(cnams,"pc")
    ncols <- ncols+1
  }
  LNRstarts <- matrix(nrow=dim(Di$D)[1],ncol=ncols)
  dimnames(LNRstarts) <- list(NULL,cnams)    
  # identify respone levels
  rlevs <- as.numeric(levels(dati$rcell))
  for (i in rlevs) {
    # get data for current response level
    datii <- dati[dati$rcell==i,]
    # get cells that correspond to current response level
    tmp <- as.numeric(row.names(Di$D))[Di$D$rcell==i]
    Cname <- Di$SC[1] # ASSUMES FIRST SC SCORES CORRECT
    # identify cells for correct and error responses
    cc <- tmp[as.logical(Di$D[tmp,Cname])] 
    ec <- tmp[!as.logical(Di$D[tmp,Cname])]
    # robust proportion correct 
    Pc <- (sum(as.logical(datii[,Cname]))+aprior)/(dim(datii)[1]+2*aprior)
    rt <- datii[,Di$RT]
    ter <- pmax(tertune*min(rt[rt>minter]),minter)
    crt <- datii[as.logical(datii[,Cname]),Di$RT]-ter
    lcrt <- log(crt[crt>minter])
    m <- mean(lcrt)
    v <- var(lcrt)
    mue <- m-sqrt(2*v)*qnorm(1-Pc)
    nacc=length(cc)
    cc.cols <- cbind(rep(m,nacc),rep(v,nacc),rep(ter,nacc))
    nerr=length(ec)
    ec.cols <- cbind(rep(m,nerr),rep(v,nerr),rep(ter,nerr))  
    if (any(pnams=="st0")) {
      cc.cols <- cbind(cc.cols,rep(st0,nacc))
      ec.cols <- cbind(ec.cols,rep(st0,nerr))
    }
    if (any(pnams=="pc")) {
      cc.cols <- cbind(cc.cols,rep(pc,nacc))
      ec.cols <- cbind(ec.cols,rep(pc,nerr))
    }
    LNRstarts[cc,] <- cc.cols
    LNRstarts[ec,] <- ec.cols
  }
  LNRstarts <- data.frame(LNRstarts)
  LNRstarts$r <- 0 # same on natural and -1 to 1 fitting
  M2P(LNRstarts)
}


dfun=function(t,parlist,cv=NULL) {  
  # Defective PDF for responses on node 1 
  # st0>0: adds t0 variability convolution by numerical integration.
  # All parameters but st0 have length = number of accumulators. 
  # If length(st0)>1, first value used for all, other values ignored
  # NOT CHECKED FOR DIFFERENT TER FOR EACH NODE
  
  n1PDFfixedt0=function(t,mu,sig,r) {
    # Generates defective PDF for responses on node= 1 for 2 NODE ONLY
    # called by nPDF, assumes parameter vectors all same length
  
    # dont calculate for non-positive t
    ok=t[1,]>0
    t[1,!ok]=0
    t[1,ok] <- dlnorm(t[1,ok],mu[1],sig[1])
    # dont do F where f=0 or undefined
    bad=!is.finite(t[1,ok]) | t[1,ok]==0
    t[1,ok][bad]=0 
    ok[ok][bad]=F
    # Multiply by probability that node 2 not finished
    mu2=mu[2]+r*(log(t[2,ok])-mu[1])*sig[2]/sig[1]
    sig2=sqrt((1-r^2)*sig[2]^2)  
#     tmp <- plnorm(t[2,ok],mu2,sig2)
#     if (any(!is.finite(tmp))) {
#       cat("Oh CRAP/n")
#       print(r)
#       print(sig)
#       print(sig2)
#     }
#     t[1,ok] <- t[1,ok]*(1-tmp)    
    t[1,ok] <- t[1,ok]*(1-plnorm(t[2,ok],mu2,sig2))
    t[1,]
  }
  
  parlist$v <- sqrt(parlist$v)
  if (!any(names(parlist)=="st0")) { 
    outs=n1PDFfixedt0(mu=parlist$m,sig=parlist$v,r=parlist$r[1],
      t=matrix(rep(t,each=length(parlist$ter))-parlist$ter,
          nrow=length(parlist$ter)))
  } else {
    tmpf=function(t,termain,teroffset,mu,sig,r,st0) {
      ter=termain+teroffset
      t=matrix(rep(t,each=length(ter))-ter,nrow=length(ter))
      n1PDFfixedt0(t,mu,sig,r)/st0[1]
    }
    outs=numeric(length(t))
    teroffset=parlist$ter-parlist$ter[1]
    for (i in 1:length(outs)) 
      outs[i]=integrate(f=tmpf,lower=parlist$ter[1]-parlist$st0[1]/2,
        upper=ter[1]+st0[1]/2,t=t,teroffset=teroffset,
        mu=parlist$mu,sig=parlist$v,r=parlist$r,
        st0=parlist$st0[1])$value
  }
  outs
}

# t=dati[is.finite(dati[,D$RT]),D$RT]
# doCDF=F
# parlist=list(A=p$A[isin],b=p$b[isin],v=p$v[isin],
#   sv=p$sv[isin],ter=p$ter[isin],st0=p$st0[isin])
pfun=function(t,parlist,doCDF=T) {
  # Defective CDF for responses on node 1
  # doCDF=F returns inter-t probabilities
  # assumes t is sorted in increasing order with no ties

  outs <- numeric(length(t))
  t <- t-parlist$ter[1]
  parlist$ter=parlist$ter-parlist$ter[1]
  bad <- t<=0
  if (all(bad)) return(rep(0,length(outs)))
  bounds=c(0,t[!bad]) # Must be at least two as largest t is Inf
  badoffset=sum(bad)
  for (i in 1:(length(bounds)-1)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {
        outs[i+badoffset] <- 0
        break
      }       
      tmp=try(integrate(f=dfun,lower=bounds[i],upper=bounds[i+1],
        parlist=parlist)$value,silent=T)
      if (is.numeric(tmp)) {
        outs[i+badoffset] <- pmax(1e-10,tmp)
        break
      }
      # FAILED MIGHT BE GOOD TO ADD SMART BOUNDS
      return(rep(0,length(outs)))
    }
  }
  out <- outs/sum(outs)
  if (doCDF) cumsum(outs) else outs
}


# ns=nsim;pi=p[isin,]
# onen=1e4; maxsamp=1e6
rfun=function(ns,pi,tlohi=c(0,1),onen=1e4,maxsamp=1e6,warn=T,sim.mult=NA,remove.outliers=TRUE) {
# PARAMETERS
#   ns: number of RTs for each response
#   pi: data frame with columns ter, v sv, b and A (at least), one row per accumulator
#   tlohi: lower and upper limits of uniform contaminant
#   truncdrifts: each trial to have at least one positive drift?
#   onen: number to simulate per loop
#   maxsamp: limit on number of samples to take before giving up
# RETURN list with named elements 
#   rt: matrix[ns rows per choice,contaminate=0/1,choice=1:number of accumulators,rt]
#   n: number samples for each accumulator
  
  # take 5% more than biggest per cycle
  if (max(ns)<onen) onen <- round(max(ns)*1.05)
  nacc <- dim(pi)[1]
  RT <- matrix(ncol=sum(ns),nrow=3)
  Ns <- numeric(nacc)
  nsamp <- 0
  Ni <- 0
  full2col <- 0
  if (!is.na(sim.mult)) { # simulated subjects information
    n1 <- sum(ns)/sim.mult
    sim <- vector(mode="list",length=sim.mult)
    repn <- 0
  } else sim=NULL  
  repeat {
    nsamp <- nsamp + onen # number of samples taken
    rts <- matrix(rnorm(nacc*onen,pi$v,pi$sv),nrow=nacc) # actually v
    # remove all neg v
    if (truncdrifts) rts <- rts[,apply(rts,2,function(x){any(x>0)})] 
    rts <- (pi$b-runif(nacc*Ni,0,pi$A))/rts # divide by v
    rts[rts<0] <- Inf # negative drifts imply never complete
    rts <- rts+pi$ter # will allow different ter for each response
    # pick the smallest
    win <- apply(rts,2,function(x){m=which.min(x);c(m,x[m])}) 
    if (remove.outliers) {
      bad <- win[2,] <= tlohi[1] | win[2,] >= tlohi[2]
      win <- win[,!bad]
    }       
    Ni <- dim(win)[2]  
    # contaminate
    if (any(names(pi)=="pc"))
      win <- rbind(rbinom(Ni,1,pi$pc[1]),win) else
      win <- rbind(rep(0,Ni),win)
    if (any(as.logical(win[1,]))) { # contaminated
      win[2,as.logical(win[1,])] <- sample(nacc,sum(win[1,]),T)            # choice
      win[3,as.logical(win[1,])] <- runif(sum(win[1,]),tlohi[1],tlohi[2])  # rt
    }
    # number of each choice to add
    Nsi <- as.vector(table(factor(win[2,],1:nacc))) 
    nadd <- ns-Ns # number left to reach ns
    toadd <- !logical(Ni) 
    for (i in 1:nacc) { # set toadd to T for first nadd choices
      isi <- win[2,]==i
      if (any(isi)) toadd[isi] <- c(1:sum(isi))<=nadd[i]
    }
    endcol <- full2col+sum(toadd) # insert into RT
    if (endcol>full2col) # Is there something to add?
      RT[,(full2col+1):(endcol)] <- win[,toadd]
    full2col <- endcol
    Ns <- Ns + Nsi # number inserted
    if ( !is.null(sim) && (repn <= sim.mult) ) {
      sims <- data.frame(t(win))
      names(sims) <- c("contam","r","rt")
      sims$r <- factor(sims$r,levels=1:length(ns),labels=names(ns))
      navail <- pmin(dim(sims)[1] %/% n1,length(sim)-repn)
      if (navail>0) for (i in 1:navail) {
        sim[[repn+i]] <- sims[(((i-1)*n1)+1):(i*n1),]
      }
      repn <- repn+navail
    }  
    if ( ( (!is.null(sim) & (repn >= sim.mult)) & all(Ns >= ns) ) | nsamp >= maxsamp) break
  }
  if (warn && nsamp > maxsamp)
    warning(paste("Unable to get sufficient samples. Required =",
                  paste(ns,collapse=" "),"Obtained =",paste(Ns,collapse=" ")))
  names(Ns) <- names(ns)
  list(n=Ns,rt=cbind(contaminate=RT[1,],choice=RT[2,],rt=RT[3,]),sims=sim)
}



