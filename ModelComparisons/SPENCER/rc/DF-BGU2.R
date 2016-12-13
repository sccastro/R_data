
# Distribution random and start funcitons for ballistic gamma drift uniform 
# start point two choice model
#
# Start, distribution and random funcitons for Uniform/Gamma, 
# with k=shape and s=scale parameterization
# mean = ks, variance =ks^2, skew=2/k^.5

# generic function names start=sfun PDF=dfun, CDF=pfun, random sample=rfun
# sfun called by mt.R, dfun and pfun by fitting.R, and rfun by simdat in latter

#aprior=.5;minter=.1;tertune=.9
sfun <-function(dati,pnams,pc=.01,K=2,
  aprior=.5,tertune=.9,Atune=.5,st0=.05) {
# Generic UG start routine for binary (correct/error) choice 
# Outputs standard b, A, k, s, ter parameterization by adding columns to
# design attribute of dati (data from one "rlevls" cell with both correct 
# and error data) assuming:
#  Pr(error)=(num_err+aprior)/(num_resp+2*aprior)
#  K is defualt shape parameter
#  Use correct RT - ter to solve scale for both accumulators 
#  s = var/K
#  pow <- c(2,3,4)[c(Pc<.9,Pc>=.9 & Pc<=.95,Pc>.95)]
#  To get error scale sE= sC/(2*Pc)^pow, crudely gets smaller with bigger proportion correct 
#  ter = tertune * min(rt)
#  b = 2 * (mean_correct_rt-ter)/ks  (i.e., mean "rate")
#  A = b*Atune
  
  Di <- attr(dati,"D")
  minter <- as.numeric(dimnames(Di$nc)$minmax[1])*1.01
  # Define output
  ncols=5
  cnams <- c("b","A","k","s","ter")
  if (any(pnams=="st0")) {
    cnams <- c(cnams,"st0")
    ncols <- ncols+1
  }
  if (any(pnams=="pc")) {
    cnams <- c(cnams,"pc")
    ncols <- ncols+1
  }
  LBAstarts <- matrix(nrow=dim(Di$D)[1],ncol=ncols)
  dimnames(LBAstarts) <- list(NULL,cnams)    
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
    # relevant RTs
    rt <- datii[,Di$RT]
    # heuristics for v bounded 0-1
    kC <- K; kE <- K
    if (length(rt[as.logical(datii[[Cname]])])<2) 
      crt <- rt else crt <- rt[as.logical(datii[[Cname]])]
    if (length(crt)<2)
      stop("Less than 2 RTs in cell")
    sC <- (IQR(crt)/1.2) # robust version of var(crt), for normal 1.349, adjust for =2    
    pow <- c(2,3,4)[c(Pc<.9,Pc>=.9 & Pc<=.95,Pc>.95)]
    sE <- sC/(2*Pc)^pow # gets smaller with proportion correct
    ter <- pmax(tertune*min(rt[rt>minter])-st0,minter)
    b <- ((mean(crt)-ter)*kC*sC)/(1-Atune)
    A <- b*Atune
    nacc=length(cc)
    cc.cols <- cbind(rep(b,nacc),rep(A,nacc),rep(kC,nacc),
      rep(sC,nacc),rep(ter,nacc))
    nerr=length(ec)
    ec.cols <- cbind(rep(b,nerr),rep(A,nerr),rep(kE,nerr),
      rep(sE,nerr),rep(ter,nerr))  
    if (any(pnams=="st0")) {
      cc.cols <- cbind(cc.cols,rep(st0,nacc))
      ec.cols <- cbind(ec.cols,rep(st0,nerr))
    }
    if (any(pnams=="pc")) {
      cc.cols <- cbind(cc.cols,rep(pc,nacc))
      ec.cols <- cbind(ec.cols,rep(pc,nerr))
    }
    LBAstarts[cc,] <- cc.cols
    LBAstarts[ec,] <- ec.cols
  }
  M2P(data.frame(LBAstarts))
}


# n1mean=function(x0max,chi,shape,scale) {
#   # Generates mean RT for responses on node #1. 
#   pc=n1CDF(Inf,x0max,chi,shape,scale)
#   fn=function(t,x0max,chi,shape,scale,st0=0,pc)
#     t*n1PDF(t,x0max,chi,shape,scale,st0)/pc
#   tmp=integrate(f=fn,lower=0,upper=100*chi,x0max=x0max,chi=chi,pc=pc,
#                 shape=shape,scale=scale,st0=st0)$value
#   list(mean=tmp,p=pc)
# }


dfun=function(t,parlist,silent=TRUE,cv=NULL) {  
  # Defective PDF for responses on node 1 
  # st0>0: adds t0 variability convolution by numerical integration.
  # If length ter or st0 >1, first value used for all, other values ignored

  ################# SINGLE UNIT FUNCTIONS

  # Density for Inverse Gamma function. Copied from MCMCpack
  # library. Saves loading the whole thing.
  dinvgamma = function (x, shape, scale = 1) {
    if (shape <= 0 | scale <= 0)
      stop("Shape or scale parameter negative in dinvgamma().\n")
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - 
      (alpha + 1) * log(x) - (beta/x)
    return(exp(log.density))
  }
  
  pinvgamma = function(x, shape, scale =1 )  1-Rgamma(shape,scale/x)
  
  # C/I/Rgamma functions copied from library ZipfR, saves loading whole thing.
  Rgamma = function (a, x, lower = TRUE, log = !missing(base), base = exp(1)) {
    if (log) {
      if (missing(base)) {
        pgamma(x, shape = a, scale = 1, lower.tail = lower, log = TRUE)
      } else {
        pgamma(x, shape = a, scale = 1, lower.tail = lower, log = TRUE)/log(base)
      }
    } else {
      pgamma(x, shape = a, scale = 1, lower.tail = lower, log = FALSE)
    }
  }
  
  Igamma=function (a, x, lower = TRUE, log = !missing(base), base = exp(1)) {
    if (log) {
      Cgamma(a, log = TRUE, base = base) + 
        Rgamma(a, x, lower = lower, log = TRUE, base = base)
    } else {
      Cgamma(a, log = FALSE) * 
        Rgamma(a, x, lower = lower, log = FALSE)
    }
  }
  
  Cgamma=function (a, log = !missing(base), base = exp(1)) {
    if (log) {
      if (missing(base)) {
        lgamma(a)
      } else {
        lgamma(a)/log(base)
      }
    } else {
      gamma(a)
    }
  } 
  
  # Original stuff below here.
  
  fptpdf=function(z,x0max,chi,phi1,phi2) {
    # E.J.'s original parameterisation in terms of shape
    # and scale (phi1 and phi2 respectively).
    if (x0max<1e-10) return(dinvgamma(z, phi1, chi/phi2))
    tmp1=z*phi2
    (phi1*phi2)/x0max * ( Rgamma(phi1+1, chi/tmp1) -
      Rgamma(phi1+1, (chi-x0max)/tmp1) )
  }
  
  fptcdf=function(z,x0max,chi,phi1,phi2) {
    # E.J.'s original parameterisation in terms of shape
    # and scale (phi1 and phi2 respectively).
    if (x0max<1e-10) return(pinvgamma(z, phi1, chi/phi2))
    tmp1=z*phi2 ; tmp2=chi/tmp1 ; tmp3=chi-x0max ; tmp4=tmp3/tmp1 ; tmp5=tmp1/x0max
    x1 = Igamma(phi1, tmp2, lower=FALSE) 
    x2 = Igamma(phi1+1, tmp2, lower=TRUE)
    x3 = Igamma(phi1, tmp4, lower=FALSE)
    x4 = Igamma(phi1+1, tmp4, lower=TRUE)
    (1/(gamma(phi1))) * ( x1*(chi/x0max) + x2*tmp5 - x3*(tmp3/x0max) - x4*tmp5 )   
  }
  
  
  n1PDFfixedt0=function(t,A,b,k,s) {
    # Generates defective PDF for responses on node #1.
    N <- length(k) # Number of responses.
    if (N>2) {
      tmp=array(dim=c(length(t),N-1))
      for (i in 2:N) tmp[,i-1] <- 
        fptcdf(z=t,x0max=A[i],chi=b[i],phi1=k[i],phi2=s[i])
      G=apply(1-tmp,1,prod)
    } else {
      G=1-fptcdf(z=t,x0max=A[2],chi=b[2],phi1=k[2],phi2=s[2])
    }
    out <- G*fptpdf(z=t,x0max=A[1],chi=b[1],phi1=k[1],phi2=s[1])
    out[t<=0 | !is.finite(out)] <- 0
    out
  }
  
    
  t <- t-parlist$ter[1]
  if (!any(names(parlist)=="st0") || parlist$st0[1]<.001) { 
    return(n1PDFfixedt0(t,
      A=parlist$A,b=parlist$b,k=parlist$k,s=parlist$s))
  } else {
    tmpf <- function(t,A,b,k,s,st0) 
      n1PDFfixedt0(t,A,b,k,s)/st0
    outs <- numeric(length(t))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=tmpf,
        lower=t[i]-parlist$st0[1],upper=t[i],
        A=parlist$A,b=parlist$b,k=parlist$k,s=parlist$s,
        st0=parlist$st0[1])$value,silent=silent)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}

# t=dati[is.finite(dati[,D$RT]),D$RT]
# doCDF=F
# parlist=list(A=p$A[isin],b=p$b[isin],v=p$v[isin],
#   sv=p$sv[isin],ter=p$ter[isin],st0=p$st0[isin])
pfun=function(t,parlist,doCDF=T,silent=TRUE) {
  # Defective CDF for responses on node 1
  # doCDF=F returns inter-t probabilities
  # assumes t is sorted in increasing order with no ties

  outs <- numeric(length(t))
  bad <- t-parlist$ter[1]<=0
  if (all(bad)) return(outs)
  bounds <- c(0,t[!bad]) # Must be at least two as largest t is Inf
  badoffset <- sum(bad)
  if (any(names(parlist)=="st0") && parlist$st0[1]<1e-6)
    parlist$st0[1] <- 0 # integral of intergal can fail for st0 small
  for (i in 1:(length(bounds)-1)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {
        outs[i+badoffset] <- 0
        break
      }       
      tmp=try(integrate(f=dfun,lower=bounds[i],upper=bounds[i+1],
        parlist=parlist)$value,silent=silent)
      if (is.numeric(tmp)) {
        outs[i+badoffset] <- tmp
        break
      }
      # Try smart lower bound.
      if (bounds[i]<=0) {
        bounds[i]=max(
          c((parlist$b-0.98*parlist$A)/(max(mean(parlist$k*parlist$s),
          (parlist$k*parlist$s+2*sqrt(parlist$k)*parlist$s)[1])))
        )
		    next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
        bounds[i+1]=0.02*max(parlist$b)/(mean(parlist$k*parlist$s)-
          2*mean(sqrt(parlist$k)*parlist$s))
		    next
      }
      return(rep(0,length(outs)))
    }
  }
  if (doCDF) cumsum(outs) else outs
}



# ns=c(nsamp,nsamp);pi=p[isin,]
# onen=1e4; maxsamp=nsamp; warn=T
rfun=function(ns,pi,tlohi=c(0,1),onen=1e4,maxsamp=1e6,warn=T,sim.mult=NA,remove.outliers=TRUE) {
  # PARAMETERS
  #   ns: number of RTs for each response
  #   pi: data frame with columns ter, v sv, b and A (at least), one row per accumulator
  #   tlohi: lower and upper limits of uniform contaminant
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
    rts <- matrix(rgamma(nacc*onen,shape=pi$k,scale=pi$s),nrow=nacc)
    Ni <- dim(rts)[2] 
    rts <- (pi$b-runif(nacc*onen,0,pi$A))/rts # divide by v
    rts <- rts+pi$ter # will allow different ter for each response
    # pick the smallest
    win <- apply(rts,2,function(x){m=which.min(x);c(m,x[m])}) 
    # t0 noise
    if (any(names(pi)=="st0")) 
      win[2,] <- win[2,] + runif(Ni,0,pi$st0[1])
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
  names(Ns) <- names(ns)
  if (warn && nsamp > maxsamp)
    warning(paste("Unable to get sufficient samples. Required =",
                  paste(ns,collapse=" "),"Obtained =",paste(Ns,collapse=" ")))
  list(n=Ns,rt=cbind(contaminate=RT[1,],choice=RT[2,],rt=RT[3,]),sims=sim)
}


# sC1 = 1; Pc = .9; pow <- c(2,3,4)[c(Pc<.9,Pc>=.9 & Pc<=.95,Pc>.95)]; sE1 =sC1/(2*Pc)^pow 
# A=c(1,1); b=c(2,2); kC=c(2,2); kE=c(2,2); sE=c(sE1,sC1); sC=c(sC1,sE1); ter=c(.1,.1); st0=c(.5,.5)
# parlistC=list(A=A,b=b,k=kC,s=sC,ter=ter)
# parlistE=list(A=A,b=b,k=kE,s=sE,ter=ter)
# parlistCst=list(A=A,b=b,k=kC,s=sC,ter=ter,st0=st0)
# parlistEst=list(A=A,b=b,k=kE,s=sE,ter=ter,st0=st0)
# 
# pfun(t=Inf,parlistC)
# pfun(t=Inf,parlistC)+pfun(t=Inf,parlistE)
# pfun(t=Inf,parlistCst)
# pfun(t=Inf,parlistCst)+pfun(t=Inf,parlistEst)
# 
# x=c(1:700)/100; par(mfcol=c(2,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistE),type="l");abline(h=0)
# 
# x=c(1:700)/100; par(mfcol=c(1,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0); lines(x,dfun(t=x,parlistCst),col="red")
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0); lines(x,dfun(t=x,parlistEst),col="red")
# print(pfun(Inf,parlistC)+pfun(Inf,parlistE))
# print(pfun(Inf,parlistCst)+pfun(Inf,parlistEst))
# 
