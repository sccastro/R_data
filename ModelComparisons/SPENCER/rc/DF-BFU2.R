if (suppressWarnings(require(gsl,quietly=TRUE))) {
#    cat("Using GSL for incomplete gamma.\n")
} else {
  library(zipfR)
  gamma_inc=function(a,x,give=FALSE,strict=TRUE) Igamma(pmax(a,1e-10),x,lower=FALSE)
#    cat("Using zipfR for incomplete gamma.\n")
}
library(evd)

# Distribution random and start funcitons for ballistic frechet drift uniform 
# start point two choice model
#
# k=shape, s=scale


# generic function names start=sfun PDF=dfun, CDF=pfun, random sample=rfun
# sfun called by mt.R, dfun and pfun by fitting.R, and rfun by simdat in latter

sfun <-function(dati,pnams,pc=.01,k=2,aprior=.5,tertune=.9,Atune=.5,st0=.05) {
# Generic FU start routine for binary (correct/error) choice 
# Outputs standard b, A, k, s, ter parameterization by adding columns to
# design attribute of dati (data from one "rlevls" cell with both correct 
# and error data) assuming:
#  Pr(error)=(num_err+aprior)/(num_resp+2*aprior)
#  s set by v vs. 1-v mechanism
#  #  s = 2 produces RT'ish shapes
#  ter = tertune * min(rt)
#  b = .5 good for .75s second type means (VERY ROUGH)
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
    sC <- Pc; sE <- 1-Pc
    if (length(rt[as.logical(datii[[Cname]])])<2) 
      crt <- rt else crt <- rt[as.logical(datii[[Cname]])]
    if (length(crt)<2)
      stop("Less than 2 RTs in cell")
    kC <- k; kE <- k
    ter <- pmax(tertune*min(rt[rt>minter])-st0,minter)
    b <- .5
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

dfun=function(t,parlist,silent=TRUE,cv=NULL) {  
  # Defective PDF for responses on node 1 
  # st0>0: adds t0 variability convolution by numerical integration.
  # If length ter or st0 >1, first value used for all, other values ignored

  ################# SINGLE UNIT FUNCTIONS

  fptcdf=function(z,x0max,chi,driftrate,sddrift) {
    if (any(c(chi,chi-x0max,driftrate,sddrift)<0)) 
      return(rep(0,length(z)))#Protection for the pfrechet()
    ok <- z>0 & is.finite(z); out.value <- numeric(length(z)); out.value[z==Inf] <- 1
    if ( any(ok) ) {
      z <- z[ok]
      mew=1/driftrate; alpha=sddrift;
      min = (chi-x0max)/z; max = chi/z;
      pmax = pfrechet(max, loc=0, scale=1/mew, shape=alpha)
      pmin = pfrechet(min, loc=0, scale=1/mew, shape=alpha)
      zfrechet = (gamma_inc(1-(1/alpha),(mew*max)^(-alpha))-gamma_inc(1-(1/alpha),(mew*min)^(-alpha))) /
        (mew*(pmax-pmin))    
      term1 = ((z*zfrechet) - chi)/x0max ;term2 = (chi-x0max-(z*zfrechet))/x0max ; 
      out.value[ok] <- (1 + pmax*term1 + pmin*term2)
      out.value[z==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
      out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
    }
    return(out.value)
  }

  fptpdf=function(z,x0max,chi,driftrate,sddrift) {
    if (any(c(chi,chi-x0max,driftrate,sddrift)<0)) 
      return(rep(0,length(z))) #protection for pfrechet()
    ok <- z>0 & is.finite(z); out.value <- numeric(length(z))
    if ( any(ok) ) {
      z <- z[ok]
      mew=1/driftrate;  alpha=sddrift;
      min = (chi-x0max)/z; max = chi/z;
      Gmax  = pfrechet(max, loc=0, scale=1/mew, shape=alpha)
      Gmin  = pfrechet(min, loc=0, scale=1/mew, shape=alpha)
      D = Gmax - Gmin
      gam = gamma_inc(1-(1/alpha), (mew*max)^(-alpha))-gamma_inc(1-(1/alpha), (mew*min)^(-alpha))
      zfrechet <- gam/(mew*D)
      diffG1 = ((-chi/(z^2))*dfrechet(chi/z, loc=0, scale=1/mew, shape=alpha))
      diffG2 = ((-(chi-x0max)/(z^2))*dfrechet((chi-x0max)/z, loc=0, scale=1/mew, shape=alpha))    
      diffD = diffG1 - diffG2    
      diffgam = (-alpha*(((mew*chi)^(-alpha+1))/(z^(-alpha+2)))*exp(-(mew*chi/z)^(-alpha))) - 
                (-alpha*(((mew*(chi-x0max))^(-alpha+1))/(z^(-alpha+2)))*exp(-(mew*(chi-x0max)/z)^(-alpha)))
      diffzfrechet = (mew^(-1))*(((-D^(-2))*diffD)*gam + (diffgam*(D^(-1))))
      term1 = (Gmax - Gmin)*(zfrechet + (z*diffzfrechet))
      term2=diffG1*((zfrechet*z)-chi)
      term3 = diffG2*(chi-x0max-(zfrechet*z))
      out.value[ok] = ((term1+term2+term3)/x0max)
      out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
    }
    return(out.value)    
  }

  
  n1PDFfixedt0=function(t,A,b,k,s) {
    # Generates defective PDF for responses on node #1.
    N <- length(k) # Number of responses.
    if (N>2) {
      tmp=array(dim=c(length(t),N-1))
      for (i in 2:N) tmp[,i-1] <- 
        fptcdf(z=t,x0max=A[i],chi=b[i],driftrate=s[i],sddrift=k[i])
      G=apply(1-tmp,1,prod)
    } else {
      G=1-fptcdf(z=t,x0max=A[2],chi=b[2],driftrate=s[2],sddrift=k[2])
    }
    out <- G*fptpdf(z=t,x0max=A[1],chi=b[1],driftrate=s[1],sddrift=k[1])
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
    rts <- matrix(rfrechet(nacc*onen,loc=0,scale=pi$s,shape=pi$k),nrow=nacc)
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

# p=c(A=.5,b.sp=.8,b.acc=1.1,t0.sp=0.15,t0.acc=0.2,vc=2.5,ve=1.5,s.vc=2,s.ve=2)
# x=c(1:1500)/1000
# plot(x,fptpdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=p["s.vc"]),
#   type="l",ylab="Density",xlab="RT")
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=p["s.vc"]),
#   type="l",lty=2)
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=p["s.ve"]),
#   type="l", col="red")
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=p["s.ve"]),
#   type="l",lty=2,col="red")
# legend("topright",c("Speed Correct","Accuracy Correct","Speed Error","Accuracy Error"),
#   lty=c(1,2,1,2),col=c("black","black","red","red"))
# 
# plot(x,fptcdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=p["s.vc"]),
#   type="l",ylab="Density",xlab="RT",ylim=c(0,1))
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=p["s.vc"]),
#   type="l",lty=2)
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=p["s.ve"]),
#   type="l", col="red")
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=p["s.ve"]),
#   type="l",lty=2,col="red")
# legend("bottomright",c("Speed Correct","Accuracy Correct","Speed Error","Accuracy Error"),
#   lty=c(1,2,1,2),col=c("black","black","red","red"))

# big=.25
# fptcdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=p["s.vc"]) # 0.7281192
# fptcdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=p["s.vc"])# 0.4349157
# fptcdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=p["s.ve"]) # 0.413383
# fptcdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=p["s.ve"])# 0.1896178 
# 
# fptpdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=p["s.vc"]) #2.357421
# fptpdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=p["s.vc"])#2.489896
# fptpdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=p["s.ve"]) #2.28799
# fptpdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=p["s.ve"])#1.345095
# 
# p=c(A=.25,b.sp=.3,b.acc=.5,t0.sp=0.15,t0.acc=0.2,vc=.9,ve=.1,s.vc=2,s.ve=2)
# A=p[c("A","A")]; b=p[c("b.acc","b.acc")]; 
# sC=p[c("vc","ve")]; sE=p[c("ve","vc")];
# kC=p[c("s.vc","s.ve")]; kE=p[c("s.ve","s.vc")]
# ter=p[c("t0.acc","t0.acc")]; st0=c(.25,.25)
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
# x=c(1:1500)/1000; par(mfcol=c(2,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistE),type="l");abline(h=0)
# 
# x=c(1:1500)/1000; par(mfcol=c(1,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0); lines(x,dfun(t=x,parlistCst),col="red")
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0); lines(x,dfun(t=x,parlistEst),col="red")
# print(pfun(Inf,parlistC)+pfun(Inf,parlistE))
# print(pfun(Inf,parlistCst)+pfun(Inf,parlistEst))

