
# Distribution random and start funcitons for ballistic lognormal drift uniform 
# start point two choice model
#
# with m=normal mean and s=normal sd

# generic function names start=sfun PDF=dfun, CDF=pfun, random sample=rfun
# sfun called by mt.R, dfun and pfun by fitting.R, and rfun by simdat in latter

#aprior=.5;minter=.1;tertune=.9
sfun <-function(dati,pnams,pc=.01,s=.5,
  aprior=.5,tertune=.9,Atune=.5,st0=.05) {
# Generic LU start routine for binary (correct/error) choice 
# Outputs standard b, A, k, s, ter parameterization by adding columns to
# design attribute of dati (data from one "rlevls" cell with both correct 
# and error data) assuming:
#  Pr(error)=(num_err+aprior)/(num_resp+2*aprior)
#  K set as Pc and 1-Pc
#  Use correct RT - ter to solve scale for both accumulators 
#  s = .5 produces RT'ish shapes
#  ter = tertune * min(rt)
#  b = 1 good for .75s means (VERY ROUGH)
#  A = b*Atune
  
  Di <- attr(dati,"D")
  minter <- as.numeric(dimnames(Di$nc)$minmax[1])*1.01
  # Define output
  ncols=5
  cnams <- c("b","A","m","s","ter")
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
    mC <- Pc; mE <- 1-Pc
    if (length(rt[as.logical(datii[[Cname]])])<2) 
      crt <- rt else crt <- rt[as.logical(datii[[Cname]])]
    if (length(crt)<2)
      stop("Less than 2 RTs in cell")
    sC <- sE <- s
    ter <- pmax(tertune*min(rt[rt>minter])-st0,minter)
    b <- 1
    A <- b*Atune
    nacc=length(cc)
    cc.cols <- cbind(rep(b,nacc),rep(A,nacc),rep(mC,nacc),
      rep(sC,nacc),rep(ter,nacc))
    nerr=length(ec)
    ec.cols <- cbind(rep(b,nerr),rep(A,nerr),rep(mE,nerr),
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



dfun=function(t,parlist,silent=TRUE,cv=NULL,precision=NULL) {  
  # Defective PDF for responses on node 1 
  # st0>0: adds t0 variability convolution by numerical integration.
  # If length ter or st0 >1, first value used for all, other values ignored

  ################# SINGLE UNIT FUNCTIONS

  fptcdf=function(z,x0max,chi,driftrate,sddrift) {
    
    ok <- z>0 & is.finite(z); out.value <- numeric(length(z)); out.value[z==Inf] <- 1
    if ( any(ok) ) {
      z <- z[ok]
      mean=driftrate ; sd=sddrift
      min = (chi-x0max)/z ; max = chi/z;
      zlognorm = 
        (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) /
        (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
      term1 = ((z*zlognorm) - chi)/x0max
      term2 = (chi-x0max-(z*zlognorm))/x0max 
      pmax = plnorm(max, meanlog=mean, sdlog=sd) 
      pmin = plnorm(min, meanlog=mean, sdlog=sd)    
      out.value[ok] = (1 + pmax*term1 + pmin*term2)
      out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
    }
    return(out.value)
  }
  
  fptpdf=function(z,x0max,chi,driftrate,sddrift) {
    ok <- z>0 & is.finite(z); out.value <- numeric(length(z))
    if (any(ok)) {
      z <- z[ok]
      mean=driftrate ; sd=sddrift;
      min = (chi-x0max)/z ; max = chi/z;
      zlognorm = 
        (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) /
          (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
      Gmax =plnorm(max,meanlog=mean,sdlog=sd) 
      Gmin =plnorm(min,meanlog=mean,sdlog=sd)    
      u = (pnorm((log(max)-mean-(sd)^2)/sd)-pnorm((log(min)-mean-(sd)^2)/sd))
      v = (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))    
      udash = (((-1/(sd*z))*dnorm((log(chi/z)-mean-(sd)^2)/sd)) - 
              ((-1/(sd*z))*dnorm((log((chi-x0max)/z)-mean-(sd)^2)/sd)))
      vdash = (((-1/(sd*z))*dnorm((log(chi/z)-mean)/sd)) - 
              ((-1/(sd*z))*dnorm((log((chi-x0max)/z)-mean)/sd)))
      const = exp(mean+((sd)^2)/2)    
      diffzlognorm = ((udash*v - vdash*u)/(v^2))*const #quotient rule
      term1 = (Gmax - Gmin)*(zlognorm + (z*diffzlognorm))
      term2=((-chi/(z^2))*dlnorm(chi/z,meanlog=mean,sdlog=sd))*((zlognorm*z)-chi)    
      term3 = (chi-x0max-(zlognorm*z))*
        ((-(chi-x0max)/(z^2))*dlnorm((chi-x0max)/z,meanlog=mean,sdlog=sd))
      out.value[ok]=((term1+term2+term3)/x0max)
      out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
    }
    return(out.value)  
  }
  
  n1PDFfixedt0=function(t,A,b,m,s) {
    # Generates defective PDF for responses on node #1.
    N <- length(m) # Number of responses.
    if (N>2) {
      tmp=array(dim=c(length(t),N-1))
      for (i in 2:N) tmp[,i-1] <- 
        fptcdf(z=t,x0max=A[i],chi=b[i],driftrate=m[i],sddrift=s[i])
      G=apply(1-tmp,1,prod)
    } else {
      G=1-fptcdf(z=t,x0max=A[2],chi=b[2],driftrate=m[2],sddrift=s[2])
    }
    out <- G*fptpdf(z=t,x0max=A[1],chi=b[1],driftrate=m[1],sddrift=s[1])
    out[t<=0 | !is.finite(out)] <- 0
    out
  }
  
    
  t <- t-parlist$ter[1]
  if (!any(names(parlist)=="st0") || parlist$st0[1]<.001) { 
    return(n1PDFfixedt0(t,
      A=parlist$A,b=parlist$b,m=parlist$m,s=parlist$s))
  } else {
    tmpf <- function(t,A,b,m,s,st0) 
      n1PDFfixedt0(t,A,b,m,s)/st0
    outs <- numeric(length(t))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=tmpf,
        lower=t[i]-parlist$st0[1],upper=t[i],
        A=parlist$A,b=parlist$b,m=parlist$m,s=parlist$s,
        st0=parlist$st0[1])$value,silent=silent)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}

pfun=function(t,parlist,doCDF=T,silent=TRUE,precision=NULL) {
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
        bounds[i]=max(c((parlist$b-0.98*parlist$A)/
  	      (max(mean(parlist$m),parlist$m[1])+2*parlist$s),1e-10))
		    next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
        bounds[i+1]=max(0.02*parlist$b/
          (mean(parlist$m)-2*parlist$s))
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
    rts <- matrix(rlnorm(nacc*onen,meanlog=pi$m,sdlog=pi$s),nrow=nacc)
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

# p=c(A=.1,b.sp=7,b.acc=10,t0.sp=0.3,t0.acc=0.5,vc=3.5,ve=2.8,s.ve=0.6)
# x=c(1:1500)/1000
# plot(x,fptpdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=1),
#   type="l",ylab="Density",xlab="RT")
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=1),
#   type="l",lty=2)
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=1),
#   type="l", col="red")
# lines(x,fptpdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=1),
#   type="l",lty=2,col="red")
# legend("topright",c("Speed Correct","Accuracy Correct","Speed Error","Accuracy Error"),
#   lty=c(1,2,1,2),col=c("black","black","red","red"))
# 
# plot(x,fptcdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=1),
#   type="l",ylab="Density",xlab="RT")
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=1),
#   type="l",lty=2)
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=1),
#   type="l", col="red")
# lines(x,fptcdf(z=x,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=1),
#   type="l",lty=2,col="red")
# legend("bottomright",c("Speed Correct","Accuracy Correct","Speed Error","Accuracy Error"),
#   lty=c(1,2,1,2),col=c("black","black","red","red"))
# 
# big=.25
# fptcdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=1) # 0.5694488 
# fptcdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=1)# 0.4270609
# fptcdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=1) # 0.2997836 
# fptcdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=1)# 0.1883862
# 
# fptpdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["vc"],sddrift=1) # 1.571514 
# fptpdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["vc"],sddrift=1)# 1.569017
# fptpdf(z=big,x0max=p["A"],chi=p["b.sp"],driftrate=p["ve"],sddrift=1) # 1.390304
# fptpdf(z=big,x0max=p["A"],chi=p["b.acc"],driftrate=p["ve"],sddrift=1)# 1.07977
# 
# 
# p=c(A=.1,b.sp=.5,b.acc=1,t0.sp=0.15,t0.acc=0.2,vc=.9,ve=.1,s.vc=.5,s.ve=.5)
# A=p[c("A","A")]; b=p[c("b.acc","b.acc")]; 
# mC=p[c("vc","ve")]; mE=p[c("ve","vc")];
# sC=p[c("s.vc","s.ve")]; sE=p[c("s.ve","s.vc")]
# ter=p[c("t0.acc","t0.acc")]; st0=c(.25,.25)
# parlistC=list(A=A,b=b,m=mC,s=sC,ter=ter)
# parlistE=list(A=A,b=b,m=mE,s=sE,ter=ter)
# parlistCst=list(A=A,b=b,m=mC,s=sC,ter=ter,st0=st0)
# parlistEst=list(A=A,b=b,m=mE,s=sE,ter=ter,st0=st0)
# 
# pfun(t=Inf,parlistC)
# pfun(t=Inf,parlistC)+pfun(t=Inf,parlistE)
# pfun(t=Inf,parlistCst)
# pfun(t=Inf,parlistCst)+pfun(t=Inf,parlistEst)
# 
# x=c(1:2000)/1000; par(mfcol=c(2,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistC),type="l");abline(h=0)
# plot(x,pfun(t=x,parlistE),type="l");abline(h=0)
# 
# x=c(1:2000)/1000; par(mfcol=c(1,2))
# plot(x,dfun(t=x,parlistC),type="l");abline(h=0); lines(x,dfun(t=x,parlistCst),col="red")
# plot(x,dfun(t=x,parlistE),type="l");abline(h=0); lines(x,dfun(t=x,parlistEst),col="red")
# print(pfun(Inf,parlistC)+pfun(Inf,parlistE))
# print(pfun(Inf,parlistCst)+pfun(Inf,parlistEst))
# 
