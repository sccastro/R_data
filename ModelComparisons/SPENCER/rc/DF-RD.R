# Distribution random and start funcitons for Ratcliff diffusion 
#
# generic function names start=sfun PDF=dfun, CDF=pfun, random sample=rfun
# sfun called by mt.R, dfun and pfun by fitting.R, and rfun by simdat in latter

os <- Sys.info()["sysname"]
if (os == "Windows") {
  dyn.load("rc/rat78pdf.dll")
  dyn.load("rc/Rfastdm.dll")
} else if (os == "SunOS") {
  dyn.load("rc/rat78pdf.sol.so")
  dyn.load("rc/Rfastdm.sol.so")
} else if (os == "Linux") {
#   dyn.load("rc/rat78pdf.so")
  dyn.load("rc/Rfastdm.so")
  dyn.load("rc/Rfastdm2.so")
  dyn.load("rc/meths.so")  
} else {
  dyn.load("rc/rat78pdf.osx.so")
  dyn.load("rc/Rfastdm.osx.so")
}
rm(os)

# t0tune=.9;aprior=.5;st=.05;sz_on_a=.2;sv_on_v=.2;pc=.01;pg=.5
sfun <-function(dati,pnams,t0tune=.9,aprior=.5,
  st=.05,sz_on_a=.2,sv_on_v=.2,pc=.01,pg=.5) {
  # Generic RD start routine for binary (correct/error) choice
  
  EZav = function(VRT,Pc,s=.1) {
    if (is.na(VRT)) VRT <- .01
    s2 <- .01
    L <- qlogis(Pc) 
    x <- L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
    v <- sign(Pc-0.5)*s*x^(1/4)
    if (v==0) 
      a <- (24*VRT*s2^2)^.25 else 
      a <- s2*qlogis(Pc)/v 
    c(a=a, v=v)
  }
  
  Di <- attr(dati,"D")
  mint0 <- as.numeric(dimnames(Di$nc)$minmax[1])*1.01
  # Define output
  cnams <- c("a","v","t0","z","sz","sv","st")
  ncols=length(cnams)
  if (any(pnams=="pc")) {
    cnams <- c(cnams,"pc")
    ncols <- ncols+1
  }
  if (any(pnams=="pg")) {
    cnams <- c(cnams,"pg")
    ncols <- ncols+1
  }
  RDstarts <- matrix(nrow=dim(Di$D)[1],ncol=ncols)
  dimnames(RDstarts) <- list(NULL,cnams)    
  # identify respone levels
  rlevs <- as.numeric(levels(dati$rcell))
  for (i in rlevs) {
    # get data for current response level
    datii <- dati[dati$rcell==i,]
    # get cells that correspond to current response level
    rli <- as.numeric(row.names(Di$D))[Di$D$rcell==i]
    # Use to get Pc and variance of correct RT
    Cname <- Di$SC[1] # ASSUMES FIRST SC SCORES CORRECT
    # ter
    rt <- datii[,Di$RT]
    tmp <- rt[rt>mint0]
    if (length(tmp)==0) t0 <- mint0 else 
      t0 <- pmax(t0tune*min(tmp),mint0)
    # a and v from EZ on robust p correct (Pc) and correct variance
    lcrt <- datii[as.logical(datii[,Cname]),Di$RT]-t0
    if (length(lcrt)<2) lcrt <- datii[,Di$RT]-t0  ###############################
    lcrt <- log(lcrt[lcrt>0])
    om <- exp(var(lcrt))
    e2mu <- exp(2*mean(lcrt))
    av <- EZav(VRT=e2mu*om*(om-1),
      Pc=(sum(as.logical(datii[,Cname]))+aprior)/(dim(datii)[1]+2*aprior))
    RDstarts[rli,"a"] <- rep(av[1],2)
    RDstarts[rli,"v"] <- rep(av[2],2)
    RDstarts[rli,"t0"] <- rep(t0,2)    
  }
  RDstarts[,"sv"] <- abs(RDstarts[,"v"])*sv_on_v      
  RDstarts[,"z"] <- RDstarts[,"a"]/2 # unbiaased      
  RDstarts[,"sz"] <- RDstarts[,"a"]*sz_on_a
  RDstarts[,"st"] <- st
  if (any(pnams=="pc")) 
    RDstarts[,"pc"] <- pc
  if (any(pnams=="pg")) 
    RDstarts[,"pg"] <- pg # probability guess LOWER 
#   print(RDstarts)
#   cat("\n")
  
  M2P(data.frame(RDstarts))
}

# TIMING TESTS OF DFUN
# c2 <- c(c2d10=1877.314,c2d1=4459.979,c2d.1=5907.575) # precision = 2
# c2.5 <- c(c2.5d1=4936.212,c2.5d1=8945.611,c2.5d=13893.247) # prec = 2.5

# C2 <- c(c2d10=2161.167,c2d1=5195.441,c2d.1=6022.825) # replication
# C2.5 <- c(c2.5d1=4920.598,c2.5d1=8956.897,c2.5d=15076.531)
# qtime <- 1064.592 # Rattiles prec=2.5 d=.1 old code

# round((c2+C2)/(2*qtime),1)
# c2d10  c2d1 c2d.1 
# 1.9   4.5   5.6 
# round((c2.5+C2.5)/(2*qtime),1)
# c2.5d1 c2.5d1  c2.5d 
# 4.6    8.4   13.6 

dfun=function(ts,parlist,cv=NULL,precision=2) {  
  # i =1 for lower boundary and i=2 for upper boundary
  
  density.RD <- function(t,p,i,precision=2) {
    # Call the C code
    densities <- vector(length=length(t))    
    output <- .C("dfastdm_b", 
      as.integer (length(t)),                                        # 1  IN:  number of densities
      as.vector  (p[c("a","v","ter","d","sz","sv","st","z")]), # 2  IN:  parameters
      as.vector  (t),                                                # 3  IN:  RTs
      as.double  (precision),                                        # 4  IN:  precision
      as.integer (i),                                                # 5  IN:  boundary 
      as.vector  (densities, mode="numeric")                         # 6 OUT:  densities
    )
    unlist(output[6])
  }
  
  out <- list(lower=numeric(0),upper=numeric(0))
  if (length(ts$lower)>0) out$lower <- rep(0,length(ts$lower))
  if (length(ts$upper)>0) out$upper <- rep(0,length(ts$upper))

  # deal with NA and NaN parameters
  if (any(is.na(unlist(parlist,use.names=F))) || 
      !all(is.finite(unlist(parlist,use.names=F)))) return(out)
  parlist[1,c("a","v","sz","sv","z")] <- parlist[1,c("a","v","sz","sv","z")]*10

  # LIMIT TO STOP CHECK IN C CODE FAILING
  maxsz <- min(c(1-parlist[1,"Z"],parlist[1,"Z"]))*parlist[1,"a"]*2
  if (parlist[1,"sz"]/maxsz > .9999)
    parlist[1,"sz"] <- .9999*maxsz
  if (parlist[1,"Z"] > .9999) 
    parlist[1,"Z"] <- .9999
  
  # SET SMALL VARIABILITY PARAMTERS TO ZERO AS QUICKER AND 
  # VERY SMALL SZ CAN CAUSE A CRASH (also maybe others but not observed yet)
  parlist[1,c("sz","sv","st")][parlist[1,c("sz","sv","st")]<1e-4] <- 0
    
  if (!any(names(parlist)=="d")) parlist$d <- 0 # new parameter in version 2 code
  
  # DEAL WITH HIGH t0 PROBLEMS
  t0 <- parlist$t0[1]
  ter <- t0+parlist$st[1]/2
  bad0 <- !(ts$lower - t0 > 0)
  bada <- !(ts$upper - t0 > 0)
  if (any(is.na(bad0)) || any(is.na(bada))) return(out)
  if ((length(bad0)==0) || all(bad0)) {
    nlo <- 0
    tlo <- numeric(0)
  } else {
    tlo <- ts$lower[!bad0]
    nlo <- sum(!bad0)
  }
  if ((length(bada)==0) || all(bada)) {
    nhi <- 0
    thi <- numeric(0)
  } else {
    thi <- ts$upper[!bada]
    nhi <- sum(!bada)
  }
  if (is.null(tlo) & is.null(thi)) return(out)
  # Note new code uses Z not z and sz as a proportion (i.e., actual sz/a)
  para <- c(a=parlist$a[1],v=parlist$v[1],ter=ter,d=parlist$d[1],
    sz=parlist$sz[1]/parlist$a[1],sv=parlist$sv[1],st=parlist$st[1],z=parlist$Z[1])
  if (length(out$lower)>0) out$lower <- 
    -density.RD(t=ts$lower,p=para,i=1,precision=precision)  
  if (length(out$upper)>0) out$upper <- 
    density.RD(t=ts$upper,p=para,i=2,precision=precision)
  out
}

#doCDF=F;precision=2.5
pfun=function(ts,parlist,doCDF=F,precision=2.5) {
  # Based on Voss and Voss RD cdf integrator, calling Rfastdm.so
  #   NB: $p is ZERO boundary completion probability  
  # NOTE: fastdm uses s=1 parameterization so non-ter parameters multiplied by 
  #       10 internally to fit with Ratcliff s=.1 parameterization.
  #  precision: controls precision of integration, larger numbers are MUCH slower
  #     default can be relaxed to 2 for quick exploration
  # Allows for lower bounds greater than lower quantile (returns zeros)

  out <- list(lower=numeric(0),upper=numeric(0))
  if (length(ts$lower)>0) out$lower <- rep(0,length(ts$lower)+1)
  if (length(ts$upper)>0) out$upper <- rep(0,length(ts$upper)+1)

  # deal with NA and NaN parameters
  if (any(is.na(unlist(parlist,use.names=F))) || 
      !all(is.finite(unlist(parlist,use.names=F)))) return(out)
  parlist[1,c("a","v","sz","sv","z")] <- parlist[1,c("a","v","sz","sv","z")]*10
  
  # SOME ARBITARY LIMITS TO STOP OBSERVED CRASHES OF C CODE
#  if (parlist[1,"a"]>100 | parlist[1,"sz"]>10) return(out)
  SZ <- parlist[1,"sz"]/parlist[1,"a"]
  parlist[1,"sz"][SZ>.99] <- parlist[1,"a"]*.99
  parlist[1,"sz"][SZ<.01] <- parlist[1,"a"]*.01
  
  # STOP VARIABILITY PARAMETERS FROM GETTING TOO SMALL
  parlist[1,c("a","sv","st")][parlist[1,c("a","sv","st")]<1e-6] <- 1e-6
  
  # DEAL WITH HIGH t0 PROBLEMS
  t0 <- parlist$t0[1]
  ter <- t0+parlist$st[1]/2
  bad0 <- !(ts$lower - t0 > 0)
  bada <- !(ts$upper - t0 > 0)
  if (any(is.na(bad0)) || any(is.na(bada))) return(out)
  if ((length(bad0)==0) || all(bad0)) {
    nlo <- 0
    tlo <- numeric(0)
  } else {
    tlo <- ts$lower[!bad0]
    nlo <- sum(!bad0)
  }
  if ((length(bada)==0) || all(bada)) {
    nhi <- 0
    thi <- numeric(0)
  } else {
    thi <- ts$upper[!bada]
    nhi <- sum(!bada)
  }
  if (is.null(tlo) & is.null(thi)) return(out)
  
  para <- c(parlist$a[1],parlist$v[1],ter,parlist$sz[1],
    parlist$sv[1],parlist$st[1],parlist$z[1])

#   print(para)
# print(parlist[1,"sv"])
#   print(nlo)
#   print(tlo)
#   print(nhi)
#   print(thi)
  
  tmp=try(.C("fastdmcdf",para=para,
         nhi=as.double(nhi),nlo=as.double(nlo),
         thi=thi,tlo=tlo,
         phi=numeric(nhi),plo=numeric(nlo),
         p=numeric(1),precision=precision))

  if (class(tmp)=="try-error") return(out)
  # NB numeric(0) returned in out$upper or out$lower if no observations in 
  # ts$upper or ts$lower respectivley. Should really be tmp$p, but
  # screws indexing and doenst mater in likelihood as muliplied by n=0
  
  if (any(!is.finite(tmp$plo)) | any(!is.finite(tmp$phi))) return(out)
  if (!is.finite(tmp$p) | is.na(tmp$p)) tmp$p <- 0 
  
  tmp$phi[tmp$phi>(1-tmp$p)] <- 1-tmp$p
  tmp$plo[tmp$plo>tmp$p] <- tmp$p
  
  if (!is.null(ts$lower)) {
    if (!is.null(tmp$plo)) out$lower[-length(out$lower)][!bad0] <- tmp$plo
    out$lower[length(out$lower)] <- tmp$p
    if (!doCDF) {
      out$lower <- diff(c(0,out$lower))
      out$lower[out$lower<0] <- 0
      osum <- sum(out$lower)
      if (osum>0) out$lower <- out$lower*tmp$p/osum
    }
  }
  if (!is.null(ts$upper)) {
    if (!is.null(tmp$phi)) out$upper[-length(out$upper)][!bada] <- tmp$phi
    out$upper[length(out$upper)] <- 1-tmp$p
    if (!doCDF) {
      out$upper <- diff(c(0,out$upper))
      out$upper[out$upper<0] <- 0
      osum <- sum(out$upper)
      if (osum>0) out$upper <- out$upper*(1-tmp$p)/osum
    }
  }
  out
}

# ns <- nsim; onen=1e4; warn=TRUE; roger.contaminate=FALSE
rfun=function(ns,pi,tlohi,onen=1e4,maxsamp=1e6,warn=TRUE,roger.contaminate=FALSE,
              sim.mult=NA,remove.outliers=TRUE) {
# PARAMETERS
#   ns: number of RTs for each response (in order 0 and a response)
#   pi: data frame with columns and one row
#   tlohi: lower and upper limits of uniform contaminant
#   onen: number to simulate per loop
#   maxsamp: limit on number of samples to take before giving up
# RETURN list with named elements 
#   rt: matrix [ns rows per choice,contaminate=0/1,choice=1:2,rt]
#   n: number samples for each response
  
#   nmc=10;v=.01;a=1;z=.5;sz=0;sv=0;t0=.1;st=.1;h=.001;s=.1
  mc.RD.C <- function(nmc,v,a,z=a/2,sz=0,sv=0,h=0.001,s=0.1,t0=0,st=0) {
    # nmc samples from Ratcliff diffusion by euler
    # h (step size default) OK ms accuracy except when very biased 
    # sz and st are FULL width of start point distribution
    # sv is SD of drift distribution, t0 is minimum RT
    # r=0 => 0 boundary termination, 1 -> a boundary termination
    
    zona <- z/a
    if (zona>.95 || zona<.05) h <- h/10
    sims=.C("euler",zmin=z-sz/2,zmax=z+sz/2,v=v,a=a,s=s,eta=sv,h=h,
            resp=numeric(nmc),rt=numeric(nmc),n=nmc)
    if (t0>0) sims$rt <- sims$rt + t0
    if (st>0) sims$rt <- sims$rt + runif(nmc,max=st)
    cbind(r=sims$resp,rt=sims$rt)
  }
  
  # SOME ARBITARY LIMITS TO STOP OBSERVED CRASHES OF C CODE
  SZ <- pi[1,"sz"]/pi[1,"a"]
  pi[1,"sz"][SZ>.99] <- pi[1,"a"][SZ>.99]*.99
  pi[1,"sz"][SZ<.01] <- pi[1,"a"][SZ<.01]*.01
  # STOP VARIABILITY PARAMETERS FROM GETTING TOO SMALL
  pi[1,c("a","sv","st")][pi[1,c("a","sv","st")]<1e-6] <- 1e-6
  
  nacc <- 2
  RT <- matrix(nrow=sum(ns),ncol=3)
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
    sims <- mc.RD.C(onen,
      v=pi$v,a=pi$a,z=pi$z,sz=pi$sz,sv=pi$sv,t0=pi$t0,st=pi$st)    
    if (remove.outliers) {
      bad <- sims[,"rt"] <= tlohi[1] | sims[,"rt"] >= tlohi[2]
      sims <- sims[!bad,]
    }    
    nsamp <- nsamp + dim(sims)[1] # number of samples taken    
    Ni <- dim(sims)[1]
    # contaminate
    if (any(names(pi)=="pc")) {
      contam <- rbinom(Ni,1,pi$pc[1])
      ncontam <- sum(contam)
      if (!roger.contaminate) { # random choice
        if (!any(names(pi)=="pg")) pg <- .5 else
          pg <- pi$pg[1] # lower contaminated
        sims[as.logical(contam),"r"] <- as.numeric(runif(ncontam)<pg)
      }
      sims[as.logical(contam),"rt"] <- 
        runif(ncontam,tlohi[1],tlohi[2])
    }  
    # number of each choice to add
    Nsi <- as.vector(table(factor(sims[,"r"],0:1))) 
    nadd <- ns-Ns # number left to reach ns
    toadd <- !logical(Ni) 
    for (i in 1:2) { # set toadd to T for first nadd choices
      isi <- sims[,"r"]==(i-1)
      if (any(isi)) toadd[isi] <- c(1:sum(isi))<=nadd[i]
    }
    endcol <- full2col+sum(toadd) # insert into RT
    if (endcol>full2col) {# Is there something to add?
      addrange <- (full2col+1):(endcol)
      addto <- cbind(sims[toadd,,drop=F],contam[toadd])
      if (class(try(RT[addrange,] <- addto))=="try-error") {
                print(dim(RT[addrange,,drop=F])); print(dim(addto))
                print(full2col); print(endcol); print(sum(toadd))
#                print(contam[toadd]); print(sims[toadd,])
      }       
      RT[addrange,] <- addto
    }
    full2col <- endcol
    Ns <- Ns + Nsi # number inserted
    if ( !is.null(sim) && (repn <= sim.mult) ) {
      sims <- data.frame(sims)
      names(sims) <- c("r","rt")
      sims$r <- factor(sims$r,levels=c(0,1),labels=names(ns))
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
  list(n=Ns,rt=cbind(contaminate=RT[,3],choice=RT[,1],rt=RT[,2]),sims=sim)
}


# pi=p[1,,drop=F]
# x=dat$RT[1:5]
# dfun(x,pi)
# pi=p[4,,drop=F]
# x=dat$RT[11:15]
# dfun(x,pi)

# pi <- data.frame(v=.1,a=.2,z=.1,sz=.1,sv=.1,t0=.2,st=.1,pc=0) 
# samp <- rfun(c(150000,500000),pi,c(0,3),maxsamp=3e6)
# x=c(2:400)/100
# tmp=pfun(list(lower=x,upper=x),pi,doCDF=T)
# lower=dfun(x,pi)
# pi$v=-pi$v; pi$z=pi$a-pi$z
# upper=dfun(x,pi)
# CDF by summing CUBA
# pupper=cumsum(upper)
# plower=cumsum(lower)
# # NORMALIZE to defective
# pupper=tmp$upper[length(x)+1]*pupper/pupper[length(pupper)]
# plower=tmp$lower[length(x)+1]*plower/plower[length(plower)]
# # PDF by differncing Voss ROUGH NORMALIZATION AT MAX CUBA
# dx <- (x[-length(x)]+x[-1])/2
# dtmp <- lapply(tmp,function(x){out=diff(x[-length(x)]); out=out/sum(out)})
# dtmp$upper=dtmp$upper*tmp$upper[length(x)+1]
# dtmp$lower=dtmp$lower*tmp$lower[length(x)+1]
# dtmp$upper=dtmp$upper*max(upper)/max(dtmp$upper)
# dtmp$lower=dtmp$lower*max(lower)/max(dtmp$lower)
# 
# # Compare CUBA and Voss CDF
# # Response probability estimates
# samp$n/sum(samp$n)                               # simulation
# c(tmp$lower[length(x)+1],tmp$upper[length(x)+1]) # Voss
# 
# par(mfcol=c(1,2))
# #CDF: black=Voss, red=CUBA
# plot(x,tmp$upper[-length(x)-1],type="l",xlab="t",ylab="CDF") 
# lines(x,tmp$lower[-length(x)-1],lty=2)
# abline(h=0)
# lines(x,pupper,col="red")
# lines(x,plower,lty=2,col="red")
# # #PDF: black=Cuba
# plot(x,upper,type="l",xlab="t",ylab="CDF") 
# lines(x,lower,lty=2)
# abline(h=0)
# lines(dx,dtmp$upper,col="red")
# lines(dx,dtmp$lower,col="red",lty=2)

# # COMPARE WITH SIMULATION
# par(mfcol=c(2,2))
# # CDF: black=sim, red=Voss
# p=c(1:99)/100  
# qsampu <- quantile(samp$rt[samp$rt[,"choice"]==1,"rt"],probs=p)
# qsampl <- quantile(samp$rt[samp$rt[,"choice"]==0,"rt"],probs=p)
# plot(qsampu,p*tmp$upper[length(x)+1],ylim=c(0,1),
#   type="l",xlab="t",ylab="CDF",xlim=c(.2,max(c(qsampu,qsampl))))
# lines(qsampl,p*tmp$lower[length(x)+1],lty=2)
# lines(x,tmp$upper[-length(x)-1],col="red")
# lines(x,tmp$lower[-length(x)-1],lty=2,col="red")

