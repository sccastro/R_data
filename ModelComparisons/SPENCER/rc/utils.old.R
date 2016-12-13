################################################################################
###############################  LOADING AND EXTRACTION ########################

############################## load objects ####################################

get.fitlist <- function(dir.name,name.fitlist=NULL) {
  if (is.null(name.fitlist)) name.fitlist <- "fitlist"
  tmp=load(paste(paste(dir.name,name.fitlist,sep="/"),"RData",sep="."))
  get(tmp)
}

get.dat <- function(fitlist,subject=1,hmodel=1) {
  if (is.numeric(subject)) subject <- subject.names(fitlist)[subject]
  if (is.numeric(hmodel)) hmodel <- hmodel+1
  dat <- fitlist$dat
  qp <- attr(dat,"qp")
  dat <- dat[dat[[design.names(fitlist)$S]]==subject,]
  row.names(dat) <- NULL  
  attr(dat,"D") <- fitlist[[hmodel]]$M[[subject]]
  if (!is.null(qp)) attr(dat,"qp") <- qp
  dat
}

get.pmap <- function(fitlist=NULL,dir.name=NULL,model=NA,subject=1,hmodel=1) {
  if (is.null(fitlist) && is.null(dir.name)) 
    stop("Must supply either a fitlist or a directory name")
  if (is.null(fitlist)) fitlist <- get.fitlist(dir.name)
  if (is.numeric(hmodel)) {
    hmodel <- hmodel+1
    hmodel <- names(fitlist)[hmodel]
  }
  if (is.null(dir.name)) dir.name <- attr(fitlist[[hmodel]]$MT,"dname")
  if (is.null(dir.name)) stop("If no dir.name fitlist must have dname attribute")
  if (is.na(model)) model <- length(fitlist[[hmodel]]$MT[[1]])
  if (is.numeric(model)) model <- model.names(fitlist,hmodel)[model]
  if (is.numeric(subject)) subject <- subject.names(fitlist)[subject]
  load(paste(paste(dir.name,subject,"pmaps",model,sep="/"),"RData",sep="."))
  pmap
}

save.pmap <- function(pmap,model=NA,subject=1,fl=NULL,
                      dir.name=NULL,fitlist=NULL,hmodel=1) {

  if (is.null(fitlist) && is.null(dir.name)) 
    stop("Must supply either a fitlist or a directory name")
  if (is.null(fitlist)) fitlist <- get.fitlist(dir.name)
  if (is.numeric(hmodel)) {
    hmodel <- hmodel+1
    hmodel <- names(fitlist)[hmodel]
  }
  if (is.null(dir.name)) dir.name <- attr(fitlist[[hmodel]]$MT,"dname")
  if (is.null(dir.name)) stop("If no dir.name fitlist must have dname attribute")
  if (is.na(model)) model <- length(fitlist[[hmodel]]$MT[[1]])
  if (is.numeric(model)) model <- model.names(fitlist,hmodel)[model]
  if (is.numeric(subject)) subject <- subject.names(fitlist)[subject]
  tmp <- try(save(pmap,file=
    paste(paste(dir.name,subject,"pmaps",model,sep="/"),"RData",sep=".")))
  if (class(tmp)=="try-error")
    stop(paste("Unable to save model",model," for subject",subject))
}

get.fit <- function(fitlist=NULL,dir.name=NULL,model=NA,subject=1,hmodel=1) {
  if (is.null(fitlist) && is.null(dir.name)) 
    stop("Must supply either a fitlist or a directory name")
  if (is.null(fitlist)) fitlist <- get.fitlist(dir.name)
  if (is.numeric(hmodel)) {
    hmodel <- hmodel+1
    hmodel <- names(fitlist)[hmodel]
  }
  if (is.null(dir.name)) dir.name <- attr(fitlist[[hmodel]]$MT,"dname")
  if (is.null(dir.name)) stop("If no dir.name fitlist must have dname attribute")
  if (is.na(model)) model <- length(fitlist[[hmodel]]$MT[[1]])
  if (is.numeric(model)) model <- model.names(fitlist,hmodel)[model]
  if (is.numeric(subject)) subject <- subject.names(fitlist)[subject]
  load(paste(paste(dir.name,subject,"fits",model,sep="/"),"RData",sep="."))
  fit
}

############### return information about fitlist  #############################
  
subject.names <- function(fitlist) {
  names(attr(fitlist$dat,"D"))
}
  
design.data <- function(fitlist,subject=1) {
  attr(fitlist$dat,"D")[[subject]]$D
}
 
design <- function(fitlist,hmodel=1,subject=1) {
  if (is.numeric(hmodel)) hmodel <- hmodel+1
  fitlist[[hmodel]]$M[[subject]]$D
}
  
design.names <- function(fitlist,subject=1) {
  design.names <- attr(fitlist$dat,"D")[[subject]]
  out=vector(mode="list",length=0)
  for (i in names(design.names)) 
    if (i!="D" & i!="nc") out[[i]] <- design.names[[i]]
  out
}  

design.censor <- function(fitlist,subject=1) {
  attr(fitlist$dat,"D")[[subject]]$nc
}  

model.names <- function(fitlist,hmodel=1,subject=1) {
  if (is.numeric(hmodel)) hmodel <- hmodel+1
  names(fitlist[[hmodel]]$MT[[subject]])
}

##################### return information about pmap ############################

parameter.names.pmap <- function(pmap) {
  names(pmap)  
}

start.pmap <- function(pmap) {
  attributes(pmap)$start  
}
  
best.pmap <- function(pmap) {
  attributes(pmap)$best  
}


################################################################################
############# PARAMETER MAPPING  AND TRANSFORMATION BETWEEN LINEAR MODELS ######

fitp2pi <- function(pmapi,pp) {
  # transforms fitting to natural parameters for one parameter type
  switch(pmapi$partype,positive=pmin(exp(pp),1e308)+pmapi$bound[1],
    doublebound=plogis(pp)*pmapi$bound[2]+pmapi$bound[1],pp)
}

p2fitpi <- function(pmapi,pp) {
  # transforms natural to fitting parameters for one parameter type
  switch(pmapi$partype,
    positive=ifelse(pp>(pmapi$bound[1]+1e-6),log(pp-pmapi$bound[1]),log(1e-6)),
    doublebound=ifelse((pp<(pmapi$bound[1]+1e-6)),qlogis(1e-6),
      ifelse((pp>(pmapi$bound[2]-1e-6)),qlogis((pmapi$bound[2]-1e-6)/pmapi$bound[2]),
          suppressWarnings(qlogis((pp-pmapi$bound[1])/pmapi$bound[2])))),
    pp)
}

p2fitp <- function(pmap,p) {
  # transforms natural to fitting parameters for all parameter types
  for (i in names(pmap)) 
    p[[i]] <- p2fitpi(pmap[[i]],p[[i]])
  p
}

fitp2p <- function(pmap,p) {
  # transforms fitting to natural parameters for all parameter types  
  for (i in names(pmap)) 
    p[[i]] <- fitp2pi(pmap[[i]],p[[i]])
  p
}

#p=fitORp$fit$par; pmap=fitORp$pmap; trans=F
mapfitp <- function(p,pmap,trans=F) {
  # maps parameters to design matrix
  # used during fitting and in mapping start values
  pnams <- names(pmap)
  pout <- vector(mode="list",length=length(pmap))
  names(pout) <- pnams
  for ( i in pnams ) {
     if ( is.null(pmap[[i]]$consts) ) pp <- p[pmap[[i]]$pos] else {
      pp <- pmap[[i]]$consts
      if (pmap[[i]]$n!=0) pp[is.na(pp)] <- p[pmap[[i]]$pos]
    }
    pp <- pmap[[i]]$mm %*% pp
    if ( trans ) pp <- fitp2pi(pmap[[i]],pp)
    pout[[i]] <- pp
  }
  pout
}

#fitORp=pmaps[[maps[1]]]; newpmap=pmaps[[i]]
getstart=function(fitORp,newpmap) {
  # takes old parameters fitORp (either p, all parameters for full 
  # design on natural scale which are transformed to fit scale), 
  # or fit, list with $par parameters from previous fit on fitting 
  # scale which are expanded to full design with no transformation
  # using $pmap), and maps to fitting parameters using newpmap
  pnams <- names(newpmap)
  if ( is.data.frame(fitORp) ) p <- p2fitp(newpmap,fitORp) else 
     p <- mapfitp(fitORp$fit$par,fitORp$pmap,trans=F) 
  strt <- vector(length=0)
  snams <- strt
  for (i in pnams) {
    if ( newpmap[[i]]$n!=0 ) {
      newpar <- qr.solve(newpmap[[i]]$mm,p[[i]])
      if (!is.null(newpmap[[i]]$consts)) 
        ok <- is.na(newpmap[[i]]$consts) else
        ok <- rep(T,length(newpar))
      strt <- c(strt,newpar[ok])
      dn <- dimnames(newpmap[[i]]$mm)[[2]][ok]
      snams <- c(snams,paste(i,dn,sep=".")) 
    }
  }      
  names(strt) <- snams
  strt
}


################################################################################
######################### LOAD AND SAVE MODEL TREES   ##########################
################################################################################

# NB: PMAP stuff assumes fitlist has dname attribute ...

#loads a model from disk
get_pmap <- function(mname,ss,fl,hmname) {
  tmp <- try(load(
    paste(attr(fl[[hmname]]$MT,"dname"),"/",ss,"/pmaps/",mname,".RData",sep="")
  )) 
  if (class(tmp)=="try-error")
    stop(paste("Unable to find model",mname," for subject",ss))
  get(tmp)
}

save_pmap <- function(pmap,mname,ss,fl,hmname) {
  tmp <- try(save(pmap,file=
    paste(attr(fl[[hmname]]$MT,"dname"),"/",ss,"/pmaps/",mname,".RData",sep="")
  )) 
  if (class(tmp)=="try-error")
    stop(paste("Unable to save model",mname," for subject",ss))
}


# loads a model tree for one subject from disk
get_pmaps <- function(ss,fl,hmname) {
  dname <- paste(attr(fl[[hmname]]$MT,"dname"),"/",ss,"/pmaps",sep="")
  isdir <- file.info(dname)$isdir
  if (is.na(isdir) || !isdir)
    stop("Model tree directory not present")
  fns <- dir(dname)
  mns <- unlist(strsplit(fns,".RData"))
  names(fns) <- mns
  for (i in mns) {
    fl[[hmname]]$MT[[ss]][[i]] <- 
      get.pmap(fl,model=i,subject=ss,hmodel=hmname)
  }
  fl
}

save_pmaps <- function(ss,fl,hmname) {
  dname <- paste(attr(fl[[hmname]]$MT,"dname"),"/",ss,"/pmaps",sep="")
  isdir <- file.info(dname)$isdir
  if (is.na(isdir) || !isdir)
    stop("Model tree directory not present")
  for (i in names(fl[[hmname]]$MT[[ss]]))
    save_pmap(fl[[hmname]]$MT[[ss]][[i]],i,ss,fl,hmname)
}
  
  
get_mt <- function(fl,hmname) {
  
  snames <- names(fl[[hmname]]$M)
  cat(paste("Loading model tree for participant(s):\n"))   
  for (i in snames) {
    cat(paste(i,""))
    fl <- get_pmaps(i,fl,hmname)
  }
  cat("\n")
  fl
}


remove_mt <- function(fl,mhname) {
  snams <- names(fl[[mhname]]$MT)
  fl[[mhname]]$MT <- vector("list",length(snams))
  fl
}  

################################################################################
################################   GET RESULTS   ###############################
################################################################################


reduce.MT <- function(MT,hmodel=1) {
        
  keep <- function(MTi,hmodel) {    
    mnams <- names(MTi)
    if (is.list(hmodel)) {
      keep <- names(hmodel)  
    } else {
      keep <- hmodel                
      parents <- hmodel
      repeat {
        children <- character(0)                
        for (i in parents)
          children <- c(children,MTi[[i]])
        if (length(children)==0) break
        parents <- unique(children)
        keep <- unique(c(keep,parents))
      }
    }
    mnams[mnams %in% keep]
  }
   
  out <- MT
  nmod <- unlist(lapply(out,length))
  if (is.numeric(hmodel)) hmodel <- names(MT[[1]])[hmodel]
  if (length(hmodel)==1) {
    hmodel <- as.list(rep(hmodel,length(MT)))
    names(hmodel) <- names(MT)
  }
  for (i in names(MT)) { 
    is.in <- keep(MT[[i]],hmodel[[i]])
    outi <- vector(mode="list",length=length(is.in))
    names(outi) <- is.in
    for (j in names(outi)) outi[[j]] <- MT[[i]][[j]]
    out[[i]] <- outi
  }
  nmod1 <- unlist(lapply(out,length))
  is.less <- nmod > nmod1
  if (any(is.less)) attr(out,"reduce") <-
    rbind(possible=nmod,to.fit=nmod1,not.fit=nmod-nmod1)
  out
}


# fl=DF;type="bic";reduce=T;trans=T;snames=snams
getcellpars=function(fl,mhname,type=1,reduce=F,trans=F,snames=character(0)) {
  # get design cell parameters from flhm=fl[[mhname]], where hmname is a 
  #   character string specifying a model tree in fl.
  # The model extracted is specified by type, which can equal: 
  # 1) an interger specifying the models location in the pmaps list 
  #    (default 1 = TOP model), 
  # 2) by character string name or 
  # 3) BMA based on type="bic" or type="aic"
  #
  # reduce=F gives value for each design cell (NB ALWAYS trans=T under this option)
  # reduce=T gives value for each parameter in linear model in order specified
  #  in design so can plug in to pdf/cdf calculations
  # trans=T puts back on natural scale
  # snams: do for subset of subjects
 
  
  flhm <- fl[[mhname]]
  AA1 <- names(flhm$M[[1]]$D)[1]=="rcell" # No latents
  Ds <- attr(fl$dat,"D")
  dname <- attr(flhm$MT,"dname")
  if (length(snames)==0) snames <- names(flhm$MT)
  hmnames <- names(flhm$MT[[1]])
  is.bma <- any( type %in% c("bic","aic") )
  # get pmap for top model or selected model
  if (!is.bma) { 
    load(paste(paste(dname,snames[1],"pmaps",type,sep="/"),"RData",sep=".")) 
    if (is.numeric(type)) if (!(type %in% 1:length(hmnames)))
      stop(paste("Specified model not found")) else
      type <- hmnames[type]  
    if (!(type %in% hmnames)) stop(paste("Specified model not found"))
    load(paste(paste(dname,snames[1],"pmaps",type,sep="/"),"RData",sep="."))
  } else load(paste(paste(
    dname,snames[1],"pmaps",hmnames[1],sep="/"),"RData",sep="."))
  pnames=names(pmap)
  if (is.bma) {
    pars <- cbind.data.frame(s=flhm$fits$pars.names[,1],
      hmname=factor(flhm$fits$pars.names[,"hmname"],hmnames),
      fitpar=flhm$fits$pars[,"fitpar"])
    stats <- cbind.data.frame(s=flhm$fits$stats.names[,1],
      hmname=factor(flhm$fits$stats.names[,"hmname"],hmnames),
      flhm$fits$stats)
    weights <- tapply(stats[[type]],stats$s,function(x){
      exp(-(x-min(x))/2)/sum(exp(-(x-min(x))/2))})   
    for (i in names(weights)) names(weights[[i]]) <- hmnames
  } else {
  	is.model <- flhm$fits$pars.names[,"hmname"]==type
    pars=cbind.data.frame(s=flhm$fits$pars.names[is.model,1],
      fitpar=flhm$fits$pars[is.model,"fitpar"])
  }
  # get design
  if (!AA1) 
    designs <- lapply(flhm$M,function(x){
      if (all(x$CV=="")) x$D[,-c(1:2)] else
        cbind(x$D[,-c(1:2)],x$CV)
    }) else
    designs <- lapply(flhm$M,function(x){
      if (all(x$CV=="")) x$D[,-1] else
        cbind(x$D[,-1],x$CV)  
    })
  # get names of factors  
  if (reduce) {
    vnames <- lapply(pmap,function(x){
      as.character(unclass(attr(terms(formula(x)),"variables")))[-1]}) 
  } else {
    vnames <- vector(mode="list",length=length(pnames))
    names(vnames) <- pnames
    for (i in pnames) vnames[[i]] <- dimnames(designs[[1]])[[2]] 
  }
  if (!is.bma) { # get cellpars list and reduce
    cellpars <- vector(mode="list",length=length(snames))
    names(cellpars) <- snames
    for (i in snames) {
      load(paste(paste(dname,i,"pmaps",type,sep="/"),"RData",sep="."))      
      if (reduce) 
        cellpars[[i]] <- mapfitp(pars$fitpar[pars$s==i],pmap,trans=trans) else {
           cellpars[[i]] <- fitp2modelp(pars$fitpar[pars$s==i],pmap,flhm$M[[i]])
          if (!AA1) designs[[i]] <- designs[[i]][attr(flhm$M[[i]],"reorder")$D,]
        }
    }
    if (!reduce) { 
      for (i in snames) cellpars[[i]] <- lapply(cellpars[[i]],function(y){
        cbind.data.frame(designs[[i]],y)	
      })
    } else {
      out <- vector(mode="list",length=length(pnames))
      names(out) <- pnames
      for (i in snames) {
      	for (j in pnames) {
          if (length(vnames[[j]])==0) # no design factors
            out[[j]] <- mean(cellpars[[i]][[j]]) else {
            tmp <- tapply(cellpars[[i]][[j]],
              data.frame(designs[[i]][,vnames[[j]]]),mean)
            names(dimnames(tmp)) <- vnames[[j]]      
            out[[j]] <- arr2df(tmp)
          }      	  	
      	}
      	cellpars[[i]] <- out 
      }
    }  
  } else { # take parameter BMA, get cellpar list and reduce
    cellpars <- vector(mode="list",length=length(snames))
    names(cellpars) <- snames
    out <- vector(mode="list",length=length(pnames))
    names(out) <- pnames
    tmp <- vector(mode="list",length=length(hmnames))
    names(tmp) <- hmnames
    cat(paste("Getting BMA for reduce =",reduce,"and trans = ",trans,"\n"))
    for ( i in snames )  {      
      cat(".")
      for (j in hmnames) {
        load(paste(paste(dname,i,"pmaps",j,sep="/"),"RData",sep="."))      
        tmp[[j]] <- mapfitp(pars[pars$s==i & pars$hmname==j,"fitpar"],
          pmap=pmap,trans=trans) 
      }
      design <- designs[[i]]
      for (j in pnames) {
        if (length(hmnames)==1) 
          out[[j]] <- apply(sapply(hmnames,function(x){
            tmp[[x]][[j]]}),1,sum) else
          out[[j]] <- apply(sapply(hmnames,function(x){
            tmp[[x]][[j]]*weights[[i]][[x]]}),1,sum)
        if (length(vnames[[j]])==0) # no design factors
          out[[j]] <- mean(out[[j]]) else {
            tmp1 <- tapply(out[[j]],data.frame(design[,vnames[[j]]]),mean)
            names(dimnames(tmp1)) <- vnames[[j]]      
            out[[j]] <- arr2df(tmp1)
          }
      }
      cellpars[[i]] <- out
    }
    cat("\n")
  }
  # add subject column and cumulate over subjects
  for (i in snames) { 
    if ( i==snames[1] ) {
      celldf <- lapply(cellpars[[i]],function(x){
        if (is.null(dim(x))) 
          tmp <- cbind.data.frame(s=i,y=as.numeric(x)) else 
          tmp <- cbind.data.frame(s=rep(i,dim(x)[1]),x)
        tmp$s <- as.character(tmp$s)
        tmp$y=as.vector(tmp$y)
        tmp})	
    } else {
      for ( j in names(cellpars[[i]]) ) {
        if (is.null(dim(cellpars[[i]][[j]])))
          tmp <- cbind.data.frame(s=i,y=as.numeric(cellpars[[i]][[j]])) else 
          tmp <- cbind.data.frame(s=rep(i,dim(cellpars[[i]][[j]])[1]),
            cellpars[[i]][[j]])
        celldf[[j]] <- rbind(celldf[[j]],tmp)
      }		
    }
  } 
  if (!reduce) celldf <- lapply(celldf,function(x){x[!is.na(x$y),]})   
  celldf
}

check.nest=function(fitlist,nesttol=.01) {
  MT <- fitlist[[2]]$MT[[1]]
  fits <- fitlist[[2]]$fits
  hmname <- names(MT)
  bads.names <- matrix(ncol=3,nrow=0)
  bads <- matrix(ncol=3,nrow=0)
  for (i in hmname[(length(hmname)-1):1]) {
    maps <- MT[[i]]
    for (j in maps) {
      dd <- fits$stats[,"dev"][fits$stats.names[,"hmname"]==j]-
        fits$stats[,"dev"][fits$stats.names[,"hmname"]==i]
      bad <- dd < -nesttol
      if (any(bad)) {
        names(bad) <- fits$stats.names[,"s"][fits$stats.names[,"hmname"]==j]
        badones=fits$stats.names[,"hmname"]==i & 
          fits$stats.names[,"s"] %in% names(bad)[bad]
        bads.names=rbind(bads.names,
          cbind(fits$stats.names[badones,,drop=F],rep(j,sum(bad))))
        bads=rbind(bads,
          cbind(fits$stats[,"dev"][fits$stats.names[,"hmname"]==i][bad],
          fits$stats[,"dev"][fits$stats.names[,"hmname"]==j][bad],-dd[bad]))
      }  
    }
  }
  bads=cbind.data.frame(bads.names,bads)
  names(bads)=c("s","Nesting","Nested","dNesting","dNested","dDev")
  bads <- bads[do.call(order,bads[,c("Nesting","Nested","s","dDev")]),]
  row.names(bads) <- NULL
  bads
} 


# compare.fits <- function(dname1,dname2,mname=1) {
#   fl1 <- get.fitlist(dname1,strsplit(dname1,".",fixed=T)[[1]][1])
#   fl2 <- get.fitlist(dname2,strsplit(dname2,".",fixed=T)[[1]][1])
#   mnames <- names(fl1[[2]]$MT[[1]]) 
#   if (is.numeric(mname)) mname <- mnames[mname]
#   d1 <- fl1[[2]]$fits$stats[fl1[[2]]$fits$stats.names[,"hmname"]==mname,"dev"]
#   d2 <- fl2[[2]]$fits$stats[fl2[[2]]$fits$stats.names[,"hmname"]==mname,"dev"]
#   names(d1) <- fl1[[2]]$fits$stats.names[
#     fl1[[2]]$fits$stats.names[,"hmname"]==mname,"s"]
#   names(d2) <- fl2[[2]]$fits$stats.names[
#     fl2[[2]]$fits$stats.names[,"hmname"]==mname,"s"]
#   odr <- order(d1)
#   d1 <- d1[odr]
#   d2 <- d2[odr]
#   if (!all(names(d1)==names(d2)))
#     stop("Subject names do not match")
#   d1d2 <- d1-d2
#   d1d2 <- d1d2[order(d1d2)]
#   plot(1:length(d1d2),d1d2,xlab="Subjects",ylab="Deviance(Fit1 - Fit2)")
#   abline(0,0)
# }

# nesttol=.1; topmodel=1;get.bma=F; subjects=NA; only.start=F
# dname="szgo.LBA"
get.fits=function(dname,nesttol=.1,topmodel=1,get.bma=F,subjects=NA,only.start=F) {  
  # Two modes: if "fitstatus.RData" present in dname loads hierarchy
  #   of models that have been fit, reduces from tomodel and creates
  #   summary file in dname, including bma if get.bma=T
  #   If not, ignores topmodel and loads hierarchy reduced to contain
  #     only models that have saved fit files and ignores get.bma
  # dname = names of directory containing fit files
  # nesttol: tolerance on nesting check
  # topmodel: default implies all fits in hierarchy 
  # get.bma: get bma of fits

  getonemod=function(out,sname,hmname,fiti) {   
    npar=length(fiti$par)
    pest=fiti$pest
    noSE=any(is.na(fiti$se))
    pname=NULL
    ptype=NULL
    for (j in names(pest)) {
      if (noSE) {
        pname=c(pname,rep(j,length(pest[[j]])))
        ptypej=names(pest[[j]]) 
      } else {
      	pname=c(pname,rep(j,dim(pest[[j]])[2]))
      	ptypej=dimnames(pest[[j]])[[2]]
      }
      ptype=c(ptype,ptypej)
    } 
    if (noSE) pestval=unlist(pest) else
      pestval=unlist(lapply(pest,function(x){x[2,]}))
    npest=length(pestval)
    out$pars=rbind(out$pars,cbind(fitpar=fiti$par,se=fiti$se))      
    out$pests=c(out$pests,pestval)
    out$pars.names=rbind(out$pars.names,cbind(rep(sname,npar),
      hmname=rep(hmname,npar)))
    out$pests.names=rbind(out$pests.names,cbind(rep(sname,npest),
      hmname=rep(hmname,npest),pname=pname,ptype=ptype))
    out$stat=rbind(out$stat,c(fiti$dev,fiti$bic,fiti$aic,
      length(fiti$par),fiti$convergence))
    out$stat.names=rbind(out$stat.names,c(sname,hmname))
    out 
  }
  
 
#   doall=T;reduce=F;trans=F;snams=character(0)
#   DF=fitlist;mhname=mname; snams=names(ok)[ok]
  bma.get <- function(DF,mhname,dname,doall=T,reduce=F,trans=F,snams=character(0)) {
    # get BMA summaries
    cat("\nCalculating BMA (slow for big model trees and/or many participants)\n")
    cellpars <- vector(mode="list",length=2)
    names(cellpars) <- c("bic","aic")
    tmp1 <- matrix(ncol=2,nrow=2)
    mode(tmp1) <- "list"
    dimnames(tmp1) <- list(reduce=c("TRUE","FALSE"),trans=c("TRUE","FALSE"))
    cellpars$bic <- tmp1
    cellpars$aic <- tmp1
    cat("Calculating BIC based results\n")
#    attr(DF[[mhname]]$MT,"dname") <- dname
    if (doall || (reduce & trans)) cellpars$bic["TRUE","TRUE"][[1]] <- 
      getcellpars(DF,mhname,type="bic",reduce=T,trans=T,snams)
    if (doall || (reduce & !trans)) cellpars$bic["TRUE","FALSE"][[1]] <- 
      getcellpars(DF,mhname,type="bic",reduce=T,trans=F,snams)
    if (doall || (!reduce & trans)) cellpars$bic["FALSE","TRUE"][[1]] <- 
      getcellpars(DF,mhname,type="bic",reduce=F,trans=T,snams)
    if (doall || (!reduce & !trans)) cellpars$bic["FALSE","FALSE"][[1]] <- 
      getcellpars(DF,mhname,type="bic",reduce=F,trans=F,snams)
    cat("Calculating AIC based results\n")
    if (doall || (reduce & trans)) cellpars$aic["TRUE","TRUE"][[1]] <- 
      getcellpars(DF,mhname,type="aic",reduce=T,trans=T,snams)
    if (doall || (reduce & !trans)) cellpars$aic["TRUE","FALSE"][[1]] <-
      getcellpars(DF,mhname,type="aic",reduce=T,trans=F,snams)
    if (doall || (!reduce & trans)) cellpars$aic["FALSE","TRUE"][[1]] <- 
      getcellpars(DF,mhname,type="aic",reduce=F,trans=T,snams)
    if (doall || (!reduce & !trans)) cellpars$aic["FALSE","FALSE"][[1]] <- 
      getcellpars(DF,mhname,type="aic",reduce=F,trans=F,snams)
    cellpars
  }

  make.tofit <- function(dname,fitlist) {
    MT <- fitlist[[2]]$MT
    tofit <- vector(mode="list",length=length(MT))
    names(tofit) <- names(MT)    
    for (i in names(tofit)) {
      isdone <- matrix(unlist(strsplit(
        dir(paste(dname,i,"fits",sep="/")),".",T)),nrow=2)[1,]
      if (length(isdone)==0) tofit[[i]] <- list("Empty") else {
        tmp=vector(mode="list",length=length(isdone))
        names(tmp) <- names(MT[[i]])[names(MT[[i]]) %in% isdone]
        tofit[[i]] <- lapply(tmp,function(x){character(0)})
      }
    }
    tofit
  }

  # main body
  load(paste(dname,"fitlist.RData",sep="/"))
  mname <- names(fitlist)[2]
  attr(fitlist[[mname]]$MT,"dname") <- dname
  hmname <- names(fitlist[[mname]]$MT[[1]])  

  if ( any(dir(dname)=="fitstatus.RData") ) { # hierarchical fit results   
    load(paste(dname,"fitstatus.RData",sep="/")) 
    if (is.numeric(topmodel)) topmodel <- hmname[topmodel]
    if (!(topmodel %in% hmname))
      stop("Top model not in model names")
  } else {
    # single fit results
    tofit <- make.tofit(dname,fitlist)      
    topmodel <- tofit
    get.bma=F
  }
  
  # Limit focus to completed subjects
  if (only.start) 
    snams <- names(tofit)[
      unlist(lapply(tofit,function(x){!any(x=="start")}))] else 
    snams <- names(tofit)[
      lapply(tofit,function(x){sum(unlist(lapply(x,length)))})==0]
  if (!all(is.na(subjects))) {
    if (is.numeric(subjects)) subjects <- subject.names(fitlist)[subjects]
    if (!all(subjects %in% snams))
      stop(paste("Subject names not in set of completed subjects:\n",snams))
    snams <- snams[snams %in% subjects]
  }
  if (length(snams)==0)
    stop("Fits not available for any subject")
  # Limit focus from topmodel
  fitlist[[mname]]$MT <- reduce.MT(fitlist[[mname]]$MT,topmodel)
  hmname <- names(fitlist[[mname]]$MT[[1]]) 
  out=list(
    pars=matrix(nrow=0,ncol=2,
      dimnames=list(NULL,c("fitpar","se"))),
    pars.names=matrix(nrow=0,ncol=2,
      dimnames=list(NULL,c("s","hmname"))),
    pests=vector(length=0),
    pests.names=matrix(nrow=0,ncol=4,
      dimnames=list(NULL,c("s","hmname","pname","ptype"))),
    stat=matrix(nrow=0,ncol=5,dimnames=
      list(NULL,c("dev","bic","aic","npar","convergence"))),
    stat.names=matrix(nrow=0,ncol=2,dimnames=
      list(NULL,c("s","hmname")))
  )
  cat(paste("  Getting results for participant(s):\n"))
  for (ss in snams) {
    cat(paste(ss,""))
    for (i in hmname ) {
      load(paste(paste(dname,ss,"pmaps",i,sep="/"),"RData",sep="."))
      names(attr(pmap,"best")) <- names(attr(pmap,"start"))
      save(pmap,file=paste(paste(dname,ss,"pmaps",i,sep="/"),"RData",sep="."))
      tmp=try(load(paste(paste(dname,ss,"fits",i,sep="/"),"RData",sep=".")),silent=TRUE)
      if (class(tmp)=="try-error") {
        stop(paste("Unable to load fit",paste(paste(dname,ss,"fits",i,sep="/"),"RData",sep=".")))
      }
      out <- getonemod(out=out,sname=ss,fiti=fit,hmname=i)
    }   
  }
  cat("\n")
  row.names(out$pars)=NULL
  row.names(out$pars.names)=NULL
  fitlist[[mname]]$fits[["pars"]]=out$pars
  fitlist[[mname]]$fits[["pars.names"]]=out$pars.names 
  fitlist[[mname]]$fits[["pests"]]=out$pests
  fitlist[[mname]]$fits[["pests.names"]]=out$pests.names 
  fitlist[[mname]]$fits[["stats"]]=out$stat 
  fitlist[[mname]]$fits[["stats.names"]]=out$stat.names
  # Check fitlist$fits to see if complete
  ok <- vector(length=length(names(attr(fitlist$dat,"D"))))
  names(ok) <- names(attr(fitlist$dat,"D"))
  for (i in names(ok)) {
    ok[i] <- all( hmname %in% fitlist[[mname]]$fits[["stats.names"]]
      [fitlist[[mname]]$fits[["stats.names"]][,"s"]==i,"hmname"]) 
  }
  if (all(ok)) 
    cat("\nSet of fits complete\n") else    
    cat(paste("The fits for subject(s)",
      paste(names(ok)[!ok],collapse=", "),"are not complete"))
  if ( any(ok) & get.bma ) fitlist[[mname]]$fits$cellpars <- 
    bma.get(fitlist,mname,dname,snams=names(ok)[ok])
  flname <- strsplit(dname,".",fixed=T)[[1]][1]
  assign(flname,fitlist)
  save(list=flname,file=paste(paste(dname,flname,sep="/"),"RData",sep="."))
  nestfail=check.nest(fitlist,nesttol)
  cat("\n")
  nfail <- dim(nestfail)[1]
  if (nfail==0) cat("All fits are nested.\n") else {
    cat(paste(nfail,"fits are not nested.\n"))
  }  
  nestfail 
}



# dname="myer1v.LBA"; sname=1
check.one.nest <- function(dname,sname=1,nesttol=.1) {
  if (length(sname)>1)
    stop("Can only check one sname at a time")
  flname <- strsplit(dname,".",fixed=T)[[1]][1]
  mname <- strsplit(dname,".",fixed=T)[[1]][2]
  load(paste(paste(dname,flname,sep="/"),"RData",sep="."))
  fitlist <- get(flname)
  snams <- names(fitlist[[mname]]$M)
  if (is.numeric(sname)) sname <- snams[sname]
  if (!(sname %in% snams))
    stop("sname not in fitlist")
  hmname <- names(fitlist[[mname]]$MT[[1]]) 
  for (i in hmname ) {
    load(paste(paste(dname,sname,"fits",i,sep="/"),"RData",sep="."))
    is.dev <- fitlist[[mname]]$fits[["stats.names"]][,"s"]==sname &
      fitlist[[mname]]$fits[["stats.names"]][,"hmname"]==i
    fitlist[[mname]]$fits[["stats"]][,"dev"][is.dev] <- fit$dev
  }   
  nestfail <- check.nest(fitlist,nesttol)
  nestfail[nestfail$s==sname,]
}

# dname="kHo.LBA";model="B~lR*LO & A~1 & v~CNI*C & sv~CNI*C & ter~CE & pc~1";probs=c(1:9/10)
# subjects=NA; sim.mult=1000; exclude.subjects=NA  
get.qcaf <- function(dname,model=1,probs=c(1:9)/10,subjects=NA,
                     exclude.subjects=NA,sim.mult=1000) {
  
  qcaf <- function(best,probs,pmap,dat,sim.mult=1000) {
    
    D <- attr(dat,"D")
    cnam <- D$SC[1]
    rtnam <- D$RT
    rcnams <- levels(D$D$rcell)
    nc <- D$nc
    tlohi=as.numeric(dimnames(nc)$minmax)
    p <- fitp2modelp(best,pmap,D)
    # get prdicted RT via simulation
    dat <- dat[is.finite(dat[[rtnam]]),]
    nsamp = 2*round((length(probs)+1)*sim.mult/2)
    D <- D$D
    for (i in 1:length(rcnams)) {
      # correct so first v biggest
      isin <- as.character(D$rcell) == rcnams[i] & 
        as.logical(as.character(D[[cnam]]))
      tmp <- rfun(c(2*nsamp,2*nsamp),p[isin,],tlohi,maxsamp=nsamp,warn=F)
      tmp$rt <- tmp$rt[!is.na(tmp$rt[,"rt"]),]
      cqrt <- c(quantile(dat[[rtnam]][dat$rcell==i & as.logical(dat[[cnam]])],probs),Inf)
      oc <- tapply(as.logical(dat[dat$rcell==i,cnam]),
        cut(dat[dat$rcell==i,rtnam],c(-Inf,cqrt)),mean) 
      pq <- c(quantile(tmp$rt[,"rt"][tmp$rt[,"choice"]==1],probs),Inf)
      pc <- tapply(tmp$rt[,"choice"]==1,cut(tmp$rt[,"rt"],c(-Inf,cqrt)),mean)
      if (i==1) 
        out <- cbind(D[D$rcell==i,][1,][rep(1,length(cqrt)),-c(1:4)],p=c(probs,1),
          oq=cqrt,pq=pq,oc=oc,pc=pc) else
        out <- rbind(out,cbind(D[D$rcell==i,][1,][rep(1,length(cqrt)),-c(1:4)],
          p=c(probs,1),oq=cqrt,pq=pq,oc=oc,pc=pc))
    }
    out 
  }
  
  fitlist=get.fitlist(dname)
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  if (any(is.na(subjects))) subjects <- subject.names(fitlist) else
    if (is.numeric(subjects)) subjects <- subject.names(fitlist)[subjects]
  if (!is.na(exclude.subjects)) {
    if (!all(exclude.subjects %in% subjects))
      stop("Subject(s) to exclude do not exist")
    subjects <- subjects[!(subjects %in% exclude.subjects)]
  }
  if (is.numeric(model)) model <- model.names(fitlist)[model[1]]
  cat(paste("Simulating",length(subjects),"subjects using model",model,"\n"))
  job.ids <- numeric(length(subjects))
  names(job.ids) <- subjects
  for (s in subjects) { # start all subjects
    cat(paste(s,""))
    dat=get.dat(fitlist,s)
    pmap=get.pmap(fitlist,dname,model,s)
    best <- best.pmap(pmap)
    tmp <- qcaf(best=best,probs=probs,pmap=pmap,
                dat=dat,sim.mult=sim.mult)
    tmp <- cbind(s=rep(s,dim(tmp)[1]),tmp)
    if (s==subjects[1]) out <- tmp else out <- rbind(out,tmp) 
  }
  out$s <- factor(out$s)
  cat("\n")
  out
}


################################################################################
########################### GENERAL UTILITIES ##################################

todf=function(arr,dvname="",lname="") {
  if (is.null(dim(arr))) out=data.frame(y=unlist(arr)) else {
    dn <- dimnames(arr)
    nd <- length(dn)
    if (is.list(arr[1])) {
      lnams <- names(arr[1][[1]])
      lens <- apply(arr,1:nd,function(x){length(x[[1]])})
      if (any(lens[1]!=lens))
        stop("todf cannot handle arrays of unequal length lists")
      if (lname=="") lname <- "List"
      dns <- vector(mode="list",length=nd+1)
      dns[[1]] <- lnams
      for (i in 2:(nd+1)) dns[[i]] <- dn[[i-1]]
      names(dns)=c(lname,names(dn))
      arr <- array(unlist(arr),dim=c(lens[1],unlist(lapply(dn,length))),
        dimnames=dns)
#      arr <- aperm(arr,c(2:(nd+1),1))
      dn <- dimnames(arr)
      nd <- length(dn)
    }
    if (nd==1) {
      out=cbind.data.frame(factor(dn[[1]]),arr)
      names(out)=c(names(dn),"y")
      row.names(out)=NULL
    } else {
      tmp=vector(mode="list",length=nd)
      names(tmp)=names(dn)
      k=1
      for (j in names(dn)) {
        n=length(dn[[j]])
        tmp[[j]]=gl(n,k,length(arr),dn[[j]])
        k=k*n
      }
      out=cbind(data.frame(tmp),y=as.vector(arr))
      row.names(out)=NULL
    }
  }
  if (dvname!="") names(out)[names(out)=="y"] <- dvname
  out[,c((dim(out)[2]-1):1,dim(out)[2])]
}


arr2df=function(arr,dvname="") {
  if (is.null(dim(arr))) out=data.frame(y=arr) else {
    dn=dimnames(arr)
    if ( length(dn)==1 ) {
      out=cbind.data.frame(factor(dn[[1]],dn[[1]]),arr)
      names(out)=c(names(dn),"y")
      row.names(out)=NULL
    } else {
      tmp=vector(mode="list",length=length(dn))
      names(tmp)=names(dn)
      k=1
      for (j in names(dn)) {
        n=length(dn[[j]])
        tmp[[j]]=gl(n,k,length(arr),dn[[j]])
        k=k*n
      }
      out=cbind(data.frame(tmp),y=as.vector(arr))
      row.names(out)=NULL
    }
  }
  if (dvname!="") names(out)[names(out)=="y"] <- dvname
  if (dim(out)[2]==1) {
    out <- cbind(as.factor(row.names(out)),out)
    row.names(out)=NULL
  } 
  out
}

