################################## FITTING #####################################

#strt=strt[,1]      
multi.fit <- function(strt,dat,pmap,niter=NA,max.ntry=20,precision=2.5,
  use.nlm=FALSE,nlm.first=FALSE,ndigit=10,repeat.fit=.1,parscale=TRUE) {    

  try.fit <- function(par,fn,dat,pmap,niter,pu,strtobj,ndigit,
    use.nlm,qmpe,parscale,precision=2.5) {
    # fit with try protection, defaults to simplex if nlm fails
    # if simplex fails returns fit$convergence=NA, fit$value=Inf
      
    try.nlm=function(par,fn,dat,pmap,niter,pu,strtobj,ndigit,qmpe=FALSE,precision=2.5) {
      # nlm with try protection 
      if (is.na(niter)) niter=100 # default
      fit=try(nlm(p=par,f=fn,dat=dat,pmap=pmap,precision=precision,
        iterlim=niter,ndigit=ndigit,print.level=0,
        fscale=strtobj,typsize=par,pu=pu,qmpe=qmpe),silent=T)
#       nlm(p=par,f=fn,dat=dat,pmap=pmap,precision=precision,
#           iterlim=niter,ndigit=ndigit,print.level=1,
#           fscale=strtobj,typsize=par,pu=pu,qmpe=qmpe)
      if (class(fit)!="try-error") {
        fit$value=fit$minimum
        fit$par=fit$estimate
  	    fit$convergence=fit$code
        fit$method <- "nlm"
 	    } else {
 	      fit <- vector(mode="list",length=0)
 	      fit$convergence <- NA
 	      fit$value <- Inf
 	    }  
      fit
    }
        
    try.simplex=function(par,fn,dat,pmap,niter,pu,qmpe=FALSE,parscale=TRUE,precision=2.5) {
      # simplex with try protection 
      if (is.na(niter)) niter=500 # default
      if (parscale)
        fit=try(optim(par=par,fn=fn,dat=dat,pmap=pmap,method="Nelder-Mead",
          control=list(maxit=niter,parscale=pmax(abs(par),1e-3)),
          pu=pu,qmpe=qmpe,precision=precision),silent=T) else
        fit=try(optim(par=par,fn=fn,dat=dat,pmap=pmap,method="Nelder-Mead",
          control=list(maxit=niter),pu=pu,qmpe=qmpe,precision=precision),silent=T)

#       fit=optim(par=par,fn=fn,dat=dat,pmap=pmap,method="Nelder-Mead",
#         control=list(maxit=niter,trace=10,parscale=pmax(abs(par),1e-3)),
#         pu=pu,qmpe=qmpe,precision=precision)      
#       fit=optim(par=par,fn=fn,dat=dat,pmap=pmap,method="Nelder-Mead",
#                 control=list(maxit=niter,trace=10),pu=pu,qmpe=qmpe,precision=precision)      
      
      
      if (class(fit)=="try-error") {
        fit <- vector(mode="list",length=0)
        fit$convergence <- NA
        fit$value <- Inf
      } else fit$method <- "simplex" 
      fit
    }
    
    if (use.nlm) 
      fit <- try.nlm(par,fn,dat,pmap,niter,pu,strtobj,ndigit,qmpe=qmpe,precision=precision) 
    if (!use.nlm)
      fit <- try.simplex(par,fn,dat,pmap,niter,pu,qmpe=qmpe,parscale=parscale,precision=precision)
    fit
  }  
 
  # main body of multi.fit
  pu <- as.numeric(dimnames(attr(dat,"D")$nc)$minmax)
  qmpe <- !is.null(attr(dat,"qp"))
  strtobj <- objective(strt,dat,pmap,pu,qmpe,precision=precision)
  if (nlm.first | use.nlm) first.nlm <- TRUE else first.nlm <- FALSE
  fit <- try.fit(par=strt,fn=objective,dat=dat,pmap=pmap,niter=niter,precision=precision,
    strtobj=strtobj,ndigit=ndigit,use.nlm=first.nlm,pu=pu,qmpe=qmpe,parscale=parscale)
  if (repeat.fit[1]>0) { # more than one fit
    if (is.finite(fit$value) & !any(is.na(repeat.fit))) { # try more fits
      if ( length(repeat.fit)>1 ) {  # fixed number of trys
        for (i in repeat.fit) {
          fit1 <- try.fit(par=fit$par,fn=objective,dat=dat,pmap=pmap,qmpe=qmpe,precision=precision,
                    niter=i,strtobj=fit$value,ndigit=ndigit,use.nlm=use.nlm,pu=pu,parscale=parscale)
          if (fit1$value<fit$value) fit <- fit1
        }
      } else {  # repeat until improvement < repeat.fit
        ntry <- 1
        repeat {
          fit1 <- try.fit(par=fit$par,fn=objective,dat=dat,pmap=pmap,qmpe=qmpe,precision=precision,
                  niter=niter,strtobj=fit$value,ndigit=ndigit,use.nlm=use.nlm,pu=pu,parscale=parscale)        
          if (is.null(fit1$value)) break
          dfit <- fit$value-fit1$value
          if ( !is.na(dfit) && dfit>0 ) fit=fit1
          ntry <- ntry + 1
          if ( is.na(dfit) || dfit < repeat.fit || ntry > max.ntry) break
        }
      }  
    }
  }
  if (fit$value>strtobj) { # no improvement on start
    fit$value=strtobj
    fit$par=strt
  }
  fit
}


multi.start.fit=function(strt,pmap,dat,
    niter=NA,           # NA => default, if vector repeat.fit only first
    use.nlm=F,nlm.first=F,ndigit=12, # default for nlm method
    repeat.fit=.1,       # scaler> 0 => until change less, <= 0 do onece, vector=set of interations
    start.from=NA,
    get.se=FALSE,
    max.ntry=20,
    parscale=TRUE,
    precision=2.5) {
  # fits dat minimizing objective for model pmap from each start points in
  # each column of strt, using multi.start, then 
  # adds dev/aic/bic/hessian/pest 

  finalize.fit <- function(fit,pmap,dat,get.se=T,ndigit=10) {
  # add extras to fit with all starts done
        
    mapp=function(p,pmap,se) {
    # maps parameters as a list, used to collect fit results
    # puts +/- SE bounds, ALL ON FIT SCALE
      pnams=names(pmap)
      pout=vector(mode="list",length=length(pnams))
      names(pout)=pnams
      for (i in pnams ) {
        if ( is.null(pmap[[i]]$consts) )  {
          pp=p[pmap[[i]]$pos]
          if (!any(is.na(se))) {
            lowpp=pp-se[pmap[[i]]$pos]
            hipp=pp+se[pmap[[i]]$pos] 
          }
        } else {
          pp=pmap[[i]]$consts
          if ( pmap[[i]]$n!=0 ) pp[is.na(pp)]=p[pmap[[i]]$pos]
          if (!any(is.na(se))) {
            sei=pmap[[i]]$consts
            if ( pmap[[i]]$n!=0 ) sei[is.na(sei)]=se[pmap[[i]]$pos]
      	    lowpp=pp-sei
      	    hipp=pp+sei 
       	   lowpp[!is.na(pmap[[i]]$consts)]=NA   
     	     hipp[!is.na(pmap[[i]]$consts)]=NA   
    	     }    
        }
        pout[[i]]=pp
        names(pout[[i]])=dimnames(pmap[[i]]$mm)[[2]]
        if (!any(is.na(se))) pout[[i]]=rbind(LB=lowpp,
          Estimate=pout[[i]],UB=hipp)
      }
      pout
    }

    mml <- function(dat) {
    # minimum multinomial likelihood
      sum(dat$qn*log(unlist(
        tapply(dat$qp,dat$cell,function(x){
        x[length(x)]*diff(c(0,x[-length(x)],1))  
      }))))
    }

    fit$dev <- 2*fit$value
    if (!is.null(attr(dat,"qp"))) { # zero for prefect fit
      fit$dev <- fit$dev + 2*mml(dat)
    }
    fit$aic=fit$dev+2*length(fit$par)
    fit$bic=fit$dev+length(fit$par)*log(sum(dat$qn))
    pu <- as.numeric(dimnames(attr(dat,"D")$nc)$minmax)
    if (get.se) {
      hessian=try(nlm(p=fit$par,f=objective,dat=dat,pmap=pmap,pu=pu,
        iterlim=1,hessian=T,ndigit=ndigit,qmpe=!is.null(attr(dat,"qp")),
        fscale=fit$value,typsize=fit$par)$hessian,silent=T)
      if (class(hessian)=="try-error") fit$se=rep(NA,length(fit$par)) else {
        inv=try(chol2inv(hessian),silent=T)
        if (class(inv)=="try-error") fit$se=rep(NA,length(fit$par)) else
          fit$se=sqrt(diag(inv))/sqrt(sum(dat$qn)) 
      }
    } else fit$se=rep(NA,length(fit$par))
    fit$pest=mapp(fit$par,pmap,fit$se)
    fit
  }

  if ( !is.matrix(strt) ) 
    strt=matrix(strt,ncol=1,dimnames=list(names(strt),NULL))
  if ( attr(pmap,"npar")!=dim(strt)[1] ) stop(paste(
    "Length of \"strt\" does not match length implied by \"pmap\" (",
    attr(pmap,"npar")))  
  fit=list(value=Inf)
  for (i in 1:dim(strt)[2] ) {
    fiti=multi.fit(strt=strt[,i],dat=dat,pmap=pmap,niter=niter,parscale=parscale,precision=precision,
      use.nlm=use.nlm,nlm.first=nlm.first,ndigit=ndigit,repeat.fit=repeat.fit,max.ntry=max.ntry)
    if (fiti$value < fit$value) {
      fit=fiti
      if (!any(is.na(start.from))) fit$start.from <- start.from[i]
    }
  }
  finalize.fit(fit=fit,pmap=pmap,dat=dat,ndigit=ndigit,get.se=get.se)
}



# subjects=NA;topmodel=1;max.ntry=20
# start.all=FALSE;multi.start=FALSE
# best.start=FALSE;best.or.start=TRUE;only.start=FALSE
# niter=NA;use.nlm=FALSE;nlm.first=FALSE;ndigit=10;repeat.fit=.1
# est.time=60*60*10;nesttol=.1;get.se=F;get.bma=FALSE
# do.fit=TRUE;parscale=TRUE;update.available=FALSE
# hydra.max=128;refit=FALSE;refits=NULL; report.updates=FALSE; save.fitreport=FALSE
# 
hydra.hfits=function(dname,subjects=NA,topmodel=1,max.ntry=20,
  start.all=FALSE,multi.start=FALSE,
  best.start=FALSE,only.start=FALSE,
  niter=NA,use.nlm=FALSE,nlm.first=FALSE,ndigit=10,repeat.fit=.1,
  est.time=60*60*10,nesttol=.1,get.se=F,get.bma=FALSE,
  do.fit=TRUE,parscale=TRUE,update.available=FALSE,
  hydra.max=128,refit=FALSE,refits=NULL,report.updates=FALSE,
  save.fitreport=FALSE,precision=2.5) {             
  
  # HYDRA VERSION OF hfit for multiple subjects
  # dname: Fit directory  
  # subjects: Vector of subjects to fit (NA => fit all)
  # topmodel: top model to fit from
  # start.all: use start for all models or only first?
  # multi.start: if more than one start available do in one hydra?
  # best.start: use pmap "best" attribute to start if available
  # only.start: only run fits specified as "start" (could also be a best)
  # nesttol: tolerance for nesting failure reporting
  # est.time: sent by hydra to grid (seconds, default 10hrs)
  # do.fit: to fits or just report fit numbers
  # get.se: use nlm to get hessian SE
  # parscale: simplex parameter scale option
  # update.available: look in fit directory to check what has been done 
  #   rather than trusting what is stored in fitstatus
  # hydra.max: launch only this many before collect cycle
  # refit: refit and if better save
  # refits: specific subject names to be refit  
  # report.updates: report impovement of fit each time a better fit is found
  # save fitreport.txt: ONLY WHEN DEBUGGING, CAN GET VERY LARGE FOR LARGE FITS 
  #                     AS RESAVES ALL OF job.lookup AFTER EACH COLLECT 
  # precision=2.5 RD integrator precision
#   fpmap=pmap; start.best=best.start;start.best.or.start=best.or.start
  
  
  startmat <- function(startms,fpmap,i,fitlist,
    start.best=FALSE,start.best.or.start=FALSE) {
    # make start point matrix
    start <- NULL
    for (k in startms) { # loop over available starts
      if (k=="start") {
        if ( !start.best ) strt <- attr(fpmap,"start") else { 
          tmp <- attr(fpmap,"best")
          if (!is.null(tmp)) strt <- tmp else
            strt <- attr(fpmap,"start")  
        }    
      } else {  # get start point from previous fit
        fit <- vector(mode="list",length=2)
        names(fit) <- c("pmap","fit")
        fit$pmap <- get_pmap(k,i,fitlist,hmname)
        fit$fit$par <- attr(fit$pmap,"best")
        strt <- getstart(fit,fpmap)
      }
      if (k==startms[1])
        start <- matrix(strt,ncol=1,dimnames=list(names(strt),NULL)) else
        start <- cbind(start,strt) # more than one start for fit
    }
    start
  }
  
  report <- function() {
  # change job.lookup to formatted data frame
    out <- lapply(job.lookup,function(x){
      x$startms <- x$startms[1];x$start <- as.character(x$start); x})
    out=t(matrix(unlist(out,use.names=F),ncol=length(out)))
    dimnames(out)=list(names(job.lookup),names(job.lookup[[1]]))
    out=data.frame(out)
    out$done <- as.logical(out$done)
    out$time <- as.numeric(as.character(out$time))/60 # minutes
    out
  }
    
  # Main body
  # Check model directory
  if (only.start) hydra.max=Inf
  isdir <- file.info(dname)$isdir
  if (is.na(isdir) || !isdir)
    stop("Model tree directory not present") 
  
  # fitlist
  load(paste(dname,"/fitlist.RData",sep=""))
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  mnams <- names(fitlist[[hmname]]$MT[[1]])
  
  # topmodel
  if (is.numeric(topmodel)) topmodel <- mnams[topmodel]
  if (!(topmodel %in% mnams))
    stop("Top model not in model names")
  fitlist[[hmname]]$MT <- reduce.MT(fitlist[[hmname]]$MT,topmodel)  
  
  # Add in "start" as starting point model in list of models to be fit
  tofit1 <- fitlist[[hmname]]$MT
  if ( start.all | only.start ) { # Use "start" model with all models
    tofit1 <- lapply(tofit1,function(x){
      lapply(x,function(y){c("start",y)})}) 
  } else { # Use start with only simplest model
    for (i in names(tofit1)) tofit1[[i]][[length(tofit1[[i]])]] <- 
      c("start",tofit1[[i]][[length(tofit1[[i]])]])
  }
  
  # Subjects
  if (any(is.na(subjects))) ss <- names(fitlist[[hmname]]$M) else
    if (is.numeric(subjects)) ss <- names(fitlist[[hmname]]$M)[subjects] else
      ss <- subjects
  if (!any(ss %in% names(fitlist[[hmname]]$M)))
    stop("Subjects not present in fitlist")
  if (length(names(fitlist[[hmname]]$M))!=length(ss))
    cat(paste("Fitting a subset of participants:",
              length(ss),"out of",length(names(fitlist[[hmname]]$M))))
  
  # Make a list of available starting models assuming no fits yet
  available1 <- as.vector(rep("start",length(tofit1)),mode="list")
  names(available1) <- names(tofit1) 
  
  # Load/update fitstatus
  if ( file.exists(paste(dname,"/","fitstatus.RData",sep="")) ) { 
    # Loading gets tofit and available objects
    load(paste(dname,"/","fitstatus.RData",sep=""))
    if ( update.available ) {
      for ( i in names(available) )
        available[[i]] <- c("start",
          unlist(strsplit(dir(paste(dname,"/",i,"/fits",sep="")),".RData")) )
    }
    # remove from full set any that are already fit
    for ( i in names(tofit1) ) {
      keep <- !sapply( names(tofit1[[i]]),function(x){
        x %in% available[[i]]} 
      ) # models not in available, so keep in set to fit
      if ( !( refit & (i %in% ss) ) | ( i %in% refits ) ) { # remove fit
        for ( j in names(tofit1[[i]]) ) 
          if ( !keep[j] ) { # fully available 
            tofit1[[i]][[j]] <- character(0) 
          } else { # check if partly done (still starts to use)
            if (j %in% names(tofit[[i]])) # start not used
              tofit1[[i]][[j]] <- tofit[[i]][[j]] # so copy it back in
          }
      }
    }
    if (refit) for (i in ss) available[[i]] <- available1[[i]]
    for (i in refits) available[[i]] <- available1[[i]]
  } else available <- available1 # as not obtained from fitstatus
  tofit <- tofit1 # tofit1 has current work requried  

  if ( update.available ) {
    cat("Updated fitstatus.RData\n")
    save(available,tofit,file=paste(dname,"/","fitstatus.RData",sep=""))
  }
  # add in subjects requested now but not previously
#   sadd <- ss[!(ss %in% names(tofit))]
#   for (i in sadd) {
#     tofit[[i]] <- fitlist[[hmname]]$MT[[i]]
#     if (start.all) { 
#       tofit[[i]] <- lapply(tofit[[i]],function(y){c("start",y)}) 
#     } else {
#       tofit[[i]][[length(tofit[[i]])]] <- "start"
#     }
#     available[[i]] <- "start"
#   }
  tostart <- tofit # entries in tostart will be removed as fitting occurs 
  
  # Report amount of work to be done
  todo <- 0
  for (i in ss) todo <- todo + sum(unlist(lapply(tofit[[i]],length)))
  if (todo==0) {
    do.fit <- F
    warning("No fits left to do!") 
  } else {
    nfit <- numeric(length(ss))
    names(nfit) <- ss
    for (i in ss) 
      nfit[i] <- sum(unlist(lapply(tofit[[i]],function(x){length(x)!=0})))
    cat(paste("\nThere are",sum(nfit),"models left to fit requiring",
      todo,"fits\n\n"))
    cat(paste("Largest model set (with",
      nfit[which.max(nfit)],"models) to fit\n"))
    tmp <- unlist(lapply(
      tofit[[names(nfit[which.max(nfit)])]],
        function(x){length(x)!=0}))
    print(names(tmp)[tmp])    
  }
  
  if (do.fit) {
    # Initialize
    job.lookup <- vector(mode="list",length=0)
    jobs.collected=vector(length=0)
    if (only.start) {
      nfit <- 0
      for (i in ss) {
        for (j in names(tostart[[i]])) 
          nfit <- nfit + 
            length(tostart[[i]][[j]][tostart[[i]][[j]] %in% available[[i]]])
      }
      todo <- nfit
      cat("\nONLY STARTS BEING DONE,",nfit,"remaining to do\n")   
    }      
    n.running <- 0
    repeat {
      # Launch new jobs if any
      for (i in ss) if ( sum(unlist(lapply(tostart[[i]],length)))>0 & (hydra.max > n.running) ) { 
        # Get subject specific stuff
        D <- fitlist[[hmname]]$M[[i]]      
        for ( j in names(tostart[[i]]) ) { # Check all models
          # Are there any unusued starts available?
          startms <- tostart[[i]][[j]][tostart[[i]][[j]] %in% available[[i]]]
          if ( (length(startms)>0) & (hydra.max > n.running) ) { # some fits to start
            # Get model specific stuff 
            pmap <- get.pmap(fitlist,model=j,subject=i,hmodel=hmname)
            dat <- get.dat(fitlist,subject=i,hmodel=hmname) 
            start <- startmat(startms,pmap,i,fitlist,best.start,best.or.start)
            if ( !is.null(start) ) {# Launch hydra job
              if ( multi.start ) {
                n.running  <- n.running+1  
                job.id <- hydraRun(multi.start.fit,dat=dat,pmap=pmap,max.ntry=max.ntry,
                  strt=start,repeat.fit=repeat.fit,niter=niter,
                  use.nlm=use.nlm,nlm.first=nlm.first,start.from=startms,
                  est.time=est.time,get.se=get.se,parscale=parscale,precision=precision)
                #record job details
                job.lookup[[as.character(job.id)]] <- list(s=i,mnam=j,
                  startms=startms,nstart=length(startms),start=Sys.time(),
                  time=NA,done=F)
              } else {
                for (k in 1:dim(start)[2]) {
                  n.running  <- n.running+1
                  job.id <- hydraRun(multi.start.fit,dat=dat,pmap=pmap,max.ntry=max.ntry,
                    strt=start[,k,drop=F],repeat.fit=repeat.fit,parscale=parscale,
                    niter=niter,use.nlm=use.nlm,nlm.first=nlm.first,precision=precision,
                    start.from=startms[k],est.time=est.time,get.se=get.se)

#               fits=multi.start.fit(dat=dat,pmap=pmap,precision=precision,
#                  strt=start[,k,drop=F],repeat.fit=repeat.fit,niter=niter,
#                  use.nlm=use.nlm,nlm.first=nlm.first,start.from=startms,
#                  get.se=get.se,parscale=parscale,max.ntry=max.ntry)                                  

                  # record job details
                  job.lookup[[as.character(job.id)]] <- list(s=i,mnam=j,
                    startms=startms[k],nstart=1,start=Sys.time(),
                    time=NA,done=F)              
                }
              }
            }
            # Remove fit from tostart candidates
            tostart[[i]][[j]] <- 
              tostart[[i]][[j]][!(tostart[[i]][[j]] %in% startms)]
          }     
        } # all starts launched for current subject        
      } # all subjects launched
      # report
      if ( !only.start ) todo <- 
        sum(unlist(lapply(lapply(tostart,unlist),length))[ss])
      done <- unlist(lapply(job.lookup,function(x){x$done}))
      if (is.null(done)) done <- 0
      nfit <- sum(unlist(lapply(job.lookup,function(x){
        if (!x$done) x$nstart else 0
      })))
      nmod=0
      nsubj=0
      for (i in ss) {
        nmod <- nmod+length(unlist(tofit[[i]]))
        nsubj <- nsubj+as.numeric(length(unlist(tofit[[i]]))>0)
      }
      
      if ( !all(done) ) {
        if ( only.start ) 
          cat(paste("\nRUNNING",sum(!done),"jobs (",nfit,
            "fits )",sum(!done),"fits remaining to run")) else
          cat(paste("\nRUNNING",sum(!done),"jobs (",nfit,
            "fits )",todo,"fits remaining to run"))
        if (todo==0) cat("\n\n") else
          cat(paste(" (for ",nsubj,"subjects)\n\n"))
        # Update job time
        ctime <- Sys.time() 
        job.lookup <- lapply(job.lookup,function(x){
          if (!x$done) x$time <- 
            as.numeric(difftime(ctime,x$start,units="secs")/x$nstart); x
        })
        # collect results
        repeat {
          results <- hydraCollect()      
          # pause if no results available  
          if (length(results) == 0) 
            Sys.sleep(0.2) else break
        }
                
        # process results
        if ( length(results)>0 ) for (n in names(results) ) if (!is.null(job.lookup[[n]]$s)) {
          # get info about job
          job.id <- as.integer(n)
          i <- job.lookup[[n]]$s
          j <- job.lookup[[n]]$mnam
          startms <- job.lookup[[n]]$startms          
          # remove corresponding start model
          tofit[[i]][[j]] <- tofit[[i]][[j]][!(tofit[[i]][[j]] %in% startms)]
          # Add fits made available by this fit and not already in available
          if (!only.start) if ( length(tofit[[i]][[j]])==0 & !any(available[[i]]==j) ) 
            available[[i]] <- c(available[[i]],j)
          # record time to run job
          job.lookup[[n]]$time <- 
            as.numeric(difftime(Sys.time(),
            job.lookup[[n]]$start,units="secs")/
            job.lookup[[n]]$nstart) 
          # save results if better 
          fitname <- paste(dname,"/",i,"/fits/",j,".RData",sep="")
          if ( file.exists(fitname) ) load(fitname) else 
            fit <- list(value=Inf,start.from="NA")
          if ( any(names(results[[n]])=="value") && results[[n]]$value < fit$value ) { 
            if (report.updates) 
              cat(paste(fitname,"saved with new value",round(results[[n]]$value),
                "<",round(fit$value),"\n"))             
            fit <- results[[n]]
            pmap <- get.pmap(fitlist,model=j,subject=i,hmodel=hmname) 
            attr(pmap,"best") <- fit$par
            save_pmap(pmap,j,i,fitlist,hmname)
            save(fit,file=fitname)                         
          }
          # update status
          n.running <- n.running - 1
          if (only.start) todo <- todo - 1
        } # results all processed 
        # remove completed jobs
        job.lookup <- job.lookup[-c(1:length(job.lookup))[
          names(job.lookup) %in% names(results)]]
        
        # Save the fit-list
        save(available,tofit,
             file=paste(dname,"/","fitstatus.RData",sep=""))
        # Write a tab-delimited report
        if (save.fitreport) {    
          tmp <- report()
          write.table(tmp,quote=F,sep="\t",row.names=F,
            file=paste(dname,"/","fitreport.txt",sep=""))
        }
        
      } # end !all.done
      # Time to quit?
      if (!only.start) todo <- sum(unlist(lapply(lapply(tofit,unlist),length))[ss])
      if (todo==0) break
    }       
  }
}



# nestfail=NA;topmodel=1;nesttol=.1;subjects=NA;only.start=F
# niter=NA;repeat.fit=.1;ndigit=10;use.nlm=F;use.simplex=T;nlm.first=F;parscale=T
# hydra.max=128;max.ntry=20
# 
hydra.fix.fits <- function(dname,nestfail=NA,topmodel=1,nesttol=.1,subjects=NA,only.start=F,
  niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F,parscale=T,
  hydra.max=128,max.ntry=20,get.bma=TRUE,precision=2.5) {
  
  
  getobj <- function(strt,dat,pmap,pu,qmpe) {
    strtobj <- 2*objective(strt,dat,pmap,pu,qmpe)  
    if (qmpe) { # zero for prefect fit
      strtobj <- strtobj + 2*sum(dat$qn*log(unlist(
        tapply(dat$qp,dat$cell,function(x){
          x[length(x)]*diff(c(0,x[-length(x)],1))  
        }))))
    }
    strtobj
  }
  
  
  fitlist=get.fitlist(dname)
  if ( any(is.na(nestfail)) ) {
    if ( any(is.na(subjects)) ) {
      nestfail <- get.fits(dname,topmodel=topmodel,
        only.start=only.start,nesttol=nesttol)
      snames <- as.character(unique(nestfail[,"s"])) 
    } else {
      snames <- subject.names(fitlist)
      if (is.numeric(subjects)) subjects <- snames[subjects]
      if (!all(subjects %in% snames))
        stop("Subject names not in fitlist")
      nestfail <- get.fits(dname,subjects=subjects,topmodel=topmodel,
          only.start=only.start,nesttol=nesttol)
    }
  } else   snames <- as.character(unique(nestfail[,"s"]))
  if (dim(nestfail)[1]==0) return(nestfail)  
  cat(paste("\n\nRUNNING fix for",dim(nestfail)[1],"fits\n\n"))
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  mnames <- names(fitlist[[hmname]]$MT[[1]])
  lost.causes <- matrix(nrow=0,ncol=3,
    dimnames=list(NULL,c("s","Nesting","Nested")))
  is.running <- lost.causes
  job.ids <- numeric(0)
  repeat {
    nfail <- dim(nestfail)[1]
    if ( nfail>0 ) {
      for (i in 1:( min(c(hydra.max,nfail)) )) {
        subject <- as.character(nestfail[i,"s"])
        model <- as.character(nestfail[i,"Nesting"]) 
        start.model <- as.character(nestfail[i,"Nested"]) 
        dat <- get.dat(fitlist,subject)
        pmap <- get.pmap(fitlist,dname,model,subject)
        # Start model fit
        fitname.start <- paste(dname,"/",subject,"/fits/",start.model,".RData",sep="")
        load(fitname.start)
        fit.start <- fit
        # Old fit
        fitname <- paste(dname,"/",subject,"/fits/",model,".RData",sep="")
        load(fitname)
        fit.old <- fit
        # Get start point model parameters and map
        spmap <- get.pmap(fitlist,dname,start.model,subject)
        strt <- best.pmap(spmap)
        new.fit <- vector(mode="list",length=2)
        names(new.fit) <- c("pmap","fit")
        new.fit$pmap <- spmap
        new.fit$fit$par <- strt
        strt <- getstart(new.fit,pmap)
        # Get deviance for new start point
        pu <- as.numeric(dimnames(attr(dat,"D")$nc)$minmax)
        qmpe <- !is.null(attr(dat,"qp"))
        strtobj <- getobj(strt,dat,pmap,pu,qmpe)          
        if ( (strtobj-fit.start$dev) < nesttol ) spmap <- NULL  # dont refit start
        job.ids <- c(job.ids,hydraRun(fix.fit.hydra,dat=dat,parscale=parscale,max.ntry=max.ntry,
          pmap=pmap,strt=strt,spmap=spmap,sstrt=fit.start$par,start.model=start.model,
          fit.start.value=fit.start$value,fit.old.value=fit.old$value,precision=precision,
          fitname.start=fitname.start,fitname=fitname,model=model,subject=subject,
          nesttol=nesttol,niter=niter,repeat.fit=repeat.fit,ndigit=ndigit,
          use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)) 
        is.running <- rbind(is.running,c(subject,model,start.model))
        #         
#                 tmp=fix.fit.hydra(dat=dat,parscale=parscale,,max.ntry=max.ntry,
#                 pmap=pmap,strt=strt,spmap=spmap,sstrt=fit.start$par,start.model=start.model,
#                 fit.start.value=fit.start$value,fit.old.value=fit.old$value,precision=precision,
#                 fitname.start=fitname.start,fitname=fitname,model=model,subject=subject,
#                 nesttol=nesttol,niter=niter,repeat.fit=repeat.fit,ndigit=ndigit,
#                 use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first) 
        # 
      }
      nestfail <- nestfail[-c(1:(min(c(nfail,hydra.max)))),]
    }
    repeat { # collect results
      results <- hydraCollect()      
      if (length(results) == 0) Sys.sleep(0.2) else break
    }
    for (n in names(results)) {
      if (!is.null(results[[n]]$spmap)) {
        save_pmap(results[[n]]$spmap,results[[n]]$start.model,
                  results[[n]]$subject,fitlist,hmname)
        fit <- results[[n]]$fit.start
        save(fit,file=results[[n]]$fitname.start)
      }
      if (is.null(results[[n]]$pmap)) 
        lost.causes <- rbind(lost.causes,
                             c(results[[n]]$subject,results[[n]]$model,
                               results[[n]]$start.model)) else {
                                 save_pmap(results[[n]]$pmap,results[[n]]$model,
                                           results[[n]]$subject,fitlist,hmname)
                                 fit <- results[[n]]$new.fit
                                 save(fit,file=results[[n]]$fitname) 
                               }
      is.done <- is.running[,"s"]==results[[n]]$subject & 
        is.running[,"Nesting"]==results[[n]]$model &
        is.running[,"Nested"]==results[[n]]$start.model
      is.running <- is.running[!is.done,,drop=F]
    }
    if (dim(is.running)[1]!=0) {
      cat(paste("\n\nRUNNING fix for",dim(is.running)[1],"fits\n\n"))    
      sdone <- unique(unlist(lapply(results,function(x){x$subject})))
#      sdone <- snames[!(snames %in% unique(is.running[,"s"]))]
    } else sdone <- character(0)
    if (length(sdone)>0) {
      for (i in sdone) {
        nestfaili <- check.one.nest(dname,i)
        if (dim(nestfaili)[1]==0) snames <- snames[snames!=i]
        if (i==sdone[1])  nestfail <- nestfaili else
          nestfail <- rbind(nestfail,nestfaili)
      }
      if (dim(nestfail)[1]>0) {
        bad <- as.character(nestfail$s) %in% lost.causes[,"s"]  &
          as.character(nestfail$Nesting) %in% lost.causes[,"Nesting"] &
          as.character(nestfail$Nested) %in% lost.causes[,"Nested"]
        nestfail <- nestfail[!bad,]
      }
      if (dim(nestfail)[1]>0)
        cat(paste("\n\nPrepared for RUNNING ",dim(nestfail)[1],"new fits\n\n"))
    }
    if (dim(is.running)[1]==0 & dim(nestfail)[1]==0) break
  }    
  get.fits(dname,only.start=only.start,topmodel=topmodel,get.bma=get.bma)  
}


# model=1;subjects=NA;exclude.subjects=NA;maxsamp=1e6;sim.mult=100
# qps=c(.1,.3,.5,.7,.9); bests=vector(length=0)

# sim.mult=10000; qps=c(1:99)/100;maxsamp=1e7

hydra.simdat <- function(dname,model=1,subjects=NA,
  exclude.subjects=NA,maxsamp=1e6,
  sim.mult=100,qps=c(.1,.3,.5,.7,.9),bests=vector(length=0)) {

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
#    tmp <- vector(mode="list",length=length(subjects)); names(tmp) <- subjects  
  for (s in subjects) { # start all subjects
    cat(paste(s,""))
    dat=get.dat(fitlist,s)
    pmap=get.pmap(fitlist,dname,model,s)
    if (length(bests)==0) 
      best <- best.pmap(pmap) else best <- bests[[s]]
    job.ids[s] <- hydraRun(simdat,best=best,pmap=pmap,
       dat=dat,sim.mult=sim.mult,maxsamp=maxsamp,qps=qps)    
#       tmp[[s]] <- simdat(best=best,pmap=pmap,
#          dat=dat,sim.mult=sim.mult,qps=qps)
  }
  repeat {
    repeat { # collect results
      results <- hydraCollect()      
      if (length(results) == 0) Sys.sleep(0.2) else break
    }
    for (n in names(results)) {
      job.id <- as.integer(n)
      s <- names(job.ids)[job.ids==job.id]
      sim <- results[[n]]$dat
      save(sim,file=paste(dname,"/",s,"/sims/",model,".RData",sep="")) 
      if (!is.null(results[[n]]$sim)) {
        sims <- results[[n]]$sim
        save(sims,file=paste(dname,"/",s,"/sims/",model,"sims.RData",sep="")) 
      }
      job.ids[s] <- 0
      cat(paste("Collected participant",s,
        "(",sum(job.ids!=0),"left to go )\n"))
    }
    if (all(job.ids==0)) break
  }
  cat("\n")
}

# model=1
# subjects=NA; est.time=60*60*10;sim.mult=10
sim.dat <- function(dname,model=1,subjects=NA,
  est.time=60*60*10,sim.mult=100) {

  fitlist=get.fitlist(dname)
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  if (any(is.na(subjects))) subjects <- subject.names(fitlist) else
    if (is.numeric(subjects)) subjects <- subject.names(fitlist)[subjects]
  if (is.numeric(model)) model <- model.names(fitlist)[model[1]]
  cat(paste("Simulating",length(subjects),"subjects using model",model,"\n"))
  for (s in subjects) { 
    dat=get.dat(fitlist,s)
    pmap=get.pmap(fitlist,dname,model,s)
    best <- best.pmap(pmap)
    if (is.null(best))
      cat(paste("Fit not available for subject",s,"\n")) else {
      cat(paste("Simulating subject",s,"\n"))
      datsim <- simdat(best=best,pmap=pmap,dat=dat,sim.mult=sim.mult)
      sim <- datsim$dat
      save(sim,file=paste(dname,"/",s,"/sims/",model,".RData",sep=""))
      if (!is.null(datsim$sim)) {
        sims <- datsim$sim
        save(sims,file=paste(dname,"/",s,"/sims/",model,"sims.RData",sep="")) 
      }
    }
  }
  cat("\n")
}


# pmap=fpmap;dat=datj; ndigit=12
fit.sim <- function(strt,pmap,dat,
                    niter=NA,           # NA => default, if vector repeat.fit only first
                    use.nlm=F,nlm.first=F,ndigit=12, # default for nlm method
                    repeat.fit=-1, # scaler> 0 => until change less, <= 0 do onece, vector=set of interations
                    max.ntry=1,
                    parscale=TRUE,
                    save.sim=TRUE){
  
  mml <- function(dat) {
    # minimum multinomial likelihood
    sum(dat$qn*log(unlist(
      tapply(dat$qp,dat$cell,function(x){
        x[length(x)]*diff(c(0,x[-length(x)],1))  
      }))))
  }
  
  fit=multi.fit(strt=strt,dat=dat,pmap=pmap,niter=niter,parscale=parscale,
    use.nlm=use.nlm,nlm.first=nlm.first,ndigit=ndigit,repeat.fit=repeat.fit,
    max.ntry=max.ntry)
  fit$dev <- 2*fit$value
  if (!is.null(attr(dat,"qp"))) { # zero for prefect fit
    fit$dev <- fit$dev + 2*mml(dat)
  }
  names(fit$par) <- names(strt)
  if (save.sim) sim <- sim1dat(best=fit$par,pmap=pmap,dat=dat) else sim=NULL
  list(fit=fit,sim=sim)
}


# subjects=NA;max.ntry=1;repeat.fit=-1;niter=NA;use.nlm=FALSE;nlm.first=FALSE;ndigit=10;
# est.time=60*60*10;do.fit=TRUE;parscale=TRUE;update.available=TRUE;hydra.max=8;refit=FALSE             
# save.sim=TRUE
# 
# refit <- TRUE;
# smodel="a~E & v~E*W & sv~E*W & Z~E & SZ~E & t0~E & st~E & pc~1"
# fmodel="a~E & v~E*W & sv~E*W & Z~E & SZ~E & t0~E & st~E & pc~1"
 
hydra.sfits=function(dname,subjects=NA,
                     smodel=1,fmodel=1,        # simulation data model, model to fit name
                     max.ntry=1,repeat.fit=-1, # default fit only once
                     niter=NA,use.nlm=FALSE,nlm.first=FALSE,ndigit=10,
                     est.time=60*60*10,
                     do.fit=TRUE,parscale=TRUE,update.savailable=TRUE,
                     hydra.max=128,refit=FALSE,save.sim=TRUE) {             
  
  # FIT SIMULATED DATA FROM MODELNAMEsims.RData FILES IN EACH SUBJCTS SIMS DIRECTORY
  # dname: Fit directory  
  # subjects: Vector of subjects to fit (NA => fit all)
  # est.time: sent by hydra to grid (seconds, default 10hrs)
  # do.fit: to fits or just report fit numbers
  # parscale: simplex parameter scale option
  # update.available: look in fit directory to check what has been done 
  #   rather than trusting what is stored in fitstatus
  # hydra.max: launch only this many before collect cycle
  # refit: ignore stored solutions and refit
  
  #   fpmap=pmap; start.best=best.start;start.best.or.start=best.or.start
  
  startmat <- function(startms,fpmap,i,fitlist,
                       start.best=F,start.best.or.start=F) {
    # make start point matrix
    start <- NULL
    for (k in startms) { # loop over available starts
      if (k=="start") {
        if (!(start.best | start.best.or.start) )
          strt <- attr(fpmap,"start") else { 
            tmp <- attr(fpmap,"best")
            if (!is.null(tmp)) strt <- tmp
            if (start.best.or.start & is.null(tmp)) 
              strt <- attr(fpmap,"start")  
          }    
      } else {  # get start point from previous fit
        fit <- vector(mode="list",length=2)
        names(fit) <- c("pmap","fit")
        fit$pmap <- get_pmap(k,i,fitlist,hmname)
        fit$fit$par <- attr(fit$pmap,"best")
        strt <- getstart(fit,fpmap)
      }
      if (k==startms[1])
        start <- matrix(strt,ncol=1,dimnames=list(names(strt),NULL)) else
          start <- cbind(start,strt) # more than one start for fit
    }
    start
  }
  
 
  # Main body
  # Check model directory
  isdir <- file.info(dname)$isdir
  if (is.na(isdir) || !isdir)
    stop("Model tree directory not present") 
  
  # fitlist
  load(paste(dname,"/fitlist.RData",sep=""))
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  mnams <- names(fitlist[[hmname]]$MT[[1]])
  
  # topmodel
  if (is.numeric(smodel)) smodel <- mnams[smodel]
  if (is.numeric(fmodel)) fmodel <- mnams[fmodel]
  if (!(smodel %in% mnams))
    stop("Simulated model not in model names")
  if (!(fmodel %in% mnams))
    stop("Model to fit not in model names")
  
  # Subjects
  if (any(is.na(subjects))) ss <- names(fitlist[[hmname]]$M) else
    if (is.numeric(subjects)) ss <- names(fitlist[[hmname]]$M)[subjects] else
      ss <- subjects
  if (!any(ss %in% names(fitlist[[hmname]]$M)))
    stop("Subjects not present in fitlist")
  if (length(names(fitlist[[hmname]]$M))!=length(ss))
    cat(paste("Fitting a subset of participants:",
              length(ss),"out of",length(names(fitlist[[hmname]]$M))))  
  # make a list of fits to do
  tofit <- vector(mode="list",length=length(ss))
  names(tofit) <- ss
  sf <- tofit; ff <- tofit
  for (i in ss) {
    sf[[i]] <- paste(dname,"/",i,"/","sims/",smodel,"sims.RData",sep="")
    if (file.exists(sf[[i]])) load(sf[[i]]) else
      stop(paste(sf[[i]],"simulated data file missing"))
    tofit[[i]] <- vector(mode="list",length=length(sims))
    nfit <- length(sims[[1]])
    tofit[[i]] <- 1:nfit
    ff[[i]] <- paste(dname,"/",i,"/","sims/",fmodel,"=",smodel,".RData",sep="")
    if ( !refit && file.exists(ff[[i]]) ) {
        load(ff[[i]])
#         
#         print(i)
#         print(unlist(lapply(sims,length)))
#         
        todo <- is.na(attr(sims,"dev"))
        tofit[[i]] <- c(1:nfit)[todo]
    } else { # make one ready for fits
      sims <- lapply(sims,function(x){lapply(x,function(y){NULL})})
      attr(sims,"dev") <- rep(NA,length(sims[[1]]))
      attr(sims,"best") <- vector(mode="list",length=nfit)
      save(sims,file=ff[[i]])     
    }
  }
  cat(paste(sum(unlist(lapply(tofit,length))),"fits left to do for",
            sum(unlist(lapply(tofit,length))>0),"subjects"))
  if ( do.fit & sum(unlist(lapply(tofit,length)))>0 ) {
    # Initialize
    job.lookup <- vector(mode="list",length=0)
    n.running <- 0
    tostart <- tofit
    current.i <- 1
    repeat { # launch subjects and collect
      if ( current.i <= length(ss) ) repeat { # launch fits 
        i <- ss[current.i]
        # Get subject specific stuff
        D <- fitlist[[hmname]]$M[[i]]      
        spmap <- get.pmap(fitlist,model=smodel,subject=i,hmodel=hmname)
        fpmap <- get.pmap(fitlist,model=fmodel,subject=i,hmodel=hmname)
        dat <- get.dat(fitlist,subject=i,hmodel=hmname) 
        load(sf[[i]])  
        for ( j in tostart[[i]] ) { 
          datj <- makedat(sim=lapply(sims,function(x){x[[j]]}),
            D=D,qp=attr(dat,"qp"),ss=ss,ssi=i)
          attr(datj,"D") <- D
          attr(datj,"qp") <- attr(dat,"qp")
          tmp <- vector(mode="list",length=2)
          names(tmp) <- c("pmap","fit")
          tmp$pmap <- spmap
          tmp$fit$par <- attr(tmp$pmap,"best")
          strt <- getstart(tmp,fpmap)
#           start <- matrix(strt,ncol=1,dimnames=list(names(strt),NULL))
          n.running  <- n.running+1
          job.id <- hydraRun(fit.sim,strt=strt,pmap=fpmap,dat=datj,
                             niter=niter,use.nlm=use.nlm,nlm.first=nlm.first,
                             repeat.fit=repeat.fit,max.ntry=max.ntry,
                             parscale=parscale,save.sim=save.sim)          
#           fits=fit.sim(strt=strt,pmap=fpmap,dat=datj,
#                        niter=niter,use.nlm=use.nlm,nlm.first=nlm.first,
#                        repeat.fit=repeat.fit,max.ntry=max.ntry,
#                        parscale=parscale,save.sim=save.sim)                                                           
          # record job details
          job.lookup[[as.character(job.id)]] <- list(s=i,mnam=j)              
          # Remove fit from tostart candidates
          tostart[[i]] <- tostart[[i]][tostart[[i]]!=j]
          if ( n.running >= hydra.max ) break        
        } # all starts launched for current subject 
        if ( n.running >= hydra.max ) break
        repeat { # get next subject
          current.i <- current.i+1
          if ( current.i>length(ss) || length(tostart[[current.i]]!=0)  ) break
        }
        if ( current.i > length(ss) ) break # all done
      } # finish launch      
      # report
      cat(paste("\nRUNNING",sum(unlist(lapply(tofit,length)))-sum(unlist(lapply(tostart,length))),
       "fits,",sum(unlist(lapply(tostart,length))),"fits remaining to run\n")) 
      # collect results
      repeat {
        results <- hydraCollect()      
        # pause if no results available  
        if (length(results) == 0) 
          Sys.sleep(0.2) else break
      }
      
#       # process results
      for ( n in names(results) ) {
        n.running <- n.running - 1
        # get info about job
        job.id <- as.integer(n)
        i <- job.lookup[[as.character(job.id)]]$s
        j <- job.lookup[[as.character(job.id)]]$mnam
        # update done
        tofit[[i]] <- tofit[[i]][tofit[[i]]!=j]
        # save results  
        load(ff[[i]])
        attr(sims,"dev")[j] <- results[[n]]$fit$dev        
        attr(sims,"best")[[j]] <- results[[n]]$fit$par
        if (save.sim) for (k in 1:length(results[[n]]$sim))
          sims[[k]][[j]] <- results[[n]]$sim[[k]][[1]]
        save(sims,file=ff[[i]])
      } # results all processed 
      if ( sum(unlist(lapply(tofit,length)))==0 ) break
    } # all collected  
  } # end did fits
}


# niter=NA;ndigit=12;repeat.fit=.1;max.ntry=20;use.nlm=T;use.simplex=T;nlm.first=T;parscale=T; max.ntry=20
time.fit <- function(strt,dat,pmap,niter=NA,ndigit=12,repeat.fit=.1,max.ntry=20,
  use.nlm=T,use.simplex=T,nlm.first=T,parscale=T,precision=2.5) {

  repeat.fit <- sort(repeat.fit,decreasing=T)
  out <- matrix(ncol=4,nrow=length(repeat.fit)*(use.nlm+use.simplex+nlm.first))
  dimnames(out) <- list(NULL,c("type","tol","t","Objective"))
  fit <- vector(mode="list",length=0)
  fit$par <- strt
  rown <- 0
  if ( use.nlm ) {
    fit.nlm=fit
    for (i in 1:length(repeat.fit)) {
      nlm.t <- system.time({
        fit.nlm <- multi.start.fit(fit.nlm$par,pmap,dat,max.ntry=max.ntry,
          niter=niter,ndigit=ndigit,repeat.fit=repeat.fit[i],
          use.nlm=T,nlm.first=T,parscale=parscale,precision=precision)
      })[3]
      if (i>1) nlm.t <- nlm.t + as.numeric(out[rown,"t"])
      rown <- rown+1
      out[rown,] <- c("nlm",repeat.fit[i],nlm.t,fit.nlm$value) 
    }
    bestfit <- fit.nlm
  }
  if ( use.simplex ) {
    fit.simplex=fit
    for (i in 1:length(repeat.fit)) {
      simplex.t <- system.time({
        fit.simplex <- multi.start.fit(fit.simplex$par,pmap,dat,max.ntry=max.ntry,
          niter=niter,ndigit=ndigit,repeat.fit=repeat.fit[i],
          use.nlm=F,nlm.first=F,parscale=parscale,precision=precision)
      })[3]
      if (i>1) simplex.t <- simplex.t + as.numeric(out[rown,"t"])
      rown <- rown+1
      out[rown,] <- 
        c("simplex",repeat.fit[i],simplex.t,fit.simplex$value) 
    }
    if (!use.nlm) bestfit <- fit.simplex else 
      if (fit.simplex$value < bestfit$value) bestfit <- fit.simplex
  }
  if ( nlm.first ) {
    fit.ns=fit
    for (i in 1:length(repeat.fit)) {
      ns.t <- system.time({
        fit.ns <- multi.start.fit(fit.ns$par,pmap,dat,max.ntry=max.ntry,
          niter=niter,ndigit=ndigit,repeat.fit=repeat.fit[i],
          use.nlm=F,nlm.first=nlm.first,parscale=parscale,precision=precision)
      })[3]
      nlm.first=F
      if (i>1) ns.t <- ns.t + as.numeric(out[rown,"t"])
      rown <- rown+1
      out[rown,] <- 
        c("nlm1",repeat.fit[i],ns.t,fit.ns$value) 
    }
    if (!use.nlm & !use.simplex) bestfit <- fit.ns else
      if (fit.ns$value < bestfit$value) bestfit <- fit.ns
  }
  attr(bestfit,"times") <- out
  bestfit
}

# subject=1;model=NA;start=NA
# niter=NA;repeat.fit=.1;ndigit=10;use.nlm=F;use.simplex=T;nlm.first=F
# model = 1; subject="21"
fit.one <- function(dname,subject=1,model=NA,start=NA,parscale=T,
    niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F) {
  
  fitlist=get.fitlist(dname)
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  mnames <- names(fitlist[[hmname]]$MT[[1]])
  if (is.na(model)) model <- mnames[length(mnames)]
  if (is.numeric(model)) model <- mnames[model]
  if (!(model %in% mnames))
    stop("Model not in model names")
  if (any(is.na(subject)) | length(subject)>1)
    stop("Specify a single subject (name or number)")
  if (is.numeric(subject)) 
    subject <- subject.names(fitlist)[subject]
  dat=get.dat(fitlist,subject)
  pmap=get.pmap(fitlist,dname,model,subject)
  if (any(is.na(start))) strt <- start.pmap(pmap) else {
    if (length(start)>1) strt <- start else {
      if (is.numeric(start.model)) start <- mnames[start]
      if (!(start %in% mnames))
        stop("Start model not in model names")
      fit <- vector(mode="list",length=2)
      names(fit) <- c("pmap","fit")
      spmap <- get.pmap(fitlist,dname,start,subject)
      strt <- best.pmap(spmap)
      if (is.null(strt)) 
        stop("Start model has not had a best fit recorded")
      fit$pmap <- spmap
      fit$fit$par <- strt
      strt <- getstart(new.fit,pmap) 
    }
  }
  pu <- as.numeric(dimnames(attr(dat,"D")$nc)$minmax)
  qmpe <- !is.null(attr(dat,"qp"))
  strtobj <- objective(strt,dat,pmap,pu,qmpe)  
  time.fit(strt,dat,pmap,niter=niter,ndigit=ndigit,
    repeat.fit=repeat.fit,parscale=parscale,
    use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)  
}

fix.fit <- function(dname,model=1,start.model=2,subject=1,parscale=T,
  niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F,
  save.fit=T,show.par=F,nesttol=.1,topmodel=1) {

  # Main Body
  fitlist=get.fitlist(dname)
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  mnames <- names(fitlist[[hmname]]$MT[[1]])
  if (any(is.na(subject)) | length(subject)>1)
    stop("Specify a single subject (name or number)")
  if (is.numeric(subject)) 
    subject <- subject.names(fitlist)[subject]
  if (is.numeric(model)) 
    model <- model.names(fitlist)[model]
  if (!(model %in% mnames))
    stop("Model not in model names")
  if (is.numeric(start.model)) 
    start.model <- model.names(fitlist)[start.model]
  if (!(start.model %in% mnames))
    stop("Start model not in model names")
  dat=get.dat(fitlist,subject)
  pmap=get.pmap(fitlist,dname,model,subject)
  # Start model fit
  fitname.start <- paste(dname,"/",subject,"/fits/",start.model,".RData",sep="")
  if (file.exists(fitname.start)) load(fitname.start) else 
    stop("No old fit to fix")
  fit.start <- fit
  cat(paste("     Start model deviance:",fit.start$dev,"\n")) 
  # Old fit
  fitname <- paste(dname,"/",subject,"/fits/",model,".RData",sep="")
  if (file.exists(fitname)) load(fitname) else 
    stop("No old fit to fix")
  cat(paste("    Model to fix deviance:",fit$dev,"\n"))
  # only fix if some previous fix has not already dealt with 
  if (fit$dev > fit.start$dev) { 
    fit.old <- fit
    # Get start point model parameters and map
    spmap <- get.pmap(fitlist,dname,start.model,subject)
    strt <- best.pmap(spmap)
    if (is.null(strt)) 
      stop("Start model has not had a best fit recorded")
    # Get new fit start point from start point model
    new.fit <- vector(mode="list",length=2)
    names(new.fit) <- c("pmap","fit")
    new.fit$pmap <- spmap
    new.fit$fit$par <- strt
    strt <- getstart(new.fit,pmap)
    # Get deviance for new start point
    pu <- as.numeric(dimnames(attr(dat,"D")$nc)$minmax)
    qmpe <- !is.null(attr(dat,"qp"))
    strtobj <- 2*objective(strt,dat,pmap,pu,qmpe)  
    if (qmpe) { # zero for prefect fit
      strtobj <- strtobj + 2*sum(dat$qn*log(unlist(
        tapply(dat$qp,dat$cell,function(x){
        x[length(x)]*diff(c(0,x[-length(x)],1))  
      }))))
    }
    cat(paste("Deviance from start model:",strtobj,"\n"))   
    # Check if stored and calculated deviance correspond
    err <- strtobj-fit.start$dev
    if (err > nesttol) {
      cat("\nThere appears to be an error in the start model's deviance, refitting ...\n")    
      fit <- time.fit(fit.start$par,dat,spmap,parscale=parscale,
        niter=niter,ndigit=ndigit,repeat.fit=repeat.fit,
        use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)
      no.better <- fit.start$dev<=fit$dev
      if (!no.better) {
        # Save refit
        attr(spmap,"best") <- fit$par
        save_pmap(spmap,start.model,subject,fitlist,hmname)
        save(fit,file=fitname.start) 
        cat(paste("     Refit start model deviance:",fit$dev,"\n")) 
        # Get new fit start point from start point model 
        new.fit$fit$par <- fit$par
        strt <- getstart(new.fit,pmap)
        # Get deviance for new start point
        old.strtobj <- strtobj
        strtobj <- 2*objective(strt,dat,pmap,pu,qmpe)  
        if (qmpe) { # zero for prefect fit
          strtobj <- strtobj + 2*sum(dat$qn*log(unlist(
            tapply(dat$qp,dat$cell,function(x){
            x[length(x)]*diff(c(0,x[-length(x)],1))  
          }))))
        }    
        cat(paste("Deviance from refit start model:",strtobj,"\n"))
      } else cat(paste("Deviance from refit start model no better.\n"))
    }
    cat(paste("Trying 10 times harder on a refit.\n"))
    new.fit <- time.fit(strt,dat,pmap,niter=niter,ndigit=ndigit,
      repeat.fit=repeat.fit/10,parscale=parscale,
      use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)
    cat(paste("     Deviance for new fit:",new.fit$dev,"\n")) 
    if (show.par)
      print(rbind(Old.Fit=fit.old$par,New.Fit=new.fit$par))
    lost.cause <- new.fit$value >= fit$value
    if (lost.cause) cat(paste("Well that didn't help.\n"))
    if (save.fit & !lost.cause) { # save results if better
      cat(paste("Problem solved!!! Saving improved fit.\n"))
      fit <- new.fit
      attr(pmap,"best") <- fit$par
      save_pmap(pmap,model,subject,fitlist,hmname)
      save(fit,file=fitname) 
    }
  } else lost.cause <- F
#   out <- get.fits(dname,nesttol=nesttol,get.bma=F,topmodel=topmodel)
#   names(lost.cause) <- model
#   attr(out,lost.cause=list.cause)
#   out
  lost.cause
}


# nesttol=.1;until.done=T;topmodel=1
# niter=NA;repeat.fit=.1;ndigit=10;use.nlm=F;use.simplex=T;nlm.first=F
fix.fits <- function(dname,nestfail,nesttol=.1,until.done=T,topmodel=1,
  niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F) {
  nfail <- dim(nestfail)[1]
  repeat {
    if (nfail>0) {
      lost.causes <- logical(nfail)
      for (i in 1:nfail) {
        cat("Fixing subject",nestfail[i,"s"],"model",
          as.character(nestfail[i,"Nesting"]),"\n")
        lost.causes[i] <- fix.fit(dname,subject=as.character(nestfail$s)[i],
          model=as.character(nestfail$Nesting)[i],
          start.model=as.character(nestfail$Nested)[i],
          niter=niter,repeat.fit=repeat.fit,ndigit=ndigit,
          use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first,
          topmodel=topmodel)  
      }
      tmp <- get.fits(dname,nesttol=nesttol,get.bma=F,topmodel=topmodel)
      if (any(lost.causes))
        cat(paste("Giving up on ",sum(lost.causes),"lost causes.\n"))
      tmp <- tmp[!lost.causes,] 
      if (dim(tmp)[1]!=0) cat("More fits to fix, run fix.fits again\n")
    } else cat("No fits to fix\n")
    if (!until.done || (nfail==0) || dim(tmp)[1]==0) break
    nestfail <- tmp
    nfail <- dim(nestfail)[1]
    cat(paste("\n\n",nfail,"fits remain to be fixed\n\n"))
    print(nestfail)
  }  
}

# fit.old.value=fit.old$value; sstrt=fit.start$par
fix.fit.hydra <- function(dat,pmap,strt,spmap,sstrt,fit.old.value,
  fitname,fitname.start,subject,start.model,model,fit.start.value,parscale=T,max.ntry=20,
  niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F,
  nesttol=.1,precision=2.5) {
  
  if ( !is.null(spmap) ) { # refit start
    fit.start <- time.fit(sstrt,dat,spmap,parscale=parscale,precision=precision,
        niter=niter,ndigit=ndigit,repeat.fit=repeat.fit,max.ntry=max.ntry,
        use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)
    if (fit.start.value < fit.start$value) {
      attr(spmap,"best") <- fit.start$par
      new.fit <- vector(mode="list",length=2)
      names(new.fit) <- c("pmap","fit")
      new.fit$pmap <- spmap
      new.fit$fit$par <- fit.start$par
      strt <- getstart(new.fit,pmap)
    } else {
      fit.start <- NULL
      spmap <- NULL
    }
  } else fit.start <- NULL
  new.fit <- time.fit(strt,dat,pmap,niter=niter,ndigit=ndigit,parscale=parscale,
    max.ntry=max.ntry,precision=precision,
    repeat.fit=repeat.fit/10,use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)
  lost.cause <- new.fit$value >= fit.old.value
  if (!lost.cause) attr(pmap,"best") <- new.fit$par else {
    new.fit <- NULL
    pmap <- NULL
  }
  list(new.fit=new.fit,pmap=pmap,fit.start=fit.start,spmap=spmap,
       fitname=fitname,fitname.start=fitname.start,subject=subject,
       start.model=start.model,model=model)
}


# dname="g1q.LBA"; model="B~1 & A~1 & v~b*C & sv~1 & ter~1 & pc~1"
# start.model="B~1 & A~1 & v~C & sv~1 & ter~1 & pc~1"; repeat.fit=c(10)
# subjects=NA; niter=NA; ndigit=10; use.nlm=F; use.simplex=T; nlm.first=F
# est.time=60*60*10;save.fits=T
hydra.fits <- function(dname,model=1,start.model=NA,subjects=NA,
  niter=NA,repeat.fit=.1,ndigit=10,use.nlm=F,use.simplex=T,nlm.first=F,
  est.time=60*60*10,save.fits=T,parscale=T) {

  # Main Body
  fitlist=get.fitlist(dname)
  hmname <- names(fitlist)[2]
  attr(fitlist[[hmname]]$MT,"dname") <- dname
  if (any(is.na(subjects))) subjects <- subject.names(fitlist) else
    if (is.numeric(subjects)) subjects <- subject.names(fitlist)[subjects]
  if (is.numeric(model)) model <- model.names(fitlist)[model[1]]
  cat(paste("Timing",length(subjects),"subjects\n"))
  fits <- vector(mode="list",length=length(subjects))
  names(fits) <- subjects
  job.ids <- numeric(length(subjects))
  names(job.ids) <- subjects
  for (s in subjects) { # start all subjects
    cat(paste(s,""))
    dat=get.dat(fitlist,s)
    pmap=get.pmap(fitlist,dname,model,s)
    if (any(is.na(start.model))) 
      strt <- start.pmap(pmap) else {
      for (i in start.model) {
        spmap <- get.pmap(fitlist,dname,i,s)
        strti <- best.pmap(spmap)   
        if (is.null(strti)) 
          stop("Start model has not had a best fit recorded") else {
            fit <- vector(mode="list",length=2)
            names(fit) <- c("pmap","fit")
            fit$pmap <- spmap
            fit$fit$par <- strti
            if (i==start.model[1]) 
              strt <- getstart(fit,pmap) else
              strt <- cbind(strt,getstart(fit,pmap))
          }
      }
    }
    job.ids[s] <- hydraRun(time.fit,strt=strt,dat=dat,pmap=pmap,ndigit=ndigit,
        repeat.fit=repeat.fit,niter=niter,est.time=est.time,parscale=parscale,
        use.nlm=use.nlm,use.simplex=use.simplex,nlm.first=nlm.first)
  }
  times <- matrix(ncol=5,nrow=0)
  dimnames(times) <- list(NULL,c("s","type","tol","t","Objective"))
  repeat {
    repeat { # collect results
      results <- hydraCollect()      
      if (length(results) == 0) Sys.sleep(0.2) else break
    }
    for (n in names(results)) {
      job.id <- as.integer(n)
      s <- names(job.ids)[job.ids==job.id]
      fits[[s]] <- results[[n]]
      times <- rbind(times,cbind(s=rep(s,dim(attr(results[[n]],"times"))[1]),
        attr(results[[n]],"times")))
      job.ids[s] <- 0
      cat(paste("Collected participant",s,
        "(",sum(job.ids!=0),"left to go )\n"))
      if (save.fits) { # save results if better
        fitname <- paste(dname,"/",s,"/fits/",model,".RData",sep="")
        if (file.exists(fitname)) load(fitname) else 
          fit <- list(value=Inf,start.from="NA")
        if (any(names(results[[n]])=="value") && results[[n]]$value < fit$value) { 
          fit <- results[[n]]
          pmap <- get.pmap(fitlist,model=model,subject=s,hmodel=hmname) 
          start <- fit$par
          names(start) <- names(attr(pmap,"start"))
          attr(pmap,"best") <- start
          save_pmap(pmap,model,s,fitlist,hmname)
          save(fit,file=fitname) 
        }
      }
    }
    if (all(job.ids==0)) break
  }
  cat("\n")
  out=data.frame(times,row.names=NULL)
  out$t <- as.numeric(as.character(out$t))
  out$Objective <- as.numeric(as.character(out$Objective))
  cat("\nAverage time (sec)\n")
  print(tapply(out$t,list(type=out$type,tol=out$tol),mean))
  cat("\nAverage objective\n")
  print(tapply(out$Objective,list(type=out$type,tol=out$tol),mean))
  attr(fits,"times") <- out
  invisible(fits)
}

