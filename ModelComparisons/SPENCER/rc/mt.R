####################  MODEL HIERARCHY GENERATION   #############################
#
# FIRST LATENT VARIABLE SPECIFIED (IF ANY) MUST CORRESPOND TO RESPONSE
# 


# Latents=NULL;Cpars=NA;SRpars=NULL;rawdat=NA
# stopFormulas=NULL;Contrasts=NULL;Constants=NULL;Bounds=NULL;Constraints=NULL
# showmodels=FALSE;dname="MT";saturated=TRUE;update=TRUE 
# 
# dname="hayes";hmodelname="LBAe";fitlist=vs
#          Latents=list(lR=c("present","absent"));Cpars=c("aV","bV","cV","sv");showmodels=F
#          Formulas=list(
#            aa="~1",ba="~1",ca="~1",
#            ab="~lR",bb="~1",cb="~1",
#            aV="~C",bV="~PA*DS*C",cV="~1",
#            sv="~C",
#            aT="~1",bT="~1",cT="~1",
#            pc="~1")
#          Bounds=list(
#            aa=c(0,1),ba=c(0,1),ca=c(0,Inf),
#            ab=c(0,Inf),bb=c(0,Inf),cb=c(0,Inf),
#            aV=c(-Inf,Inf),bV=c(0,Inf),cV=c(0,Inf),
#            sv=c(0,Inf),
#            aT=c(0,Inf),bT=c(0,Inf),cT=c(0,Inf),
#            pc=c(0,1))
#          stopFormulas=list(aV="~C")
#          Constants=list(
#            sv=c(I=0),
#            pc=c(I=qlogis(.00001)))
#          saturated=T

# dname="baker";hmodelname="LBA";fitlist=baker;
#          Latents=list(lR=c("lure","target"));Cpars=c("v","sv");showmodels=F;
#          Formulas=list(B="~lR*LH",A="~1",v="~S*F*D*LH*C",sv="~C",ter="~1",pc="~1");
#          stopFormulas=list(v="~C");
#          Bounds=list(B=c(0,Inf),A=c(0,Inf),v=c(-Inf,Inf),sv=c(0,Inf),ter=c(0,Inf),pc=c(0,1));
#          Constants=list(sv=c(I=0),pc=c(I=qlogis(.00001)));
#          saturated=T

# dname="stroop1vsvA";hmodelname="LBA";fitlist=stroop1
#          Latents=list(lR=c("red","green","blue","yellow"))
#          SRpars=c("v","sv2");showmodels=FALSE
#          Formulas=list(B="~lR*IS",A="~1",v="~RS*SR",sv2="~SR",ter="~1",pc="~1")
#          Bounds=list(B=c(0,Inf),A=c(0,Inf),v=c(-Inf,Inf),sv2=c(0,Inf),ter=c(0,Inf),pc=c(0,1))
#          stopFormulas=list(v="~SR")
#          Constants=list(sv2=c(I=0),v=c(I=0),pc=c(I=qlogis(.00001)))
#          Contrasts=list(v=list(SR=contr),sv2=list(SR=contr))
#          saturated=T

# dname="stroop2vsvAa";hmodelname="LBA";fitlist=stroop2a
#          Latents=list(lR=c("MC.1","MC.2","MC.3","MC.4","MI.1","MI.2","MI.3","MI.4"))
#          SRpars=c("v","sv2");showmodels=FALSE
#          Formulas=list(B="~lR",A="~1",v=c("~SR","~MSR"),sv2="~1",ter="~1",pc="~1")
#          Bounds=list(B=c(0,Inf),A=c(0,Inf),v=c(-Inf,Inf),sv2=c(0,Inf),ter=c(0,Inf),pc=c(0,1))
#          stopFormulas=list(v="~SR")
#          Constants=list(sv2=c(I=0),v=c(I=0),pc=c(I=qlogis(.00001)))
#          Contrasts=list(v=list(SR=contr,MSR=contr2),
#                      sv2=list(SR=contr,MSR=contr2))
#          saturated=T
 
model.rc=function(hmodelname,fitlist,Formulas,Latents=NULL,Cpars=NA,SRpars=NULL,rawdat=NA,
  stopFormulas=NULL,Contrasts=NULL,Constants=NULL,Bounds=NULL,Constraints=NULL,
  showmodels=FALSE,dname="MT",saturated=TRUE,update=TRUE) {
  # Makes M (model matrix, D augmented with latents) and MT (model tree,
  # includig start points), adds fitlist$hmodelname with M and MT (names
  # only), saves saves it and seperate model tree file for each subject (
  # named with subject names) to MT

  hforms.old <- function(f,f0="~1",mapnumbers=T,saturated=F) {
  # f is formula, f0 is stop formula
  # returns flist, a list of formulas with attribute "maps" specifying 
  # nesting relationships between formulas
  # if mapnumbers=F maps specified by names, if T by numbers (required
  # for use by makehmaps)

    # returns true if character vector x in nestlist is hierarchical     
    is.h <- function(x,nlist){
      okj <- T
      for ( j in x[length(x):1] )  {
        if ( !is.null(nlist[[j]]) ) 
          if ( !all(nlist[[j]] %in% x[x!=j]) ) {
            okj <- F
            break
          }      	
      }
      okj	
    }
    
    if ( class(f)!="formula" ) {
      if (substr(f,1,1) != "~") {
        noform = T
      } else {
        f <- formula(f)
        if (class(f0)!="formula") f0 <- formula(f0)
        noform <- F
      }
    } else noform <- F
    if ( noform ) {
      facnams <- strsplit(f,",")[[1]]
      if (f0=="~1") facnams <- c("1",facnams) # otherwise first element of f
      flist <- vector(mode="list",length=length(facnams))
      names(flist) <- facnams[length(facnams):1]
      for (i in facnams) flist[[i]] <- formula(paste("~",i,sep=""))
      attr(flist[[length(facnams)]],"maps") <- logical(0)
      for (i in 1:(length(facnams)-1)) 
        attr(flist[[i]],"maps") <- names(flist)[i+1]
    } else {
      if ( f==formula("~1") || f==f0 ) {       # only one formula
        flist <- list('1'=f)
        attr(flist[[1]],"maps") <- logical(0)
      } else {                                 # create full list of formulas
        # expand to full list of terms strings, e.g. a*b -> a b a:b
        trms <- attr(terms(f),"term.labels")
        is.intercept <- as.logical(attr(terms(f),"intercept")) 
        ntrms <- length(trms)
        
        # create list to check nesting of interaciton terms
        nestlist <- vector(mode="list",ntrms)
        names(nestlist) <- trms[ntrms:1]
        for ( i in names(nestlist) ) {
          facs <- strsplit(i,":",fixed=T)[[1]]
          nfacs <- length(facs)        
          if ( nfacs>1 ) for ( j in (nfacs-1):1 )   
            nestlist[[i]] <- c(nestlist[[i]],
                               apply(combn(nfacs,j),2,function(x){paste(facs[x],collapse=":")}))
        }      
        
        # create formulas list
        fvec <- character(0)
        if (ntrms==1) fvec <- trms else for ( i in 1:ntrms ) {
          # generate all combinations of the elements of trms taken i at a time
          fmat <- matrix(apply(combn(ntrms,i),2,function(x){trms[x]}),nrow=i)	    
          # keep hierarchical columns
          fvec <- c(fvec,apply(matrix(fmat[,apply(fmat,2,is.h,nlist=nestlist)],
                                      nrow=dim(fmat)[1]),2,paste,collapse="+"))
        }
        if (saturated) fvec <- fvec[sapply(fvec,is.saturated)]
        
        flist <- sapply(c("1",fvec)[(length(fvec)+1):1],
                        function(x){formula(paste("~",x,sep=""))},simplify=F)
        labs <- lapply(flist,function(x){attr(terms(x),"term.labels")})
        lens <- unlist(lapply(labs,length)) # number of terms
        
        # make a list to contain maps for each model
        maps <- lapply(vector(length=length(flist),mode="list"),
                       function(x){vector(length=0)})
        
        # for all but bottom formula get formulas that are nested
        for ( i in 1:(length(flist)-1) ) {
          lensi <- lens[i] # number of terms in i'th formula
          j0 <- i+1        # numer of next lowest order model
          # repeat until formula with less terms as only want immediate nests
          repeat { 
            if ( lensi>lens[j0] ) break
            j0 <- j0+1
          }
          lensj <- lens[j0] # number of terms in first immediate nests
          # add nested number for first immediate and following immediate nests
          for ( j in (j0):length(flist) ) {
            if ( lensj>lens[j] ) break # stop as down two steps
            # add in nested model number if all terms in model above
            if ( all(sapply(labs[j][[1]],function(x){any(x==labs[i][[1]])})) ) 
              maps[[i]] <- c(maps[[i]],names(flist)[j])
          }  
        }
        
        # put nested formula numbers into attribute
        for ( i in 1:length(flist) ) attr(flist[[i]],"maps") <- maps[[i]]
        
        # make stopped flist  
        if ( !(is.null(f0) || f0==formula("~1")) ) {     
          # find f0 in flist names
          nlist <- strsplit(names(flist),"+",fixed=T)
          f0n <- attr(terms(f0),"term.labels")
          f0n <- names(flist)[unlist(lapply(nlist,function(x){
            all(f0n %in% x) & length(f0n)==length(x)}))]
          if (length(f0n)!=1) 
            stop("Could not find formula matching stop formula")
          f0 <- formula(paste("~",f0n))
          # construct reverse maps
          usedby <- vector(mode="list",length=length(flist)-1)
          names(usedby) <- names(flist)[-length(flist)]
          for ( i in names(usedby)[length(usedby):1] )
            usedby[[i]] <- names(flist)[unlist(lapply(flist,function(x){
              any(attr(x,"maps")==i)}))]
          # construct stopped flist 
          slist <- list()
          snams <- paste(attr(terms(f0),"term.labels"),collapse="+")
          slist[[snams]] <- flist[[snams]]
          attr(slist[[snams]],"maps") <- logical(0)
          repeat {
            for (i in snams) {
              usevec=usedby[[i]]
              for (j in usevec) {
                slist[[j]] <- flist[[j]]
                attr(slist[[j]],"maps") <- attr(slist[[j]],"maps")[
                  attr(slist[[j]],"maps") %in% names(slist)]	        
              }
            }
            snams <- usevec
            if ( any(snams==names(flist)[1]) ) break           
          }
          flist <- slist
          # reverse order of flist
          rnams <- names(slist)[length(slist):1]
          for (i in 1:length(flist)) flist[[i]] <- slist[[rnams[i]]]
          names(flist) <- rnams
        }
        if (!is.intercept) {
          flist <- lapply(flist,function(x){
            out <- formula(paste(c(as.character(x),"+0"),collapse=""))
            maps <- attr(x,"maps")
            maps <- sapply(maps,function(x){paste(x,"+0",sep="")})
            names(maps) <- NULL
            attr(out,"maps") <- maps
            out
          })
          attr(flist[[length(flist)]],"maps") <- logical(0)
          names(flist) <- lapply(flist,function(x){
            paste(strsplit(as.character(x)[-1]," ")[[1]],collapse="")})
        }
      }      
    }
    if ( mapnumbers ) {
      mnams <- names(flist)
      mnums <- 1:length(mnams)
      names(mnums) <- mnams
      for ( i in mnams[-length(mnams)] ) 
        attr(flist[[i]],"maps") <- mnums[attr(flist[[i]],"maps")]
    }
    flist
  }


  #design=fitlist[[hmodelname]]$M[[1]]$D;Formulas=fi
  makemap <- function(design,Formulas,
    Contrasts=NULL,Constants=NULL,Bounds=NULL,Constraints=NULL) {
    # called by addhmodel, to make a single entry in the list of hierarchical
    # models made by makehmaps or can be called directly

    makedf <- function(facs) {
      fn <- unlist(lapply(facs,length))
      nr <- prod(fn)
      fnams <- names(facs)
      k <- 1
      dd <- data.frame(gl(n=fn[1],k=k,length=nr,labels=facs[[fnams[1]]]))
      k <- k*fn[1]
      for (i in fnams[-1] ) {
        dd <- cbind(dd,gl(n=fn[i],k=k,length=nr,labels=facs[[i]]))
        k <- k*fn[i]
      }
      names(dd) <- names(facs)
      dd
    }
   
    get.constraints <- function(pmap,Constaints) {
      nrow <- dim(pmap[[1]]$mm)[1]
      ci <- vector(length=0)
      for (i in names(Constraints)) {  
        if (length(Constraints[[i]])>1)
          stop(paste("Constraint for parameter",i,"is not a scalar"))
        if ( !is.na(Constraints[[i]]) ) 
    	    ci <- c(ci,rep(Constraints[[i]],nrow)) else
    	    ci <- c(ci,rep(-Inf,nrow))
      }
      ci
    }
      
    # design = data frame of factors, one row per design cell starting  
    # with special factors C (correct choice) and R (response = accumulator)
    # in columns 1 and 2 (so can have arbitary names)
    pnams <- names(Formulas)
    out <- vector(length=length(pnams),mode="list")
    names(out) <- pnams
    first <- 1 # position in fitting vector
    mapmat <- matrix(nrow=dim(design)[1],ncol=0)
    for ( i in  pnams ) {
      Formulas[[i]] <- formula(Formulas[[i]]) # in case not already a formula
      # get factor names
      ffnams <- dimnames(attr(terms(Formulas[[i]]),"factors"))[[1]]     
      # Some checking
      for ( j in ffnams ) if ( !any(j==names(design)) )
        stop(paste("Parameter \"",i,"\" formula has factor name \"",j,
          "\" not in data",sep=""))    
      # Make model matrix and put it in output list with formula
      # supress warning when when factors in Contrasts not in Formulas 
      if (all(is.na(Contrasts[[i]]))) 
        contrasts <- NULL else
        contrasts <- Contrasts[[i]]
      suppressWarnings({out[[i]] <- list(formula=Formulas[[i]],
        mm=model.matrix(Formulas[[i]],data=design,contrasts.arg=contrasts) )})      
      if ( attr(terms(Formulas[[i]]),"intercept")==0 ) {
        mm <- out[[i]]$mm
        con <- attr(mm,"contrasts")[[1]]
        ncol=dim(con)[1]
        newmm <- mm[,2:ncol,drop=F]
        for (j in 2:ncol) newmm[,j-1] <- con[,(j-1)]*apply(mm[,1:j],1,sum)
        mm[,2:ncol] <- newmm
        mm <- mm[,-1]
        tmp <- attributes(out[[i]]$mm)
        tmp$dim[2] <- tmp$dim[2]-1
        tmp$dimnames[[2]] <- tmp$dimnames[[2]][-1]
        tmp$assign <- tmp$assign[-1]
        attributes(mm) <- tmp
        out[[i]]$mm <- mm 	
      } else dimnames(out[[i]]$mm)[[2]][1] <- "I" # "(Intercept)" by detault
      
      # update fitting vector
      
      if ( is.null(names(Constants[[i]])) ) {  # No constants
        out[[i]]$n <- dim(out[[i]]$mm)[2]   # n = number fit parameters
        last <- first+out[[i]]$n-1          # position of last element 
        out[[i]]$pos <- first:last          # pos = position range
        first <- last+1                     # start of next part of fit vector
      } else {                          # Constant parameter value
        if ( any(!is.numeric(Constants[[i]])) )
          stop(paste("Constant",i,"must be numeric"))
        dn <- dimnames(out[[i]]$mm)[[2]]
        consts <- rep(NA,length(dn))
        names(consts) <- dn
        in.mm <- dn %in% names(Constants[[i]])
#         consts[in.mm] <- Constants[[i]][in.mm]
        in.Constants <- names(Constants[[i]]) %in% dn
        consts[in.mm] <- Constants[[i]][in.Constants]

        out[[i]]$consts <- consts # Constant values
        out[[i]]$n <- sum(is.na(consts))      
        if (out[[i]]$n>0) {
          last <- first+out[[i]]$n-1          # position of last element 
          out[[i]]$pos <- first:last          # pos = position range
          first <- last+1
        }
      }
       
      # Bounds checks
      if ( length(Bounds[[i]])!=2) 
        stop(paste("Bounds for parameter \"",i,"\" must be length=2",sep=""))
      if ( !all(is.numeric(Bounds[[i]])) ) 
        stop(paste("Bounds for parameter \"",i,"\" must be numeric",sep=""))
      lbound <- Bounds[[i]][1]
      ubound <- Bounds[[i]][2]
      if ( ubound<=lbound ) 
        stop(paste("Parameter \"",i,"\" lower >= upper bound",sep=""))
      # Get bound type and store values (double a bottom and range)
      if ( lbound==-Inf & ubound==Inf )         # no bound
        out[[i]]$partype <- "unbounded" else
      if ( is.finite(lbound) & ubound==Inf ) { # postive by log
        out[[i]]$partype <- "positive" 
        out[[i]]$bound <- Bounds[[i]]
      } else {
        out[[i]]$partype <-"doublebound"   # double bound by logistic
        out[[i]]$bound <- c(Bounds[[i]][1],diff(Bounds[[i]]))
      }
    }
    attr(out,"npar") <- last # length of fitting parameter vector
    if (!is.null(Constraints)) 
      attr(out,"ci") <- get.constraints(out,Constraints)  
    for ( i in  pnams ) # otherwise expands save massively
      attr(out[[i]]$formula,".Environment") <- NULL 
    out
  }

#   f=formula(Formulas[[i]]); f0=formula(stopFormulas[[i]]); mapnumbers=T
  hforms <- function(f,f0="~1",mapnumbers=T,saturated=F) {
  # f is formula, f0 is stop formula
  # returns flist, a list of formulas with attribute "maps" specifying 
  # nesting relationships between formulas
  # if mapnumbers=F maps specified by names, if T by numbers (required
  # for use by makehmaps)

    # Used to expand saturated formulas
    starform <- function(starstr) {
      formula(paste("~",paste(attr(terms(formula(paste("~",starstr,sep="")
        )),"term.labels"),collapse="+"),sep=""))
    }

    # returns true if character vector x in nestlist is hierarchical     
    is.h <- function(x,nlist){
      okj <- T
      for ( j in x[length(x):1] )	{
        if ( !is.null(nlist[[j]]) ) 
          if ( !all(nlist[[j]] %in% x[x!=j]) ) {
            okj <- F
            break
          }      	
      }
      okj	
    }
    
    if ( class(f)!="formula" ) {
      if (substr(f,1,1) != "~") {
        noform = T
      } else {
        f <- formula(f)
        if (class(f0)!="formula") f0 <- formula(f0)
        noform <- F
      }
    } else noform <- F
    if ( noform ) {
      facnams <- strsplit(f,",")[[1]]
      if (f0=="~1") facnams <- c("1",facnams) # otherwise first element of f
      flist <- vector(mode="list",length=length(facnams))
      names(flist) <- facnams[length(facnams):1]
      for (i in facnams) flist[[i]] <- formula(paste("~",i,sep=""))
      attr(flist[[length(facnams)]],"maps") <- logical(0)
      for (i in 1:(length(facnams)-1)) 
        attr(flist[[i]],"maps") <- names(flist)[i+1]
    } else {
      if ( f==formula("~1") || f==f0 ) {       # only one formula
        flist <- list('1'=f)
        attr(flist[[1]],"maps") <- logical(0)
      } else {                                 # create full list of formulas
        # expand to full list of terms strings, e.g. a*b -> a b a:b
        trms <- attr(terms(f),"term.labels")
        is.intercept <- as.logical(attr(terms(f),"intercept")) 
        ntrms <- length(trms)
        
        # create list to check nesting of interaciton terms
        nestlist <- vector(mode="list",ntrms)
        names(nestlist) <- trms[ntrms:1]
        for ( i in names(nestlist) ) {
          facs <- strsplit(i,":",fixed=T)[[1]]
          nfacs <- length(facs)        
          if ( nfacs>1 ) for ( j in (nfacs-1):1 )   
            nestlist[[i]] <- c(nestlist[[i]],
                               apply(combn(nfacs,j),2,function(x){paste(facs[x],collapse=":")}))
        }      
                
        if (saturated) {
          fnams <- strsplit(as.character(f)[[2]]," * ",fixed=T)[[1]]
          fvec <- character(0)
          for (i in 1:length(fnams)) {
            fmat <- combn(fnams,i)
            fvec <- c(fvec,apply(fmat,2,paste,collapse="*"))
          }
          flist <- sapply(fvec[length(fvec):1],starform,simplify=F)
          names(flist) <- unlist(lapply(flist,function(x){
            paste(attr(terms(x),"term.labels"),collapse="+")}))
          flist <- c(flist,'1'=formula("~1"))        
        } else {
          # create formulas list
          fvec <- character(0)
          if (ntrms==1) fvec <- trms else for ( i in 1:ntrms ) {
          # generate all combinations of the elements of trms taken i at a time
            fmat <- matrix(apply(combn(ntrms,i),2,function(x){trms[x]}),nrow=i)	    
          # keep hierarchical columns
            fvec <- c(fvec,apply(matrix(fmat[,apply(fmat,2,is.h,nlist=nestlist)],
                                      nrow=dim(fmat)[1]),2,paste,collapse="+"))
          }
#            if (saturated) fvec <- fvec[sapply(fvec,is.saturated)]
          flist <- sapply(c("1",fvec)[(length(fvec)+1):1],
                        function(x){formula(paste("~",x,sep=""))},simplify=F)
        }
        labs <- lapply(flist,function(x){attr(terms(x),"term.labels")})
        lens <- unlist(lapply(labs,length)) # number of terms
        
        # make a list to contain maps for each model
        maps <- lapply(vector(length=length(flist),mode="list"),
                       function(x){vector(length=0)})
        
        # for all but bottom formula get formulas that are nested
        for ( i in 1:(length(flist)-1) ) {
          lensi <- lens[i] # number of terms in i'th formula
          j0 <- i+1        # numer of next lowest order model
          # repeat until formula with less terms as only want immediate nests
          repeat { 
            if ( lensi>lens[j0] ) break
            j0 <- j0+1
          }
          lensj <- lens[j0] # number of terms in first immediate nests
          # add nested number for first immediate and following immediate nests
          for ( j in (j0):length(flist) ) {
            if ( lensj>lens[j] ) break # stop as down two steps
            # add in nested model number if all terms in model above
            if ( all(sapply(labs[j][[1]],function(x){any(x==labs[i][[1]])})) ) 
              maps[[i]] <- c(maps[[i]],names(flist)[j])
          }  
        }
        
        # put nested formula numbers into attribute
        for ( i in 1:length(flist) ) attr(flist[[i]],"maps") <- maps[[i]]
        
        # make stopped flist  
        if ( !(is.null(f0) || f0==formula("~1")) ) {     
          # find f0 in flist names
          nlist <- strsplit(names(flist),"+",fixed=T)
          f0n <- attr(terms(f0),"term.labels")
          f0n <- names(flist)[unlist(lapply(nlist,function(x){
            all(f0n %in% x) & length(f0n)==length(x)}))]
          if (length(f0n)!=1) 
            stop("Could not find formula matching stop formula")
          f0 <- formula(paste("~",f0n))
          # construct reverse maps
          usedby <- vector(mode="list",length=length(flist)-1)
          names(usedby) <- names(flist)[-length(flist)]
          for ( i in names(usedby)[length(usedby):1] )
            usedby[[i]] <- names(flist)[unlist(lapply(flist,function(x){
              any(attr(x,"maps")==i)}))]
          # construct stopped flist 
          slist <- list()
          snams <- paste(attr(terms(f0),"term.labels"),collapse="+")
          slist[[snams]] <- flist[[snams]]
          attr(slist[[snams]],"maps") <- logical(0)
          usevec <- character(0)
          repeat {
            for (i in snams) {
              usevec=c(usevec,usedby[[i]])
              for (j in usevec) {
                slist[[j]] <- flist[[j]]
                attr(slist[[j]],"maps") <- attr(slist[[j]],"maps")[
                  attr(slist[[j]],"maps") %in% names(slist)]	        
              }
            }
            snams <- unique(usevec)
            if ( any(snams==names(flist)[1]) ) break           
          }
          flist <- slist
          # reverse order of flist
          rnams <- names(slist)[length(slist):1]
          for (i in 1:length(flist)) flist[[i]] <- slist[[rnams[i]]]
          names(flist) <- rnams
        }
        if (!is.intercept) {
          flist <- lapply(flist,function(x){
            out <- formula(paste(c(as.character(x),"+0"),collapse=""))
            maps <- attr(x,"maps")
            maps <- sapply(maps,function(x){paste(x,"+0",sep="")})
            names(maps) <- NULL
            attr(out,"maps") <- maps
            out
          })
          attr(flist[[length(flist)]],"maps") <- logical(0)
          names(flist) <- lapply(flist,function(x){
            paste(strsplit(as.character(x)[-1]," ")[[1]],collapse="")})
        }
      }      
    }
    
    if ( mapnumbers ) {
      mnams <- names(flist)
      mnums <- 1:length(mnams)
      names(mnums) <- mnams
      for ( i in mnams[-length(mnams)] ) 
        attr(flist[[i]],"maps") <- mnums[attr(flist[[i]],"maps")]
    }
    flist
  }


  defaultp=function(p,Formulas,dval=NA) {
    default=Formulas
    for (i in names(Formulas)) default[[i]]=dval
    if (!is.null(p)) for (i in names(p)) {
      if (!(i %in% names(default)))
        stop("Parameter name not in Formulas arguement")   
      default[[i]]=p[[i]]
    }
    default
  }
  
  # fitlist=tmp; model.name=hmodelname
  get.start <- function(fitlist,model.name,pnams,Bounds) {
    # get generic (all cells) start points
    M <- fitlist[[model.name]]$M
    snams <- names(M)
    starts <- vector(mode="list",length=length(snams))
    names(starts) <- snams
    dat <- fitlist$dat
    for (i in snams) {
      Di <- M[[i]]
      dati <- dat[dat[[Di$S]]==i & is.finite(dat[,Di$RT]),]
      attr(dati,"D") <- Di
      starts[[i]] <- start.model(dati,pnams,Bounds)
    }
    starts
  }
    
  augment.design <- function(D,Latents,Cpars,SRpars=NULL) {
    if ( !is.null(Latents) ) {
      rdf <- is.data.frame(Latents[[1]]) # Non-factorial mapping
      if ( rdf ) {
        L.df <- Latents[[1]]
        Latents[[1]] <- as.character(levels(L.df[[2]]))
      }
      # Augment D with Latents to produce M    
      nL <- lapply(Latents,length)
      indx <- Latents[[1]]
      nfac <- length(Latents)
      if ( nfac>1 ) {
        for (i in 2:nfac) 
          indx <- outer(indx,Latents[[i]],"paste")
        out <- t(sapply(indx,function(x){strsplit(x," ")[[1]]}))   
      } else out <- matrix(indx,ncol=1)
      ncell <- dim(out)[1]
      dimnames(out) <- list(1:ncell,names(Latents))
      L <- data.frame(out[,dim(out)[2]:1,drop=F])
      nD <- dim(D[[1]]$D)[1]
      L <- L[rep(1:ncell,each=nD),,drop=F]
      for ( i in names(D) ) { # participants
#         
#         print(i)
#         
        D[[i]]$L=names(Latents)
        # Extract and augment
        Di <- D[[i]]$D
        for (j in 2:ncell) Di <- rbind(Di,D[[i]]$D)
        if ( !all(D[[i]]$CV=="") ) {
          CVi <- D[[i]]$CV
          for (j in 2:ncell) CVi <- rbind(CVi,D[[i]]$CV)          
        }
        Di$lcell=factor(rep(1:nD,times=ncell))
        Di <- cbind.data.frame(Di[,c("lcell","rcell",D[[i]]$SC)],
          L,Di[,c(D[[i]]$R,D[[i]]$F)])
        # latent fastest then response then facs
        reorder <- c(D[[i]]$L,D[[i]]$R,D[[i]]$F)
        reorder <- reorder[length(reorder):1] 
        tmp=Di
        for (j in names(tmp)) tmp[[j]] <- as.numeric(tmp[[j]]) 
        odr <- do.call(order,tmp)
        Di <- Di[odr,]
        dimnames(Di)[[1]]=1:dim(Di)[1]
        if ( !all(D[[i]]$CV=="") ) {
          CVi <- CVi[odr,,drop=F]
          dimnames(CVi)[[1]]=1:dim(CVi)[1]       
        }
        if ( rdf ) {
          ok <- logical(dim(Di)[1])
          snam <- names(L.df)[1]
          rnam <- names(L.df)[2]
          lrnam <- names(Latents)[1]
          for (j in 1:length(ok))
            ok[j] <- as.character(Di[[lrnam]][j]) %in% as.character(L.df[
              as.character(L.df[[snam]])==as.character(Di[[snam]][j]),rnam])
          Di <- Di[ok,]
          row.names(Di) <- 1:dim(Di)[1]
          if ( !all(D[[i]]$CV=="") ) {
            CVi <- CVi[ok,,drop=F]
            row.names(CVi) <- 1:dim(CVi)[1]
          }
        } else ok=rep(T,dim(Di)[2])
        D[[i]]$D <- Di
        if (!all(D[[i]]$CV=="")) D[[i]]$CV <- CVi
        # Model specific definition from seperate script (e.g., LBA, RD etc.)
        attr(D[[i]],"reorder") <- get.reorder(D=D[[i]],Cpars,ok,SR=SRpars)
      }    
    }
    D
  }
  
  make.star <- function(fveci) {
    if (is.logical(fveci)) return(logical(0))
    facs <- dimnames(attr(terms(formula(paste("~",fveci))),"factors"))[[1]]
    if (is.null(facs)) facs <- "1"
    if (length(facs)==1) return(paste("~",facs,sep=""))
    paste("~",paste(facs,collapse="*"),sep="")
  }

  make.all.star <- function(mnam) {
    vnam <- sapply(strsplit(mnam," & ")[[1]],function(x){
      strsplit(x,"~")[[1]][1]})
    snam <- sapply(strsplit(mnam," & ")[[1]],function(x){
      make.star(strsplit(x,"~")[[1]][2])})  	
  	paste(paste(vnam,snam,sep=""),collapse=" & ")
  }

  is.saturated <- function(fveci) {
    if (fveci == "1") return(T)
    facs <- attr(terms(formula(paste("~",fveci))),"term.labels")
    if (length(facs)==1) return(T)
    all(attr(terms(formula(paste("~",paste(facs,collapse="*")))),
    "term.labels") %in% facs)
  }    

  is.all.saturated <- function(mnam) {
    tmp <- strsplit(mnam," & ")[[1]]
    tmp <- sapply(tmp,function(x){strsplit(x,"~")[[1]][2]})
    all(unlist(sapply(tmp,is.saturated)))
  }  

  name.maps <- function(MT) {
    # changes map numbers back to names  
    add.names <- function(mlist,mnames) {mnames[mlist]}
  
    mnames=names(MT[[1]])
    MT <- lapply(MT,function(x){
      lapply(x,add.names,mnames=mnames)
    })
    MT
  }

  # Main Body

  if (!is.null(attr(fitlist$dat,"qp")) && any(is.na(rawdat)))
    stop("Must supply fitlist for raw data when creating models based on quantile data")
  fitlist[[hmodelname]] <- 
    list(M=augment.design(D=attr(fitlist$dat,"D"),Latents=Latents,Cpars=Cpars,SRpars=SRpars)) 	  
  stopFormulas <- defaultp(stopFormulas,Formulas,"~1")
  Contrasts <- defaultp(Contrasts,Formulas,dval=NA)
  Constants <- defaultp(Constants,Formulas,dval=NA)
  Bounds <- defaultp(Bounds,Formulas,dval=c(-Inf,Inf))
  if (!is.null(Constraints))
    Constraints <- defaultp(Constraints,Formulas,dval=NA)
  pnams <- names(Formulas) 
  npar <- length(pnams)
  if ( npar==1 ) stop("Function \"makehmaps\" needs > 1 parameter")
  flists <- vector(mode="list",length=npar)
  names(flists) <- pnams
  f <- flists  
  # hierarchy of formulas and immediate nesting maps in flists
  cat(paste("\nCreating hierarchies for each parameter",ifelse(saturated,
    "(SATURATED models only)","(Including ADDITIVE models)"),"\n"))
 
  for (i in pnams) {
#     print(i)
    flists[[i]] <- hforms(Formulas[[i]],stopFormulas[[i]],saturated=saturated)
  }
  # make matrix indx with rows = product of inidividual parameter
  # models and a column for each parameter type
  nout <- prod(unlist(lapply(flists,length)))
  flens <- lapply(flists,length)
  tmp <- lapply(flens,function(x){1:x})
  indx <- outer(tmp[[1]],tmp[[2]],"paste")
  if ( npar>2 ) 
    for ( i in pnams[-c(1,2)] ) indx <- outer(indx,tmp[[i]],"paste")
  indx <- t(sapply(as.vector(indx),
    function(x){as.numeric(strsplit(x," ")[[1]])}))
  dimnames(indx) <- list(1:nout,pnams)  
  # make out, a list with one entry for each parameter type
  # containing output of makemap  
  out <- vector(mode="list",length=nout)
  fs <- out
  cat(paste("\nMaking set of",nout,"models\n\n"))
  # initially make for just first subject
  for (i in 1:nout) {
    fi <- f
    for (j in pnams) fi[[j]] <- flists[[j]][[indx[i,j]]]
    fs[[i]] <- fi
    des <- fitlist[[hmodelname]]$M[[1]]$D
    if (!all(fitlist[[hmodelname]]$M[[1]]$CV==""))
      des <- cbind(des,fitlist[[hmodelname]]$M[[1]]$CV)
    out[[i]] <- makemap(design=des,Formulas=fi,
      Contrasts=Contrasts,Constants=Constants,Bounds=Bounds,
      Constraints=Constraints)   
  }
  cat(paste("TOP model has",attr(out[[1]],"npar"),"parameters\n\n")) 
  # get maps stored by hforms for each parameter
  mapslist <- lapply(fs,function(x){lapply(x,function(y){attr(y,"maps")})})
  # fill in last maps entry in ouput, indicates has no maps
  attr(out[[nout]],"maps") <- logical(0)
   
  # compute maps for each all parameter model in out
  cat("Creating start point maps (slow for large model trees)\n")
  repi <- pmax(1,nout%/%10)
  pcdone=10

  if ( nout>1 ) { # compute mapping if more than one model
    for (i in (nout-1):1) { # start from model above simplest
      okmap <- mapslist[[i]]
      for (j in pnams)  # store number of formulas for each parameter type
        if ( length(okmap[[j]])==0 ) okmap[[j]] <- flens[[j]]
      mindx <- outer(okmap[[1]],okmap[[2]],"paste")
      if ( npar>2 ) for ( j in pnams[-c(1,2)] ) 
        mindx <- outer(mindx,okmap[[j]],"paste")
      mindx <- t(sapply(as.vector(mindx),
        function(x){as.numeric(strsplit(x," ")[[1]])}))
      nm <- dim(mindx)[1]
      maps <- vector(length=0)
      compm <- matrix(indx[-i,],nrow=nout-1) # indices of models to compare to
      dimnames(compm) <- list(dimnames(indx)[[1]][-i],NULL)
      comp <-indx[i,]                       # current model index
      for (j in 1:nm) {
        mi <- mindx[j,]
        for (k in 1:npar) {
          tmp <- comp
          tmp[k] <- mi[k]
          maps <- c(maps,as.numeric(dimnames(compm)[[1]])[
            apply(compm,1,function(x){all(x==tmp)})])
        }
      }
      attr(out[[i]],"maps") <- sort(unique(maps))
      if ( (nout>100) && i%%repi==0 ) {
        cat(paste("   ",pcdone,"% done\n"))
        pcdone=pcdone+10
      }
    }  
  }
  
  # add names to entries in out, pnam~formula seperated by &
  tmp <- lapply(fs,paste)
  outnams <- vector(length=0)
  for (i in 1:length(out)) {
    tmp1 <- tmp[[i]]
    for (j in 1:length(pnams))
      tmp1[j] <- paste(pnams[j],
        paste(strsplit(tmp1[j]," ")[[1]],collapse=""),sep="")
    tmp1 <- paste(tmp1,collapse=" & ") 
    outnams <- c(outnams,tmp1)
  }
  star <- sapply(outnams,is.all.saturated)
  outnams[star] <- sapply(outnams[star],make.all.star) 
  names(out)=outnams
  cat("\nCreating model trees for subjects:\n")
  
  # Add in numeric starting points
  if ( is.null(attr(fitlist$dat,"qp")) ) 
    starts <- get.start(fitlist,model.name=hmodelname,pnams=pnams,Bounds=Bounds) else {
    tmp <- fitlist
    tmp$dat <- rawdat$dat
    starts <- get.start(tmp,model.name=hmodelname,pnams=pnams,Bounds=Bounds)
  }
  snams <- names(starts)
  MT <- vector(mode="list",length=length(snams))
  names(MT) <- snams
  mnams <- names(out)
  mlist <- vector(mode="list",length=length(mnams))
  names(mlist) <- mnams
  dname <- paste(dname,hmodelname,sep=".")   
  dir.create(dname)
  for (i in snams) {
    suppressWarnings({
      dir.create(paste(dname,"/",i,sep=""))
      dir.create(paste(dname,"/",i,"/fits",sep=""))      
      dir.create(paste(dname,"/",i,"/sims",sep=""))      
    })
    MT[[i]] <- mlist
    cat(paste(i,""))
    suppressWarnings({dir.create(paste(dname,"/",i,"/pmaps",sep=""))})
    p <- starts[[i]]      
    for (j in 1:nout) {
      fi <- f
      for (k in pnams) fi[[k]] <- flists[[k]][[indx[j,k]]]
      des <- fitlist[[hmodelname]]$M[[i]]$D
      if (!all(fitlist[[hmodelname]]$M[[i]]$CV==""))
        des <- cbind(des,fitlist[[hmodelname]]$M[[i]]$CV)
      pmap <- makemap(des,Formulas=fi,
        Contrasts=Contrasts,Constants=Constants,Bounds=Bounds,
        Constraints=Constraints)
      attr(pmap,"maps") <- attr(out[[j]],"maps")
      MT[[i]][[j]] <- attr(out[[j]],"maps")
      if (!all(names(pmap) %in% names(p)))
        stop("Parameter names in formula not present in model")
      attr(pmap,"start") <- getstart(p,pmap)
      pmap.dir <- paste(dname,"/",i,"/","pmaps","/",sep="")
      pmap.name <- paste(mnams[j],".RData",sep="")
      if (!update || (update & !any(dir(pmap.dir)==pmap.name)))
          save(pmap,file=paste(pmap.dir,pmap.name,sep=""))  
    } 
  }
  cat("\n")
  MT <- name.maps(MT) 
#  attr(MT,"dname") <- dname
  attr(MT,"stopFormulas") <- stopFormulas
  attr(MT,"Contrasts") <- Contrasts
  attr(MT,"Constants") <- Constants
  attr(MT,"Bounds") <- Bounds
  attr(MT,"Constraints") <- Constraints
  attr(MT,"saturated") <- saturated 
  fitlist[[hmodelname]]$MT <- MT
  # Print out model names and associated start point names
  cat("\n")
  nmods <- length(MT[[1]])
  mnams <- names(MT[[1]])
  nfits <- 1
  for ( i in 1:nmods ) {
    show <- showmodels | i==1 | i==length(MT[[1]])
    if (show) {
      cat(paste("Model",i,"\n"))
      cat(paste(mnams[i],"\n"))
    }
    if (length(MT[[1]][[i]])==0) {
      if (show) cat("No nested models\n")
    } else {
      if (show) cat("   Start points from: \n")
      for (j in MT[[1]][[i]]) {
        if (show) cat(paste("   ",j,"\n"))
        nfits=nfits+1
      }
    }
    if (show) cat("\n")
    if (!showmodels && i==2) cat("...\n\n")
  }
  cat(paste("NUMBER OF FITS REQUIRED PER SUBJECT =",nfits,"\n"))
  cat(paste("NUMBER OF FITS REQUIRED IN TOTAL =",nfits*length(MT),"\n"))
  save(fitlist,file=paste(dname,"/","fitlist.RData",sep=""))
  if ( any(dir(dname)=="fitstatus.RData") ) { # check if subjects have been added
    load(paste(dname,"fitstatus.RData",sep="/"))
    # Add in "start" as starting point model in list of models to be fit
    tofit1 <- fitlist[[hmodelname]]$MT
    for (i in names(tofit1)) tofit1[[i]][[length(tofit1[[i]])]] <- 
      c("start",tofit1[[i]][[length(tofit1[[i]])]])
    # Make a list of available starting models assuming no fits yet
    available1 <- as.vector(rep("start",length(tofit1)),mode="list")
    names(available1) <- names(tofit1) 
    for (i in names(available)) available1[[i]] <- available[[i]]
    for (i in names(tofit)) tofit1[[i]] <- tofit[[i]]
    available <- available1
    tofit <- tofit1
    save(available,tofit,file=paste(dname,"/","fitstatus.RData",sep=""))
  } 
}



# value=0;add.to.best=TRUE; constant=TRUE;constant.target=FALSE;subjects=NA
# add.bests <- function(dsource,dtarget,fac,value=0,add.to.best=TRUE,
#                       constant=TRUE,constant.target=FALSE,subjects=NA) {
#   # dsource and dtarget are names of source and target model tree directories
#   # fac = name of factor to be added
#   # add.to.best otherwise add to start
#   # value = value of start point for fac (on TRANSFORMED scale)
#   # constant = fac is in source model as a constant (or not in source model at all)
#   # constant.target = fac in target model to be made a constant (i.e., removed)
#   
#   # remove white space around factor specifications
#   stripwhite <- function(tms) {
#     for (i in names(tms) ) for (j in 1:length(tms[[i]]) )
#       if (j==1) tms[[i]][j] <- strsplit(tms[[i]][j]," ")[[1]][1] else
#         tms[[i]][j] <- strsplit(tms[[i]][j]," ")[[1]][2]
#     tms
#   }
#   
#   # find model(s) in sms matching tvec
#   tname <- function(tvec,sms) {
#     names(sms)[unlist(lapply(sms,function(x){
#       !any(is.na(sapply(x,match,table=tvec)))
#     }))]      
#   }
#   
#   fac1 <- paste(fac,"~1",sep="")
#   facI <- paste(fac,".I",sep="")
#   load(paste(dsource,"fitlist.RData",sep="/"))
#   sfl <- fitlist
#   load(paste(dtarget,"fitlist.RData",sep="/"))
#   tfl <- fitlist
#   if (any(is.na(subjects)))  subjects <- names(sfl[[2]]$MT)
#   subjects <- as.character(subjects)
#   if (!all(subjects %in% names(sfl[[2]]$MT)))
#       stop("Selected subjects not present in source")
#   for (s in subjects) {
#     if (!(s %in% names(tfl[[2]]$MT)) )
#       stop(paste(s,"not in target fitlist"))
#     cat(paste("Processing participant",s,"\n"))
#     # source model names
#     sms <- strsplit(names(sfl[[2]]$MT[[s]]),"&")
#     names(sms) <- names(sfl[[2]]$MT[[s]])
#     sms <- stripwhite(sms)
#     # target model names
#     tms <- strsplit(names(tfl[[2]]$MT[[s]]),"&")
#     names(tms) <- names(tfl[[2]]$MT[[s]])
#     tms <- stripwhite(tms)
#     if (!constant) for (m in names(tms)) {
#       is.fac <- fac1 == tms[[m]]
#       if (!any(is.fac)) 
#         stop(paste(fac,"not in model",m)) else 
#           tms[[m]] <- tms[[m]][!is.fac]
#     }
#     for (m in names(tms) ) { # loop over target models
#       sm <- tname(tms[[m]],sms)
#       if (length(sm)!=1)
#         stop(paste("Could not find unique model",m,"\n"))
#       strt <- attributes(get.pmap(dir.name=dsource,model=sm,subject=s))$best
#       if (is.null(strt))
#         stop("Best fitting parameters not in source")
#       pmap <- get.pmap(dir.name=dtarget,model=m,subject=s)
#       start <- attributes(pmap)$start
#       if ( !all(names(strt %in% names(start))) ) {
#         stop(paste("All source start names not in target start"))
#       } else {
#           if (constant.target) start <- strt[names(start)] else {
#             start[names(strt)] <- strt
#             is.facI <- names(start)==facI
#             if (!any(is.facI))
#               stop(paste("Target start vector does not contain",facI))
#             start[is.facI] <- value
#           }
#           if (add.to.best)
#             attributes(pmap)$best <- start else
#             attributes(pmap)$start <- start  
#       }
#       save.pmap(pmap,dir.name=dtarget,model=m,subject=s)
#     }
#   }  
# }



# value=0;add.to.best=TRUE; constant=TRUE;constant.target=FALSE;subjects=NA
# 
# dsource="hayesbv.LBAe";dtarget="hayes.LBAe"; constant=TRUE; add.to.best=TRUE
# fac <- c("bT","cT"); value=c(log(1e-6),log(1e-6))

add.bests <- function(dsource,dtarget,fac,value=0,add.to.best=TRUE,
                      constant=TRUE,constant.target=FALSE,subjects=NA) {
  # dsource and dtarget are names of source and target model tree directories
  # fac = name of factor to be added (can be vector)
  # add.to.best otherwise add to start (must be same length as fac)
  # value = value of start point for fac (on TRANSFORMED scale)
  # constant = fac is in source model as a constant (or not in source model at all)
  # constant.target = fac in target model to be made a constant (i.e., removed)
  
  # remove white space around factor specifications
  stripwhite <- function(tms) {
    for (i in names(tms) ) for (j in 1:length(tms[[i]]) )
      if (j==1) tms[[i]][j] <- strsplit(tms[[i]][j]," ")[[1]][1] else
        tms[[i]][j] <- strsplit(tms[[i]][j]," ")[[1]][2]
    tms
  }
  
  # find model(s) in sms matching tvec
  tname <- function(tvec,sms) {
    names(sms)[unlist(lapply(sms,function(x){
      !any(is.na(sapply(x,match,table=tvec)))
    }))]      
  }
  
  fac1 <- paste(fac,"~1",sep="")
  facI <- paste(fac,".I",sep="")
  load(paste(dsource,"fitlist.RData",sep="/"))
  sfl <- fitlist
  load(paste(dtarget,"fitlist.RData",sep="/"))
  tfl <- fitlist
  if (any(is.na(subjects)))  subjects <- names(sfl[[2]]$MT)
  subjects <- as.character(subjects)
  if (!all(subjects %in% names(sfl[[2]]$MT)))
      stop("Selected subjects not present in source")
  for (s in subjects) {
    if (!(s %in% names(tfl[[2]]$MT)) )
      stop(paste(s,"not in target fitlist"))
    cat(paste("Processing participant",s,"\n"))
    # source model names
    sms <- strsplit(names(sfl[[2]]$MT[[s]]),"&")
    names(sms) <- names(sfl[[2]]$MT[[s]])
    sms <- stripwhite(sms)
    # target model names
    tms <- strsplit(names(tfl[[2]]$MT[[s]]),"&")
    names(tms) <- names(tfl[[2]]$MT[[s]])
    tms <- stripwhite(tms)
    if ( !constant ) for (m in names(tms)) {
      is.fac <-  tms[[m]] %in% fac1 # sapply(fac1,function(x){x==tms[[m]]})
      if ( !any(is.fac) ) 
        stop(paste(fac,"not in model",m)) else 
          tms[[m]] <- tms[[m]][!is.fac]
    }
    for (m in names(tms) ) { # loop over target models
      sm <- tname(tms[[m]],sms)
      if (length(sm)!=1)
        stop(paste("Could not find unique model",m,"\n"))
      strt <- attributes(get.pmap(dir.name=dsource,model=sm,subject=s))$best
      if (is.null(strt))
        stop("Best fitting parameters not in source")
      pmap <- get.pmap(dir.name=dtarget,model=m,subject=s)
      start <- attributes(pmap)$start
      if ( !all(names(strt %in% names(start))) ) {
        stop(paste("All source start names not in target start"))
      } else {
          if ( constant.target ) start <- strt[names(start)] else {
            start[names(strt)] <- strt
            is.facI <- names(start) %in% facI
            if ( !any(is.facI) )
              stop(paste("Target start vector does not contain",facI))
            start[is.facI] <- value
          }
          if (add.to.best)
            attributes(pmap)$best <- start else
            attributes(pmap)$start <- start  
      }
      save.pmap(pmap,dir.name=dtarget,model=m,subject=s)
    }
  }  
}

# dname="simT.LBA"; mname="B~SA & A~1 & v~C & sv~C & ter~SA & pc~1"; start="sim.LBA"
add.starts <- function(dname,mname,start) { 
  start1 <- start
  for ( s in dir(dname)[-grep("RData",dir(dname))] ) {
     if ( is.character(start1) ) { # add in best from another fit directory
       load(paste(paste(start1,s,"pmaps",mname,sep="/"),"RData",sep=".")) 
       start <- attributes(pmap)$best
     }
     pmap.name <- paste(dname,s,"pmaps",paste(mname,"RData",sep="."),sep="/")
     load(pmap.name)
     nams <- names(attributes(pmap)$start)
     if ( is.null(dim(start)) ) {
       names(start) <- nams
       attributes(pmap)$start <- start 
     } else {
        dimnames(start)[[2]] <- nams
        attributes(pmap)$start <- start[s,]
     }
     save(pmap,file=pmap.name)
   } 
}
