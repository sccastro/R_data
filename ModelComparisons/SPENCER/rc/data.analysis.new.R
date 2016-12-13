########################### DATA ANALYSIS UTILITIES ############################

wsAnova=function(dat,SStype=3) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  dat[,1] <- factor(dat[,1])
  dat=dat[do.call(order,dat[,-dim(dat)[2]]),]
  snams=levels(dat[,1]); ns=length(snams)
  dvnam=names(dat)[dim(dat)[2]]
  facn=names(dat)[-c(1,dim(dat)[2])]
  nifac=length(facn)
  idata=data.frame(dat[dat[,1]==snams[1],facn])
  names(idata)=facn
  for (i in facn) 
    if (i==facn[1]) ifacn=as.character(idata[,1]) else 
                    ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
  facnr=facn[nifac:1]
  e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
              ncol=length(ifacn),dimnames=list(snams,ifacn))
  out = try(summary(Anova(lm(e.mv ~ 1),
    idata=idata,type=SStype, 
    idesign=formula(paste("~",paste(facn,collapse="*")))),
    multivariate=FALSE))
}

mixedAnova=function(dat,bsfacn,wsfacn=NULL,sfac="s",SStype=3,spss=F,dvnam=NA) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  if (is.na(dvnam)) dvnam <- names(dat)[dim(dat)[2]]
  dat <- dat[,c(sfac,bsfacn,wsfacn,dvnam)]
  if (length(wsfacn)>0) dat <- dat[do.call(order,dat[,c(sfac,wsfacn)]),]
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  snams=levels(dat[,sfac]); ns=length(snams)
  nifac=length(wsfacn)
  lev1s=unlist(lapply(lapply(dat[,wsfacn],levels),function(x){x[1]}))
  bsfacs=dat[apply(dat[,wsfacn,drop=F],1,function(x){all(x==lev1s)}),bsfacn,drop=F]  
  if ( nifac>0 ) {
    idata=data.frame(dat[dat[,sfac]==snams[1],wsfacn])
    names(idata)=wsfacn
    for (i in wsfacn) 
      if (i==wsfacn[1]) ifacn=as.character(idata[,1]) else 
        ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
    e.mv=matrix(unlist(
      tapply(dat[,dvnam],dat[,wsfacn[length(wsfacn):1]],function(x){x})),
                ncol=length(ifacn),dimnames=list(snams,ifacn))
    summary(Anova(
      lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),bsfacs),
      idata=idata,type=SStype,
      idesign=formula(paste("~",paste(wsfacn,collapse="*")))),multivariate=F)
  } else {
    e.mv <- cbind(y=dat[,dvnam]) 
    print(Anova(lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),
                   bsfacs),type=3))
  }
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),bsfacs,e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}

mneffects=function(df,elist,digits=3,err=F,vars=F,save=F,stat="mean") {
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  for (i in 1:length(elist)) {
    cat(paste(paste(elist[[i]],collapse=":"),"\n"))
    mns=tapply(df[,dvnam],df[,elist[[i]]],stat,na.rm=T)
    if (err) print(round(plogis(mns),digits)) else
      if (vars) print(round(sqrt(mns),digits))  else
        print(round(mns,digits))    
    cat("\n")
  } 
  if (save) mns
}

se=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  se=tapply(df[,dvnam],df[,facs],sd)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (ci=="SE") se/sqrt(ns) else
    qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
}

add.bars=function(mn,se,xvals=NA,len=.1,antiprobit=FALSE,col="black") {
  
  plotbars <- function(x,m,l,h,len,col="black") {
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],l[j],length=len,angle=90,col=col)
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],h[j],length=len,angle=90,col=col)    
  }
  
  if (any(is.na(xvals))) if (is.matrix(mn)) 
    xvals <- as.numeric(dimnames(mn)[[2]]) else
      xvals <- as.numeric(names(mn))
  lo <- mn-se
  hi <- mn+se
  if (antiprobit) {
    mn=pnorm(mn)
    lo=pnorm(lo)
    hi=pnorm(hi)
  }
  if (!is.matrix(mn)) 
    plotbars(xvals,mn,lo,hi,col=col,len=len) else
      for (i in 1:dim(mn)[1]) 
        plotbars(x=xvals,m=mn[i,],l=lo[i,],h=hi[i,],len=len,col=col)
}    


############################ CHECK AND IMPORT DATA #############################
#
# Data imported to an object with S3 class "rc" (rapid choice) but
# class based mechanisms not yet implemented.
# 



# makes a data.frame specifying non-factorial response assignment
# used as input to score.rc R parameter
# snam is name of an element of score.rc F parameter, alldat is data frame passed
# to score.rc and rnam is name of a column in alldat specifying responses
make.R <- function(alldat,snam="S",rnam="R") {
  out <- t(matrix(unlist(strsplit(
    outer(levels(alldat[[snam]]),levels(alldat[[rnam]]),paste),
    " ")),nrow=2))
  names(out) <- c(snam,rnam)
  ok <- logical(dim(out)[1])
  cat("First element is condition, second is repsonse\n")
  for (i in 1:length(ok)) {
    print(out[i,])
    ok[i] <- as.logical(menu(c("No","Yes"))-1)
  }
  R.df <- data.frame(out[ok,])
  names(R.df) <- c(snam,rnam)
  R.df
}


# Prepares a data frame for use by rc

# S=NULL;F=NULL;RT=NULL;CV=NULL;SC=NULL;cv=NULL;R=NULL  
# default.scores=list();autoscore=list();autoscore.list=list()
# save.trials=NA;precision.rt=.001
# 
# alldat=datPM;S="s";R="R";RT="RT";SC=c("C","SR");F=c("PM","S")
# autoscore=list(C=c(sfac="S",rfac="R"),SR=list(sfac=c("S","R"))) 
# autoscore.list=list(SR=list(
#   nonword.nonword="nn",nonword.word="nw",nonword.pm="pfa",
#   word.nonword="wn",word.word="ww",word.pm="pfa",
#   pm.nonword="pn",pm.word="pw",pm.pm="pp"))     

# alldat=dat1;S="s";R="R";RT="RT";F=c("RS","IS");CV="CI"
# SC=c("C","SR");autoscore=list(C=c(sfac="RS",rfac="R"),SR=list(sfac=c("RS","IS","R"))) 
# autoscore.list=list(SR=SR)

score.rc <- function(alldat,
  S=NULL,F=NULL,RT=NULL,CV=NULL,SC=NULL,cv=NULL,
  R=NULL, # character or R.df data frame made by R.df
  # Following 3 have can have one entry for each level of SC
  default.scores=list(), # e.g. default.scores=list(C="FALSE")) 
  autoscore=list(),      # e.g. for CV a boolean autoscore=list(C=list(sfac="S",rfac="R"))
  #  e.g., for CV a factor(s) autoscore=list(SR=list(sfac=c("S","R"))) (no rfac needed)
  autoscore.list=list(), # specify 1:many R- > S mapping, 
  #  e.g., where CV is a boolean list(C=list(word=c("HFW","LFW"),nonword=c("HNN","LNN")))
  #  e.g., where CV is a factor (NB. the names in the list must be a paste of the 
  #        levels in the columns specified in sfac in autoscore)
  #    list(SR=list(nonword.nonword="nn",nonword.word="nw",nonword.pm="pfa",
  #      word.nonword="wn",word.word="ww",word.pm="pfa",
  #      pm.nonword="pn",pm.word="pw",pm.pm="pp")) 
  save.trials=NA,      # keep trial factor(s), e.g., save.trials=c("cycle","trial")
  precision.rt=.001) {
  # alldat = data frame
  # S,R,F,RT,CV,SC: subject, response,factor,RT,covariate and score column names
  # cv names data level (per trial) covariates
  # default.scores: named vector of values to fill in if score missing
  
  add.design <- function(alldat,F,R,cv) {  
  	rdf <- is.data.frame(R)
    if (rdf) {
      R.df <- R
      R <- names(R)[2]
  	}
    nR <- length(levels(alldat[[R]]))
    faclevs <- lapply(alldat[,c(R,F)],levels)
    indx <- faclevs[[1]]
    nfac <- length(faclevs)
    for (i in 2:nfac) indx <- outer(indx,faclevs[[i]],"paste")
    out <- t(sapply(indx,function(x){strsplit(x," ")[[1]]}))
    ncell <- dim(out)[1]
    dimnames(out) <- list(1:ncell,c(R,F))
    D=data.frame(out[,dim(out)[2]:1],stringsAsFactors=FALSE)
    for (i in c(F,R)) D[[i]] <- factor(D[[i]],levels(alldat[[i]]))
    if (rdf) {
      ok <- logical(dim(D)[1])
      snam <- names(R.df)[1]
      for (i in 1:length(ok))
        ok[i] <- as.character(D[[R]][i]) %in% as.character(R.df[
          as.character(R.df[[snam]])==as.character(D[[snam]][i]),R])
      D <- D[ok,]
      ncell <- dim(D)[1]
      row.names(D) <- 1:ncell
    }
    # insert cell and rcell indicator in data
    cell=rep(NA,dim(alldat)[1])
    for (i in 1:dim(D)[1]) {  
      ok=rep(T,length=dim(alldat)[1])
      for (j in dimnames(D)[[2]]) 
         ok <- ok & as.character(alldat[,j])==as.character(D[i,j])
      cell[ok]=i
    }
    if (rdf) {
      if (dim(D)[2]==2) cnams <- as.character(D[,1]) else
        cnams <- apply(D[,-dim(D)[2]],1,paste,collapse="")
      rcell <- integer(dim(D)[1])
      ucnams <- unique(cnams)
      for (i in 1:length(ucnams)) rcell[cnams==ucnams[i]] <- i       
    } else rcell <- 1 + c(0:(ncell-1)) %/% nR
  	D=cbind.data.frame(D,rcell=factor(rcell))
    cat("Added the following manifest design\n")
    print(D)
    cat("\n")
    rcell= 1 + ((cell-1) %/% nR)
    dat <- cbind.data.frame(cell=factor(cell,row.names(D)),rcell=factor(D[cell,"rcell"]),alldat)
    if (!all(levels(dat$rcell)==levels(D$rcell)))
      stop("Some response cells have no response!")
    attr(dat,"D") <- list(D=D[,(dim(D)[2]):1],F=F,R=R,cv=cv) 
    dat
  }
  
  score <- function(dat,S,SC,CV,RT,default.scores=NA,autoscore=NA,autoscore.list=NA) {
    D=attr(dat,"D")
    D$SC=SC
    if ( is.null(CV) ) D$CV="" else {
      D$CV <- data.frame(D$D[,rep(D$F[1],length(CV))])
      names(D$CV) <- CV
      for (i in CV) D$CV[,i] <- as.character(D$CV[,i])
      for (i in 1:dim(D$D)[1]) {
        flevs <- D$D[i,]
        for (j in D$F) if (j==D$F[1]) 
          isin <- dat[,j]==as.character(flevs[,j]) else
            isin <- isin & dat[,j]==as.character(flevs[,j])        
        for (j in CV) {
          lev <- unique(as.character(dat[isin,j]))
          if (length(lev)>1)
            stop(paste("Covariate",j,"does not uniquely map to design")) else
            D$CV[i,j] <- lev
        }
      }
      for (i in CV) D$CV[,i] <- factor(D$CV[,i])
    }
    D$S=S
    D$RT=RT
    ncell <- dim(D$D)[1]
    slevs <- levels(dat[[S]])
    out <- vector(mode="list",length=length(slevs))
    names(out) <- slevs
    for (i in slevs) out[[i]]=D
    for (i in SC) {
      cat(paste("SCORING",i,"\n\n"))
      if ( is.null(autoscore.list[[i]]) ) 
        okscore <- sort(unique(dat[,i])) else
        okscore <- sort(unlist(unique(autoscore.list[[i]])))           
      boolean.score <- is.logical(dat[,i])
      for (j in slevs) {
        cat(paste("Scoring subject",j,"\n"))
        Di=out[[j]]
        dati <- dat[dat[[S]]==j,]
        
        # Get scoring (NA if not defined) and check is unambiguous
        score <- vector(length=0)
        for ( k in 1:ncell ) { 
          scorek <- dati[dati$cell==k,i]
          if ( length(scorek)==0 ) score=c(score,"NA") else {
            if ( length(unique(scorek))!=1 )
              stop(paste("Ambiguous score",i,"for subject",j,
                "must be fixed in data")) else
              score <- c(score,as.character(unique(scorek))) 
          }  
        }	
        score[any(is.na(score))] <- "NA"
        
        Di$D <- cbind(Di$D,score) # score=as.character(score)
        names(Di$D)[dim(Di$D)[2]] <- i
        badscore <- score == "NA"
        if ( any(badscore) ) {
          Di$D[,i] <- factor(as.character(Di$D[,i]),levels=okscore)
          if ( !is.null(default.scores[[i]]) ) {
            print(Di$D[is.na(Di$D[[i]]),])
            Di$D[[i]][badscore] <- default.scores[i] 
            cat(paste("\nDefault score added to subject",j,"in cell(s)",
              paste(c(1:ncell)[badscore],collapse=" "),"\n\n"))
            cat("\n")
          } else if ( !any(is.null(autoscore[[i]])) ) {
            cat(paste("Autoscoring",dim(Di$D[is.na(Di$D[[i]]),])[1],
              "cells for participant",j,"\n\n"))
            if ( !any(is.null(autoscore.list[[i]])) ) {
              if ( boolean.score ) {
                for ( rlev in levels(Di$D[[ autoscore[[i]][["rfac"]] ]]) ) {
                  isin <- Di$D[[autoscore[[i]][["rfac"]]]]==rlev
                    Di$D[is.na(Di$D[[i]]) & isin,i] <- 
                      Di$D[is.na(Di$D[[i]]) & isin,autoscore[[i]][["sfac"]]] %in% 
                        autoscore.list[[i]][[rlev]]
                } 
              } else {
                Di$D[is.na(Di$D[[i]]),i] <- unlist(autoscore.list[[i]])[
                  apply(Di$D[ is.na(Di$D[[i]]),autoscore[[i]][["sfac"]] ],
                      1,paste,collapse=".")]
              }
            } else {
              if ( boolean.score ) {
                Di$D[is.na(Di$D[[i]]),i] <- 
                  Di$D[ is.na(Di$D[[i]]),autoscore[[i]][["sfac"]] ]==
                    Di$D[is.na(Di$D[[i]]),autoscore[[i]][["rfac"]]]
              } else
                stop("Must provide an autoscore.list to score a non-boolean CV")
            }
          } else { 
            cat(paste("Problem in scoring participant",j,"\n\n"))
            print(Di$D[is.na(Di$D[[i]]),])
            for (l in c(1:ncell)[badscore] ) {
              cat(paste("\nScore for cell",l,"\n"))
              Di$D[l,i] <- as.character(okscore[menu(okscore)])
            }
          }  
        } 	
        Di$D[[i]]=factor(as.character(Di$D[[i]]))
        out[[j]]=Di
      }
    }
    attr(dat,"D")=out
    dat 	
  }    
  
  if (is.null(F)) 
    stop("Must specify design factors as character vector \"F\"")
  if (!all(F %in% names(alldat)))
    stop("Some factor column names doe not exist in data")
  if (!is.null(CV) && !all(CV %in% names(alldat)))
    stop("Some covariate column names do not exist in data")
  if (!is.null(cv) && !all(cv %in% names(alldat)))
    stop("Some data covariate column names do not exist in data")
  if (!any(names(alldat)==S))
    stop("Subject column name does not exist in data")
  if (!any(names(alldat)==RT))
    stop("Response time column name does not exist in data")
  if (is.data.frame(R)) {
    if (!any(names(R)[1] %in% F))
      stop("Response scoping column name not in F")
    if (!any(names(R)[2] %in% names(alldat)))
      stop("Response column name does not exist in data")
  } else if (!any(names(alldat)==R))
    stop("Response column name does not exist in data")
  if (any("rcell" %in% c(S,R,F,RT,CV,SC)))
    stop("Cannot use reserved name \"rcell\" in any names")
  if (any("cell" %in% c(S,R,F,RT,CV,SC)))
    stop("Cannot use reserved name \"cell\" in any names")
  if (any("qp" %in% c(S,R,F,RT,CV,SC)))
    stop("Cannot use reserved name \"qp\" in any names")
  if (any("qn" %in% c(S,R,F,RT,CV,SC)))
    stop("Cannot use reserved name \"qn\" in any names")
  # enforce precision
  alldat[,RT] <- round(alldat[,RT]/precision.rt)*precision.rt 
  # spread ties evenly over precision interval
  is.dup <- duplicated(alldat[,RT])
  is.dup[is.na(alldat[,RT])] <- FALSE
  dups <- unique(alldat[,RT][is.dup])
  if (any(is.dup)) cat(paste("\nSpreading",sum(is.dup)+length(dups),"of",
    dim(alldat)[1],"RTs that are ties given preceision",precision.rt,
    ".\n  ",length(dups),"have ties out of",length(unique(alldat[,RT])),"unique values\n\n"))
  if (length(dups)>0) for (i in dups) {
    is.in <- alldat[,RT]==i
    is.in[is.na(is.in)] <- FALSE
    n.dup <- sum(is.in)    
    alldat[is.in,RT] <- 
      c(i-precision.rt/2 + precision.rt*c(1:n.dup)/(n.dup+1))[order(runif(n.dup))]
  } 
  for (i in SC) 
    if ( any(is.na(alldat[[i]])) & is.null(autoscore[[i]]) ) 
      stop(paste("Cannot have \"NA\" in a score column (",i,")"))
  dat <- add.design(alldat=alldat,F=F,R=R,cv=cv)
  dat <- score(dat=dat,S=S,SC=SC,CV=CV,RT=RT,default.scores=default.scores,
    autoscore=autoscore,autoscore.list=autoscore.list)
  if (is.data.frame(R))  R <- names(R)[2]
  tmp <- dat[,c("cell","rcell",S,F,CV,cv,SC,R,RT)]
  attr(tmp,"D") <- attr(dat,"D")
  if (!(any(is.na(save.trials)))) attr(tmp,"trials") <- dat[,save.trials]
  tmp
}



# dat=seq11c.df;correct.name="C";minrt=.2;maxrt=3;qp=NA
# save.trials=T
# sremove=NULL; fremove=NULL
# show.summary=TRUE;show.censors=FALSE
# show.cell.censor=FALSE; show.cell.censors=FALSE
make.rc <- function(dat,      # data frame prepared by score.rc()
  minrt=.2,maxrt=2.5,         # Censoring limits
  qp=NA,                      # Quantile probabilities
  correct.name="C",           # Name of factor scoring correct
  save.trials=F,              # Add in an attribute to identify trials 
  fremove=NULL,               # Factors in dat to NOT use
  sremove=NULL,               # Subjects in dat to NOT use
  show.summary=TRUE,          # Overall censoring and cell size summary 
  show.censors=FALSE,         # Count of number censored for each subject
  show.cell.censors=FALSE,    #   add breakdown by cell per subject 
  show.cell.censor=FALSE) {   # Censoring per cell summed over subjects
  # takes a data frame prepared by score.rc() and makes a class rc
  # object (list with data frame as first entry)
  # NOTE
  # correct.name: column name in dat (must convert nicely with as.logical)
  # Censors <=minrt and >=maxrt, and records these limits and 
  #   cenored counts as an attribute "nc"
  # qp: If !is.na codes data as specified quantiles
 
  fix.dat <- function(dat) {  
    F=attr(dat,"D")[[1]]$F
  	R=attr(dat,"D")[[1]]$R
  	nR <- length(levels(dat[[R]]))
    faclevs <- lapply(dat[,c(R,F)],levels)
    indx <- faclevs[[1]]
    nfac <- length(faclevs)
    for (i in 2:nfac) indx <- outer(indx,faclevs[[i]],"paste")
    out <- t(sapply(indx,function(x){strsplit(x," ")[[1]]}))
    ncell <- dim(out)[1]
    dimnames(out) <- list(1:ncell,c(R,F))
    D=data.frame(out[,dim(out)[2]:1])
    # insert cell and rcell indicator in data
    cell=rep(NA,dim(dat)[1])
    for (i in 1:dim(D)[1]) {  
      ok=rep(T,length=dim(dat)[1])
      for (j in dimnames(D)[[2]]) 
         ok <- ok & as.character(dat[,j])==as.character(D[i,j])
      cell[ok]=i
    }
    rcell= 1 + ((cell-1) %/% nR)
    dat$cell <- factor(cell)
    dat$rcell <- factor(rcell)    
    dat
  }

  # Main Body
  if (length(qp)>1 & any(is.na(qp)))
    stop("qp must be either a single NA or contain no NAs")
  D <- attr(dat,"D")
  if (!is.null(fremove)) {
    faclist <- D[[1]]$F	
  	if (!all(fremove %in% faclist))
  	  stop("Some factor(s) to remove not in rc object")
  	cat("\nReducing design\n")
  	ok <- rep(T,length(names(dat)))
  	names(ok) <- names(dat)
  	ok[names(ok) %in% fremove] <- FALSE
  	dat <- dat[,ok]
  	dat$cell <- as.character(dat$cell)
  	dat$rcell <- as.character(dat$rcell)
  	for (i in names(D)) {
  	  Di <- D[[i]]
  	  faclist <- Di$F[!(Di$F %in% fremove)]
  	  Di$F <- faclist
  	  for (j in fremove) {
  	  	l1 <- levels(Di$D[[j]])[1]
  	  	ok <- !(names(Di$D)==j)
  	  	Di$D <- Di$D[Di$D[[j]]==l1,ok] 
  	  }
  	  Di$D$rcell <- factor(as.character(Di$D$rcell),sort(unique(Di$D$rcell)))
  	  levels(Di$D$rcell) <- 1:length(levels(Di$D$rcell))
  	  D[[i]] <- Di
  	}
#    if (!any(is.na(qp))) D[[i]]$qp <- qp 
  	attr(dat,"D") <- D
  	dat <- fix.dat(dat)
  }
  snams <- names(D)
  if (!is.null(sremove)) snams <- snams[!(snams %in% sremove)]
  ns <- length(snams)
  if (length(minrt)==1) {
    minrt=rep(minrt,ns)
    names(minrt) <- snams
  }
  if (length(maxrt)==1) {
    maxrt=rep(maxrt,ns)
    names(maxrt) <- snams
  }
  if (length(minrt)!=ns) 
    stop("Minimum has to be length 1 or length=number of subjects")
  if (length(maxrt)!=ns) 
    stop("Maximum has to be length 1 or length=number of subjects")
  if (any(maxrt-minrt<=0))
    stop("Maximum RT <= minimum RT")
  if (any(!is.finite(maxrt)) | any(!is.finite(minrt)))
    stop("Maximum and minimum RTs must be finite")
  Dnew <- vector(mode="list",length=length(snams))
  names(Dnew) <- snams
  for (i in snams) {
    Di <- D[[i]]
    dati <- dat[dat[[Di$S]]==i,]
    lo <- dati[[Di$RT]]<=minrt[i]
    lo[is.na(lo)] <- FALSE
    hi <- dati[[Di$RT]]>=maxrt[i]
    hi[is.na(hi)] <- FALSE
    faclist <- dati[,c(Di$R,Di$F)]
    ncmin <- tapply(lo,faclist,sum)
    ncmin[is.na(ncmin)] <- 0
    ncmax <- tapply(hi,faclist,sum)
    ncmax[is.na(ncmax)] <- 0      
    dn <- dimnames(ncmin)
    dn[["minmax"]] <- c(as.character(minrt[i]),as.character(maxrt[i]))
    Di$nc <- array(c(ncmin,ncmax),dim=c(dim(ncmin),2),dimnames=dn)
    if (show.censors) {
      cat(paste("Processed participant",i,"  Censored=",sum(Di$nc),"\n"))
      if (show.cell.censors) print(Di$nc)
    }
    Dnew[[i]] <- Di
    if (i==snams[1]) {
      datc <- dati[!lo &!hi,] 
      if (save.trials)
        trials <- attr(dat,"trials")[dat[[Di$S]]==i,][!lo &!hi,]
    } else {
      datc <- rbind.data.frame(datc,dati[!lo &!hi,])
      if (save.trials) trials <- rbind(trials,
        attr(dat,"trials")[dat[[Di$S]]==i,][!lo &!hi,])
    }
  }
  dat <- datc
  D <- Dnew
  dat[,Di$S] <- factor(as.character(dat[,Di$S]),levels=snams)
  attr(dat,"D") <- D
  cdat <- dat[as.logical(dat[[correct.name]]),]
  edat <- dat[!as.logical(dat[[correct.name]]),]
  Di <- D[[1]]
  nc <- D[[snams[1]]]$nc
  for (i in snams[-1]) nc=nc+D[[i]]$nc
  dn <- dimnames(nc)
  minmax <- dn[[length(dn)]]
  tmp=100*unlist(lapply(D,function(x){sum(x$nc)}))/
    (tapply(dat[[Di$RT]],dat[[Di$S]],length)+
     unlist(lapply(D,function(x){sum(x$nc)})))
  if (show.summary) {
    cat("\nOverall % RTs censored per subject\n")
    print(round(sort(tmp,decreasing=T),1))
  }
  if (show.cell.censor) {
    cat("\nOverall number of RTs censored per cell\n")
    print(nc)
  }
  if (show.summary) { 
    cat(paste("Overall percent censored:",
      round(100*sum(nc)/dim(dat)[1],1),"\n"))
    # min num correct
    cat(paste("\nMin/Median/Max number of correct responses per cell\n"))
    tmp <- tapply(cdat[[Di$RT]],cdat[,c(Di$S,Di$F)],length)
    tmp[is.na(tmp)] <- 0
    print(c(min(tmp),round(median(tmp)),max(tmp)))
    cat(paste("\nMin/Median/Max number of error responses per cell\n"))
    tmp <- tapply(edat[[Di$RT]],edat[,c(Di$S,Di$F)],length)
    tmp[is.na(tmp)] <- 0
    print(c(min(tmp),round(median(tmp)),max(tmp)))
  }
  # put into sorted (subject, cell, RT) quantile form
  odr <- do.call(order,dat[,c(D[[1]]$S[1],"cell",D[[1]]$RT)])
  if (save.trials) trials <- trials[odr,]
  dat <- dat[odr,]
  dat$qp <- dat[,D[[1]]$RT]
  dat$qn <- dat[,D[[1]]$RT]
  qdat <- dat[NULL,]
  cat("Processing participant: ")
  for (i in levels(dat[,D[[1]]$S[1]]) ) {
     cat(paste(i,""))
     for (j in levels(dat$cell) ) {
      datij <- dat[dat[,D[[1]]$S[1]]==i & dat$cell==j,]
      nij <- dim(datij)[1]
      if (nij>0) { # observations in cell
        if ( any(is.na(qp)) ) {  # all data
          datij$qp <- (1:nij)/(nij+1)
          datij <- rbind(datij,datij[1,])
          datij[nij+1,D[[1]]$RT] <- Inf
          datij$qp[nij+1] <- 1
          datij$qn <- nij/(nij+1)
          qdat <- rbind(qdat,datij)
        } else {               # quantile reduce data
          # Use on qp internal to empirical cdf quantiles
          eqpij <- c(1:nij)/(nij+1)
          qpij <- qp[qp>=min(eqpij) & qp<=max(eqpij)] 
          nq <- length(qpij)
          q <- quantile(datij[,D[[1]]$RT],qpij,type=6)
          if (any(duplicated(q))) {
            print(sort(datij[,D[[1]]$RT]))
            print(qpij)
            print(q)
            stop("Duplicated quantiles!")
          }
          # Note following is > lower bound and <= upper bound
          n <- diff(c(0,qpij,1))*nij
          datij <- datij[rep(1,nq+1),]
          datij$qn <- n
          datij$qp <- c(qpij,1)
          datij[,D[[1]]$RT] <- c(q,Inf)
          qdat <- rbind(qdat,datij)
        }
      }    
    }
  }
  cat("\n")
  row.names(qdat) <- NULL
  for (i in levels(qdat[,D[[1]]$S[1]])) {
    iss <- qdat[,D[[1]]$S[1]]==i
    for (j in levels(qdat$rcell)) {
      datij <- qdat[iss & qdat$rcell==j,]
      datij$cell <- factor(as.character(datij$cell))
      nij <- tapply(datij$qn,datij$cell,sum)
      pij <- nij/sum(nij)
      for (k in names(pij))
        qdat[iss & qdat$cell==k & qdat$qp==1,"qp"] <- pij[k]
    }
  }
  out <- list(dat=qdat)
  if (!any(is.na(qp))) attr(out$dat,"qp") <- qp
  class(out) <- c("rc","list")
  if (save.trials) attr(out,"trials") <- trials
  out
}

####################################################################################

get.simdat <- function(dname,mname,ss=NA,excludeS=NA) {
  load(paste(dname,"fitlist.RData",sep="/"))
  attr(fitlist[[2]]$MT,"dname") <- dname
  dat <- fitlist$dat
  D <- attr(dat,"D")[[1]]
  snams <- names(fitlist[[2]]$M)
  if (any(is.na(ss))) ss <- snams else 
    if (is.numeric(ss)) ss <- snams[ss]
  if (!all(ss %in% snams))
    stop("Some subject names not in data file")
  if (!all(is.na(excludeS))) {
    ss <- ss[!(ss %in% excludeS)]
  }
  dat <- dat[dat[[D$S]] %in% ss,]
  dat[[D$S]] <- factor(as.character(dat[[D$S]]))
  RTnam <- D$RT
  rtsimnam <- paste(RTnam,"sim",sep=".")
  RTsim <- paste(RTnam,"sim",sep=".")
  dat$lm <- dat[[RTnam]]
  dat$lc <- dat[[RTnam]]
  dat[[RTsim]] <- dat[[RTnam]]
  dat$qp.sim <- dat[[RTnam]]
  for (i in ss) {
    sdnam <- paste(dname,i,"sims",sep="/")
    fnam <- paste(mname,"RData",sep=".")
    if (!(fnam %in% dir(sdnam)))
      stop(paste("No simulation for subject",i))
    load(paste(sdnam,fnam,sep="/"))
    dat[dat[[D$S]]==i,c("lm","lc",rtsimnam,"qp.sim")] <- sim
  }
  attr(dat,"dname") <- dname
  attr(dat,"model") <- mname
  dat
}

# fl=rc; subject=NA
# model=topnam; mhname=NA;trans=T; reduce=F 
getpar.rc <- function(fl,model=NA,mhname=NA,trans=T,reduce=T,subject=NA) {
  if (is.na(mhname)) mhname <- names(fl)[2]
  if (!(mhname %in% names(fl)[-1]))
    stop("Model Hierarchy name not present in fit list")
  mnams <- names(fl[[mhname]]$MT[[1]])
  if (is.na(model)) model <- length(mnams)
  if (is.numeric(model)) model <- mnams[model]
  if (!(model %in% mnams))
    stop("Model name not present in fit list")
  out <- getcellpars(fl,mhname=mhname,reduce=reduce,trans=trans,type=model) 
  if (is.na(subject)) out else 
    lapply(out,function(x){x[x[,1]==subject,]})
}


####################################################################################

# fl=g1st0;mhname=NA; type="bic"; trans=T;showS=T;showbestmodels=2;digits=3
# Sshow=NA; showtest=F;show.BMA=F;show.TOP=F; excludeS=NA; output=""
# topname="B~b*lR & A~b & v~b*C*wnw & sv~1 & ter~1 & pc~1 & st0~1" 

# TOP="B~E*CT & A~1 & v~BW*C & sv~1 & ter~1 & pc~1"
# TOP="B~lR*E*CT & A~1 & v~BW*C & sv~1 & ter~1 & pc~1"


combine.fitlist <- function(fl1,fl2,exclude.subjects1=NA,exclude.subjects2=NA) {
  out <- fl1
  snam <- attributes(fl1$dat)$D[[1]]$S
  sn1 <- names(fl1[[2]]$M)
  if (!is.na(exclude.subjects1))
    sn1 <- sn1[!(sn1 %in% exclude.subjects1)]
  sn2 <- names(fl2[[2]]$M)
  if (!is.na(exclude.subjects2))
    sn2 <- sn2[!(sn2 %in% exclude.subjects2)]
  out$dat <- rbind(fl1$dat[fl1$dat[[snam]] %in% sn1,],
                   fl2$dat[fl2$dat[[snam]] %in% sn2,])
  out$dat[[snam]] <- factor(as.character(out$dat[[snam]]))
  newlist <- vector(mode="list",length=length(sn1)+length(sn2))
  names(newlist) <- c(sn1,sn2)
  out[[2]]$M <- newlist
  out[[2]]$MT <- newlist
  for (i in sn1) {
    out[[2]]$M[[i]] <- fl1[[2]]$M[[i]]
    out[[2]]$MT[[i]] <- fl1[[2]]$MT[[i]]
  }
  for (i in sn2) {
    out[[2]]$M[[i]] <- fl2[[2]]$M[[i]]
    out[[2]]$MT[[i]] <- fl2[[2]]$MT[[i]]
  }
  fn <- matrix(names(out[[2]]$fits),nrow=2)
  for (i in 1:3) {
    isin1 <- fl1[[2]]$fits[[fn[2,i]]][,1] %in% sn1
    isin2 <- fl2[[2]]$fits[[fn[2,i]]][,1] %in% sn2
    out[[2]]$fits[[fn[2,i]]] <-  rbind(out[[2]]$fits[[fn[2,i]]][isin1,],
      fl2[[2]]$fits[[fn[2,i]]][isin2,])  
    if (i != 2) out[[2]]$fits[[fn[1,i]]] <- rbind(
      out[[2]]$fits[[fn[1,i]]][isin1,],fl2[[2]]$fits[[fn[1,i]]][isin2,]) else
      out[[2]]$fits[[fn[1,i]]] <- c(
        out[[2]]$fits[[fn[1,i]]][isin1],fl2[[2]]$fits[[fn[1,i]]][isin2])    
  }  
  out
}


# mhname=NA;type="aic";trans=T;showS=F;Sshow=NA;showbestmodels=1;do.pars=T;show.BMA=F
# show.TOP=F;top.summary=TRUE;showtest=F;topname=NA;only.top=F;TOP=1;digits=0;
# excludeS=NA;output=""
# 
# fl=lba
#
# digits=3; show.BMA=T
# fl=lba; output="pars"

summary.rc=function(fl, # fitlist object
  mhname=NA,            # Name of model hierarchy in fitlist, default is second name
  type="aic",           # bic or aic
  trans=TRUE,           # natural parameterizaiton?
  showS=FALSE,          # Show results for subjects
  Sshow=NA,             # Subset of subjects to show (!is.na => dont show average)
  showbestmodels=1,     # How many of best models to show BIC etc.
  do.pars=TRUE,         # Get parameter estimates?
  show.BMA=FALSE,       # Show BMA parameters? Also switches output="pars"
  show.TOP=FALSE,       # Show TOP model paraemters?
  top.summary=TRUE,     # Print top model summary?
  showtest=FALSE,       # Show ANOVA tests on parameters (if show.BMA or show.TOP)
  topname=NA,           # reduce heirarchy with this as new top model
  only.top=FALSE,       # reduce to ONLY the topname model (to pick out single model results)
  TOP=1,                # model name or number to treat as top AFTER reduction
  digits=0,             # Accuracy of parameter printout
  excludeS=NA,          # Subjects to exclude from analysis
  output="") {          # "pars" => return parameters (BMA and TOP)
                        # "stats" => dev, bic, aic npar convergence
  
  reduce_fl <- function(flhm,mapsin,show.BMA,reduce=T,trans=F) {
    mlist <- vector(mode="list",length=length(mapsin))
    names(mlist) <- mapsin
    for (i in names(flhm$MT)) {
      MTi <- mlist
      fitsi <- mlist
      for (j in names(mlist)) {
        MTi[[j]] <- flhm$MT[[i]][[j]]
        fitsi[[j]] <- flhm$fits[[i]][[j]]      
      }
      flhm$MT[[i]] <- MTi
    }
    ok <- flhm$fits$pars.names[,"hmname"] %in% mapsin
    flhm$fits$pars.names <- flhm$fits$pars.names[ok,]
    flhm$fits$pars <- flhm$fits$pars[ok,]
    ok <- flhm$fits$pests.names[,"hmname"] %in% mapsin
    flhm$fits$pests.names <- flhm$fits$pests.names[ok,]
    flhm$fits$pests <- flhm$fits$pests[ok]
    ok <- flhm$fits$stats.names[,"hmname"] %in% mapsin
    flhm$fits$stats.names <- flhm$fits$stats.names[ok,]
    flhm$fits$stats <- flhm$fits$stats[ok,]
    if (show.BMA) flhm$fits$cellpars <- 
      getbma(flhm,doall=F,reduce=reduce,trans=trans)
    flhm
  }  
    
  if (any(is.na(mhname))) mhname <- names(fl)[2]
  if (!(type %in% c("aic","bic")))
    stop("Only analyses of type \"bic\" and \"aic\" supported") 
  # reduce if new top specified
  if (length(topname)>1) stop("Only one top model at a time!")
  if (!is.na(topname)) {
    mnams <- names(fl[[mhname]]$MT[[1]])
    if (!(topname %in% mnams)) stop("Top model not in model tree")
    maps <- lapply(fl[[mhname]]$MT[[1]],function(x){x})
    if (any(is.na(unlist(maps))))
      stop("The full model tree is not present")
    mapsin <- topname
    if (!only.top) {
      mapsi <- maps[[topname]]
      mapsin <- c(mapsin,mapsi)
      repeat {
        mapsi1 <- vector(length=0)
        for (i in mapsi) mapsi1 <- c(mapsi1,maps[[i]])
        mapsi1 <- unique(mapsi1)
        mapsin <- c(mapsin,mapsi1)
        mapsi <- mapsi1
        if (length(mapsi)==0) break
      }
    }
    fl[[mhname]] <- reduce_fl(fl[[mhname]],mapsin,show.BMA)
  }
  mnams <- names(fl[[mhname]]$MT[[1]])
  # collect results
  stats <- cbind.data.frame(s=fl[[mhname]]$fits$stats.names[,1],
    hmname=factor(fl[[mhname]]$fits$stats.names[,"hmname"],mnams),
    fl[[mhname]]$fits$stats)
  ptypes <- paste(fl[[mhname]]$fits$pests.names[,"pname"],
    fl[[mhname]]$fits$pests.names[,"ptype"],sep=".")
  ptype <- ptypes[fl[[mhname]]$fits$pests.names[,"hmname"]==
    levels(stats$hmname)[1] & fl[[mhname]]$fits$pests.names[,"s"]==
    levels(stats$s)[1]]
  pests <- cbind.data.frame(s=fl[[mhname]]$fits$pests.names[,1],
    hmname=factor(fl[[mhname]]$fits$pests.names[,"hmname"],mnams),
    pname=fl[[mhname]]$fits$pests.names[,"pname"],
    ptype=factor(ptypes,ptype),pests=fl[[mhname]]$fits$pests)     
  # Get parameters for model TOP
  if (is.numeric(TOP)) TOP <- mnams[TOP]
  if (!any(TOP==mnams))
    stop("TOP not in model tree")
  if (do.pars) cellpars <- 
    getcellpars(fl,mhname=mhname,reduce=T,trans=trans,type=TOP)
  
  if (!any(is.na(excludeS))) {
    stats$s <- as.character(stats$s)
    pests$s <- as.character(pests$s)
    stats <- stats[!(stats$s %in% excludeS),]
    pests <- stats[!(pests$s %in% excludeS),]
    stats$s <- factor(stats$s)
    pests$s <- factor(pests$s)
    if (do.pars) for (i in names(cellpars))
      cellpars[[i]] <- cellpars[[i]][!(cellpars[[i]]$s %in% excludeS),]  
  }
  
  # get model weights
  models <- levels(stats$hmname)
  weights <- tapply(stats[[type]],stats$s,function(x){
    exp(-(x-min(x))/2)/sum(exp(-(x-min(x))/2))})   
  for (i in names(weights)) names(weights[[i]]) <- models   
  
  if ( do.pars ) {
    # Get BMA
    if (!show.BMA) tmp <- cellpars else { 
      if (!any(names(fl[[mhname]]$fits)=="cellpars"))
        stop("A BMA has not been calculated")
      tmp <- fl[[mhname]]$fits$cellpars[[type]]["TRUE",as.character(trans)][[1]] 
    }
    for ( i in names(cellpars) ) {
      cellpars[[i]] <- cbind.data.frame(cellpars[[i]][,-dim(cellpars[[i]])[2]],
        TOP=cellpars[[i]]$y,BMA=tmp[[i]]$y)
      names(cellpars[[i]])[1]="s"
    }
    pnames <- names(cellpars) 
  }
  if (is.numeric(Sshow)) Sshow <- names(weights)[Sshow]
  # process each participant
  for (i in names(weights)) if (showS | (!is.na(Sshow) && Sshow==i)) {
    # get best model for subject i
    stat <- cbind(stats[stats$s==i,],pmodel=weights[[i]])
    row.names(stat) <- stat$hmname
    stati <- stat[,-1]
    stati[,7] <- round(stati[,7],4)
    cat(paste(rep("#",80),collapse=""))
    cat(paste("\nPARTICIPANT",i,"\n"))
    if (showbestmodels>0) {
      cat(paste("\nBest",showbestmodels,type,"models\n"))
      print(round(stati[order(weights[[i]],decreasing=T),][1:showbestmodels,-c(1,6)],digits=digits))
      cat("\n")
    }
    if (do.pars) if (show.BMA | show.TOP) for (j in names(cellpars)) {
      cols <- 1:dim(cellpars[[j]])[2]
      if (show.BMA & !show.TOP) cols <- cols[names(cellpars[[j]])!="TOP"]
      if (!show.BMA & show.TOP) cols <- cols[names(cellpars[[j]])!="BMA"]      
      cat(paste("Estimates for parameter",j,"\n"))
      print(round(cellpars[[j]][cellpars[[j]]$s==i,cols],digits=digits))
      cat("\n")
    }
  } # end of each subject analyses  
  if (showS) {
  	cat(paste(rep("#",80),collapse=""))
    cat("\nOVERALL RESULTS \n")
  }
  if ( is.na(Sshow) ) {
  	if (top.summary) cat("TOP MODEL\n")
    # Best models
    sumstats <- cbind.data.frame(
      dev=tapply(stats[["dev"]],stats$hmname,sum),
      npar=as.vector(tapply(stats[["npar"]],stats$hmname,sum)))
    sumstats$dev=as.vector(sumstats$dev)
      groupn <- dim(fl$dat)[1]
    sumstats$bic <- sumstats$dev+sumstats$npar*log(groupn)
    sumstats$aic <- sumstats$dev+2*sumstats$npar
    weight <- exp(-(sumstats[[type]]-min(sumstats[[type]]))/2)/
      sum(exp(-(sumstats[[type]]-min(sumstats[[type]]))/2))
    sumstats <- cbind(sumstats,pmodel=round(weight,4))
    if (top.summary) print(round(sumstats[1,],digits=digits))
    if (showbestmodels>0) {
      cat(paste("\n\nBest",showbestmodels,type,"models\n"))
      print(round(sumstats[order(sumstats[[type]]),][1:showbestmodels,],digits=digits))  
    }
    if ( do.pars ) {
      if ( show.TOP | show.BMA )
        cat("\nAVERAGE DESIGN PARAMETER ESTIMATES AND TESTS\n")
      for ( i in names(cellpars) ) {
        ncol <- dim(cellpars[[i]])[2]
        if ( show.BMA ) {
          cat(paste("\nAVERAGE BMA ESTIMATE(S) FOR",i,"\n"))
          if ( ncol==3 ) 
            print(round(mean(cellpars[[i]][,3]),digits=digits)) else
            print(round(tapply(cellpars[[i]][,ncol-1],
              cellpars[[i]][,-c(1,ncol-1,ncol)],mean),digits=digits))
          if (showtest & ncol>3)
            try(wsAnova(cellpars[[i]][,-ncol]))
        }
        if ( show.TOP ) {
          cat(paste("\nAVERAGE TOP ESTIMATE(S) FOR",i,"\n"))
          if ( ncol==3 ) 
            print(round(mean(cellpars[[i]][,3]),digits=digits)) else
            print(round(tapply(cellpars[[i]][,ncol],
              cellpars[[i]][,-c(1,ncol-1,ncol)],mean),digits=digits))
          if (showtest & ncol>3)
            try(wsAnova(cellpars[[i]][,-(ncol-1)]))     
          cat("\n")
        }
      }
    }
  }
  if ( do.pars ) {
    # default is to savepar TOP, renamed "y"
    if ( !show.TOP & show.BMA ) for (i in names(cellpars)) {
      cellpars[[i]] <- 
        cellpars[[i]][,names(cellpars[[i]])[names(cellpars[[i]])!="TOP"]]
    } else for (i in names(cellpars)) 
      cellpars[[i]] <- cellpars[[i]][,names(cellpars[[i]])[names(cellpars[[i]])!="BMA"]]
    for (i in names(cellpars)) names(cellpars[[i]])[dim(cellpars[[i]])[2]] <- "y"
    if (output=="pars") return(cellpars)
  }
  if (output=="stats") stats
}

# measure="accuracy";stat="mean"
# anovas=F; factors=NA;probs=.5;digits=3;sfac=0; show.unbalanced=TRUE; 
# save=F; sfac=0
# 
# rc=seq; anovas=TRUE

table.rc <- function(rc,
  # raw or probit accuracy, all/correct/error rt
  measure = c("accuracy","probit","rt","correct","error")[1],  
  stat = c("mean","quantile","length","sd","var")[1], # name of any statistic
  factors=NA,                                         # subset of factor names
  probs=.5,                                           # quantile probabilit
  show.unbalanced=TRUE,                               # na.rm=T subject averages shown
  digits=3,                                           # table rounding
  anovas=FALSE,                                       # do within subjects anova?
  save=FALSE,                                         # save data frame
  sfac=0) {  # 0=average, 1 fastest ... length(factors)+1 slowest
  
  check <- function(wsF,allF) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
  }  
  
  getstat <- function(dat,measure,stat,probs,wsF,Di,verbose=T) {
    if (!is.null(attr(dat,"qp")) & measure[1]=="accuracy") 
      dat <- dat[!is.finite(dat[,Di$RT]) & as.logical(dat[,Di$SC[1]]),] else
      dat <- dat[is.finite(dat[,Di$RT]),]
    if (stat[1]=="quantile") {
      if (verbose) cat(paste(" at p =",probs[1],"\n"))
      if (!is.null(qp)) {
        tmp <- dat[dat$qp==probs[1],]
        data <- tapply(tmp[,Di$RT],tmp[,c(Di$S,wsF)],function(x){mean(x)})
      } else data <- tapply(dat[,Di$RT],dat[,c(Di$S,wsF)],stat,probs=probs[1],type=2)
    } else {
      if (verbose) cat("\n")
      if (measure[1]=="probit") {
        data <- qnorm((tapply(dat[,Di$RT],dat[,c(Di$S,wsF)],sum)+.5)/ 
          (tapply(dat[,Di$RT],dat[,c(Di$S,wsF)],length)+1))      
       } else if (!is.null(attr(dat,"qp")) & measure[1]=="accuracy")
         data <- tapply(dat$qp,dat[,c(Di$S,wsF)],function(x){mean(x)}) else
         data <- tapply(dat[,Di$RT],dat[,c(Di$S,wsF)],stat)       
    }
    data
  }
  
  if (is.data.frame(rc)) dat <- rc else dat <- rc[[1]]
  qp <- attr(dat,"qp")
  
  Di <- attr(dat,"D")[[1]]
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
#   check(wsF,Di$F)
  measures <- c("rt","correct","error","accuracy","probit")
  if (is.na(pmatch(measure[1],measures)))
    stop("Not a supported measure (rt, correct, error, accuracy)")
  if (!is.null(qp)) {
    if (measure %in% c("probit","rt"))
      stop("probit and rt measure not available for quantile data")
    if (measure=="accuracy") {
      if (stat!="mean")
        stop("For accuracy only mean can be used with quantile data")
    } else {
      if (stat != "quantile")
        stop("Only quantile stat can be applied to quantile data")
      if (!all(probs %in% qp))
        stop(paste("Quantile probs must be in set",paste(qp,collapse=","),"\n"))
    }
  }
  dat <- switch(pmatch(measure[1],measures),
    {cat("\nAll RT: "); dat},
    {cat("\nCorrect RT: "); dat[as.logical(dat[,Di$SC[1]]),]},
    {cat("\nError RT: "); dat[!as.logical(dat[,Di$SC[1]]),]},
    {cat("\nAccuracy: "); dat[is.finite(dat[,Di$RT]),Di$RT] <- 
      as.numeric(as.logical(dat[is.finite(dat[,Di$RT]),Di$SC[1]])); dat},
    {cat("\nProbit accuracy: "); dat[,Di$RT] <- as.logical(dat[,Di$SC[1]]); dat}
  )
  cat(stat)
  data <- getstat(dat,measure,stat,probs,wsF,Di)
  snams <- dimnames(data)[[1]]
  bad <- apply(data,1,function(x){any(is.na(x))})
  if (any(bad) & show.unbalanced) {
    cat("Average including incomplete participants\n")
    print(round(apply(data,2:length(dim(data)),mean,na.rm=T),digits))
  }
  if (all(bad))
    stop("All participants have missing data cells")
  if (any(bad)) {
    cat(paste("Excluding participants",paste(snams[bad],collapse=","),
      "(",sum(!bad),"left )\n"))  
    tmp <- dat[dat[,Di$S] %in% snams[!bad],]
    tmp[,Di$S] <- factor(as.character(tmp[,Di$S]))
    data <- getstat(tmp,measure,stat,probs,wsF,Di,F) 
  }
  if (anovas) wsAnova(arr2df(data))
  cat("\n")
  if (sfac>0) {
    fodr <- order(c(sfac+.5,2:(length(wsF)+1)))
    print(round(aperm(data,fodr),digits))
  } else print(round(apply(data,2:length(dim(data)),mean),digits))
  if (save) arr2df(data)
}


# rc=myer1; model="a~E & v~E*S & sv~E*S & Z~E & SZ~E & t0~E & st~1 & pc~E"
# factors=NA; probs=c(.1,.3,.5,.7,.9); gsq=F; by.subjects=T
# include.subjects=NA; exclude.subjects=NA;exclude.incomplete=F
# 
qgof.rc <- function(rc,
    by.subjects=F,                            # output per subject
    correct=0,                                # fix zero p
    factors=NA,                               # subset of factor names
    probs=list(c=c(.1,.5,.9),e=c(.1,.5,.9)),  # quantile probability, replicated if vector 
    include.subjects=NA,                      # subjects to include in analysis NA => all
    exclude.subjects=NA,
    model=NA) {                               # specify to get bic and aic               
  # If no factors values specified sets to all factors
    
  check <- function(wsF,allF,data.fit.name) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
  }      
    
  if (!is.list(probs)) probs <- list(c=probs,e=probs)
  if (any(unlist(lapply(probs,length))<2))
    stop("Need at least two quantile probabiliteis")
  if (!all(names(probs) %in% c("e","c")))
    stop("Names of probs must be c and e")
  probs <- lapply(probs,sort)
  dat <- rc[[1]]
  qp <- attr(dat,"qp")
  D <- attr(dat,"D")
  Di <- D[[1]]
  # Get factors
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
  check(wsF,Di$F,data.fit.name)
  dwsFs <- c(wsF,Di$S)
  if (!is.null(qp)) if ( !( all(probs$c %in% qp) & all(probs$e %in% qp) ) )
      stop(paste("Quantile probs must be in set",paste(qp,collapse=","),"\n"))
  probs$c <- c(probs$c,1)
  probs$e <- c(probs$e,1)  
  snams <- levels(dat[[Di$S]])
  if (!any(is.na(exclude.subjects))) 
    snams <- snams[!(snams %in% exclude.subjects)]
  dat <- dat[dat[,Di$S] %in% snams,]
  if (!any(is.na(include.subjects))) {
    if (!all(include.subjects %in% snams))
      stop("Some specifed subjects not in set of subjects")
    cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
    dat <- dat[dat[,Di$S] %in% include.subjects,]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    snams <- levels(dat[,Di$S])
  }  
  dat <- dat[order(dat$s),] # might be sorted in numeric not name order
  # sum to 1, last is probability of response, 
  # convert to probability suming to one over rcell
  dat[["qp.sim"]] <- unlist(
    tapply(dat[["qp.sim"]],dat[,c("cell",Di$S)],function(x){
      c(x[-length(x)],1-sum(x[-length(x)]))})
    )
  if ( any(abs(tapply(dat[["qp.sim"]],dat[,c("cell",Di$S)],sum)-1)>1e-6) )
    Warning("Possible error in qp.sim alignment!")  
  
  cns <- vector(length=length(probs$c),mode="list")
  cNs <- cps <- cns
  for (i in 1:length(probs$c)) {
    # get quantiles
    if (is.null(qp)) {
      data=arr2df(tapply(dat[as.logical(dat[,Di$SC[1]]),Di$RT],
        dat[as.logical(dat[,Di$SC[1]]),dwsFs],quantile,probs=probs$c[i])) 
      row.names(data) <- apply(data[,dwsFs],1,function(x){
        paste(as.character(unlist(x)),collapse="")})
      # RT less than or equal to quantile
      ok <- as.logical(dat[,Di$SC[1]]) & dat[,Di$RT] <= 
        data[apply(dat[,dwsFs],1,function(x){
        paste(as.character(unlist(x)),collapse="")}),"y"]
      cns[[i]] <- tapply(ok,dat[,dwsFs],sum)
      cps[[i]] <- tapply(dat[["qp.sim"]][ok],dat[ok,dwsFs],sum)
      cps[[i]][is.na(cps[[i]])] <- 0
    } else { # data is quantiles
      ok <- as.logical(dat[,Di$SC[1]])
      if (probs$c[i]!=1) ok <- ok & (dat$qp <= probs$c[i]) & is.finite(dat[,Di$RT])   
      cns[[i]] <- tapply(dat$qn[ok],dat[ok,dwsFs],sum)
      cps[[i]] <- tapply(dat$qp.sim[ok],dat[ok,dwsFs],sum) 
    }
  }
  ens <- vector(length=length(probs$e),mode="list")
  eNs <- eps <- ens
  for (i in 1:length(probs$e)) {
    if (is.null(qp)) {
      # get quantiles
      data=arr2df(tapply(dat[!as.logical(dat[,Di$SC[1]]),Di$RT],
                         dat[!as.logical(dat[,Di$SC[1]]),dwsFs],quantile,probs=probs$e[i]))  
      row.names(data) <- apply(data[,dwsFs],1,function(x){
        paste(as.character(unlist(x)),collapse="")})
      # RT less than or equal to quantile
      ok <- !as.logical(dat[,Di$SC[1]]) & dat[,Di$RT] <= 
        data[apply(dat[,dwsFs],1,function(x){
          paste(as.character(unlist(x)),collapse="")}),"y"]
      ens[[i]] <- tapply(ok,dat[,dwsFs],sum)
      eps[[i]] <- tapply(dat[["qp.sim"]][ok],dat[ok,dwsFs],sum)
      eps[[i]][is.na(eps[[i]])] <- 0
    } else { # data is quantiles
      ok <-  !as.logical(dat[,Di$SC[1]])
      if (probs$e[i]!=1) ok <- ok & (dat$qp <= probs$c[i]) & is.finite(dat[,Di$RT])  
      ens[[i]] <- tapply(dat$qn[ok],dat[ok,dwsFs],sum)
      eps[[i]] <- tapply(dat$qp.sim[ok],dat[ok,dwsFs],sum) 
    }
  }
  # convert model p from cumulative and to ns
  for (i in length(probs$c):2) 
    cps[[i]] <- cps[[i]]-cps[[i-1]]
  for (i in length(probs$e):2) 
    eps[[i]] <- eps[[i]]-eps[[i-1]]
  # zero probability correction, NOT RENORMALIZED
  cps <- lapply(cps,function(x){x[x==0]<-correct/100;x})
  eps <- lapply(eps,function(x){x[x==0]<-correct/100;x})
  for (i in length(probs$c):2) 
    cNs[[i]] <- cns[[length(probs$c)]]*cps[[i]]
  for (i in length(probs$e):2) 
    eNs[[i]] <- ens[[length(probs$e)]]*eps[[i]]
  cNs[[1]] <- cns[[length(probs$c)]]*cps[[1]]
  eNs[[1]] <- ens[[length(probs$e)]]*eps[[1]]
  # convert ns from cumulative
  for (i in length(probs$c):2)
    cns[[i]] <- cns[[i]]-cns[[i-1]]
  for (i in length(probs$e):2)
    ens[[i]] <- ens[[i]]-ens[[i-1]]
  if (!by.subjects) {
    cns <- lapply(cns,function(x){apply(x,1:(length(dim(x))-1),sum)})
    ens <- lapply(ens,function(x){apply(x,1:(length(dim(x))-1),sum)})
    cNs <- lapply(cNs,function(x){apply(x,1:(length(dim(x))-1),sum)})
    eNs <- lapply(eNs,function(x){apply(x,1:(length(dim(x))-1),sum)})
    cps <- lapply(cps,function(x){apply(x,1:(length(dim(x))-1),mean)})
    eps <- lapply(eps,function(x){apply(x,1:(length(dim(x))-1),mean)})
  }
  cg2 <- cc2 <- cd <- cns
  eg2 <- ec2 <- ed <- ens
  for (i in 1:length(cps)) {
    cg2[[i]] <- 2*log(cns[[i]]/cNs[[i]])*cns[[i]]      
    cc2[[i]] <- ((cns[[i]]-cNs[[i]])^2)/cNs[[i]]
    cd[[i]] <- -2*log(cns[[i]])*cps[[i]]      
  }
  for (i in 1:length(eps)) {
    eg2[[i]] <- 2*log(ens[[i]]/eNs[[i]])*ens[[i]]      
    ec2[[i]] <- ((ens[[i]]-eNs[[i]])^2)/eNs[[i]]
    ed[[i]] <- -2*log(ens[[i]])*eps[[i]]      
  }
  if (by.subjects) {
    for (i in 2:length(cg2)) cg2[[1]] <- cg2[[1]] + cg2[[i]]; cg2 <- cg2[[1]]
    for (i in 2:length(eg2)) eg2[[1]] <- eg2[[1]] + eg2[[i]]; eg2 <- eg2[[1]]
    g2 <- eg2+cg2
    for (i in 2:length(cc2)) cc2[[1]] <- cc2[[1]] + cc2[[i]]; cc2 <- cc2[[1]]
    for (i in 2:length(ec2)) ec2[[1]] <- ec2[[1]] + ec2[[i]]; ec2 <- ec2[[1]]
    c2 <- ec2+cc2
    for (i in 2:length(cd)) cd[[1]] <- cd[[1]] + cd[[i]]; cd <- cd[[1]]
    for (i in 2:length(ed)) ed[[1]] <- ed[[1]] + ed[[i]]; ed <- ed[[1]]
    d <- ed+cd
    ncell <- length(cns)*apply(cns[[1]],length(dim(cns[[1]])),function(x){prod(dim(x))})+
      length(ens)*apply(ens[[1]],length(dim(ens[[1]])),function(x){prod(dim(x))})
    out <- cbind(c2=apply(c2,length(dim(c2)),sum),
      g2=apply(g2,length(dim(g2)),sum),
      d=apply(d,length(dim(d)),sum),
      ncell=ncell
    )
  } else 
    out <- c(c2=sum(unlist(cc2))+sum(unlist(ec2)),
      g2=sum(unlist(cg2))+sum(unlist(eg2)),
      d=sum(unlist(cd))+sum(unlist(ed)),
      ncell=length(cns)*prod(dim(cns[[1]]))+length(ens)*prod(dim(ens[[1]]))
    )
  if (!is.na(model)) {
    if (is.null(qp)) n <- table(dat[,Di$S]) else n <- tapply(dat$qn,dat[,Di$S],sum)
    is.mod <- rc[[2]]$fits$stats.names[,"hmname"]==model
    nps=cbind.data.frame(npar=as.numeric(rc[[2]]$fits$stats[,"npar"][is.mod]))
    row.names(nps) <- rc[[2]]$fits$stats.names[,"s"][is.mod]
    nps$n <- n[row.names(nps)]
    if (by.subjects) {
      nps <- nps[row.names(out),]
      out <- cbind(out,p=nps$npar,
        aic.g2=out[,"g2"]+2*nps$npar,bic.g2=out[,"g2"]+nps$npar*log(sum(nps$n)),
        aic.d=out[,"d"]+2*nps$npar,bic.d=out[,"d"]+nps$npar*log(sum(nps$n)))   
    } else {
      np=sum(nps$npar)
      out <- c(out,p=np,
        aic=out["g2"]+2*np,bic=out["g2"]+np*log(sum(nps$n)),
        aic=out["d"]+2*np,bic=out["d"]+np*log(sum(nps$n)))
    }
  }
  out
}



####################################################################################
# rawdat <- rc$dat
# boot.se <- function(rawdat,bars,facs,measure,stat,probs,n.boot=100,boot.s=TRUE) {
#   facs <- c(facs[length(facs)],facs[-length(facs)])
#   Di <- attributes(rawdat)$D[[1]]
#   snam <- Di$S
#   rtnam <- Di$RT
#   fnams <- Di$F
#   cnam <- Di$SC[1]
#   snams <- levels(rawdat[[snam]])
#   ns <- length(snams) 
#   old.rows <- 1:dim(rawdat)[1]
#   out <- vector(mode="list",length=n.boot)
#   for (n in 1:n.boot) {
#     if (boot.s) {
#       snams <- sample(snams,length(snams),TRUE)
#       srows <- vector(mode="list",length=ns); names(srows)=snams
#       rows=numeric(0)
#       for (i in snams) 
#         rows <- c(rows,old.rows[rawdat[[snam]]==i])
#       dat <- rawdat[rows,]
#       dat[,snam] <- factor(rep.int(1:ns,as.vector(table(rawdat$s)[snams])))    
#     } else dat <- rawdat
#     new.rows <- 1:dim(dat)[1]
#     for (i in levels(dat[[snam]])) for (j in levels(dat$rcell)) {
#       rowij <- new.rows[dat[[snam]]==i & dat$rcell==j & is.finite(dat[[rtnam]])] 
#       dat[rowij,c(cnam,rtnam)] <- dat[sample(rowij,length(rowij),TRUE),c(cnam,rtnam)]
#     }
#     dat <- switch(pmatch(measure[1],measures),
#       dat,
#       dat[as.logical(dat[,cnam]),],
#       dat[!as.logical(dat[,cnam]),],
#       fixholes(dat,Di,measure),
#       fixholes(dat,Di,measure),
#       fixstop(dat,Di,outc=fixholes(dat,Di,"accuracy",plot.model),plot.model)
#     )      
#     out[[n]]  <- getstat(dat,measure,stat,probs,facs,Di)
#   }  
#   save(out,file="out.RData")  
# }

# measure="inaccuracy"; stat="mean"; factors=NA
# probs=.5; xlim=NA; ylim=NA; xlab=NA
# plot.subjects=F; bars="se"; data.name="data"
# include.subjects=NA; exclude.subjects=NA; exclude.incomplete=F;
# plot.model=T;  model.name="fit"; data.fit.name="Type"; 
# plot.lines=T; cex=1; save.dat=F; do.plot=T; layout=NULL;
# 
# rc=ejldt1pc;measure="inaccuracy";plot.model=T;plot.subjects=F;bars=95
# xlab="Emphasis";model.name=model.name;factors="E";bar.type="pbmodel"

plot.rc <- function(rc,
  # raw accuracy, all/correct/error rt
  measure = c("accuracy","inaccuracy","rt","correct","error","stop")[1],  
  stat = c("mean","quantile","sd","var")[1], # name of any statistic
  factors=NA,                                         # subset of factor names
  probs=.5,                                           # quantile probability
  xlim=NA,ylim=NA,                                    # plot limits
  xlab=NA,
  plot.subjects=FALSE,                                # subjects as panels
  bars="se",                                          # if NOT subjects do "se" or % CI (if numeric)
  bar.type=c("within","between","data","model","pbmodel")[1],
    # within=Morey 2008, between=sd/sqrt(n), data=resample data (not fully implemented) 
    # model=simulate model or pbmodel=parameteric bootstrap
  include.subjects=NA,                                # subjects to include in analysis NA => all
  exclude.subjects=NA,
  exclude.incomplete=FALSE,                           # dump entire subject where any cell NA
  plot.model=FALSE,                                   # data and model as lines
  model.name="fit",
  data.name="data",
  data.fit.name="Type",                               # data/model factor name
  plot.lines=TRUE,                                    # first factor as lines (if !plot.model)
  cex=1,
  save.dat=FALSE,
  do.plot=T,
  layout=NULL) {                                      # xyplot key layout parameter
  # if no factors values specified sets to all factors
  # If plot.lines, plots first entry in factors as lines except if only one, 
  # then plots as x, otherwise second factor is x, further factors panels.
  # If not plot.lines first factor is x, subsequent factors are panels
  # if plot.subjects add subjects as last factor (i.e., panels)
  # if plot.model adds data vs. fit as first factor (i.e., lines) OVERRIDES plot.lines
  # layout length=2 as c(0,n) => panels/pages, c(x,y) x columns, y rows, c(x,y,z) z=num pages
  
  panel.se <- function(x, y, groups, subscripts, ...) {
    if (any(groups %in% c("lo","hi"))) { # doing ses
      do.se <- T; is.datafit <- !(groups=="lo" | groups=="hi")
      if (any(groups==model.name)) {
        is.lohi <- c(1:length(x))>length(x)/2
        is.data <- c(1:length(x))<=length(x)/4
      } else {
        is.lohi <- c(1:length(x))>length(x)/3
        is.data <- c(1:length(x))<=length(x)/3        
      }
      y.lohi <- y[is.lohi]
      if (bar.type!="model") 
        y.0 <- c(y[is.data],y[is.data]) else
        y.0 <- c(y[!is.data & !is.lohi],y[!is.data & !is.lohi])
      x.0 <- c(x[is.data],x[is.data]) 
      x <- x[is.datafit]; y <- y[is.datafit]; groups <- groups[is.datafit]
    } else do.se <- F
    panel.xyplot(x,y,groups=groups,subscripts=subscripts, ...)
    #      panel.text(4,8,paste(round(y.lohi,1),collapse="|"))
    if (do.se) for (i in 1:length(y.0)) if (bar.type!="model")
      panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1) else
      panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1,lty=2)
  }  
  
  check <- function(wsF,allF,data.fit.name) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
    if (any(wsF==data.fit.name))
      stop(paste("data.fit.name (",data.fit.name,,
                 "and factor name are the same, change the former"))
  }  
  
  # verbose=F;sim=T; wsF=wsFs
  getstat <- function(dat,measure,stat,probs,wsF,Di,verbose=F,sim=F) {
    pc <-c("accuracy","inaccuracy","stop")
    if (!(measure %in% pc)) # && stat != "quantile") 
      dat <- dat[is.finite(dat[,Di$RT]),]
    if (sim) {
      rtnam=paste(Di$RT,"sim",sep=".") 
      qpnam="qp.sim"
      if (!any(names(dat)==rtnam))
        stop("Simulation results have not been added to rc object")
    } else {
      rtnam <- Di$RT
      qpnam="qp"
    }
    if (stat[1]=="quantile") {
      if (verbose) cat(paste(" at p =",probs[1],"\n"))
      if (!is.null(qp)) {
        tmp <- dat[dat$qp==probs[1] & is.finite(dat[[Di$RT]]),]
        if (any(is.na(tmp[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(tmp[,rtnam],tmp[,wsF],function(x){mean(x)})
      } else {
        if (any(is.na(dat[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(dat[,rtnam],dat[,wsF],stat,probs=probs[1])
      }
    } else {
      if (verbose) cat("\n")
      if (measure[1] %in% pc) {
        if (any(is.na(dat[,qpnam])))
          stop(paste("NAs present in statistic",qpnam))
        data <- tapply(dat[,qpnam],dat[,wsF],function(x){mean(x)}) 
        #        data[is.na(data)] <- 0    
      } else {
        if (any(is.na(dat[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(dat[,rtnam],dat[,wsF],stat)
      }
    }
    names(dimnames(data)) <- wsF
    if (measure[1] %in% pc) {
       data[is.na(data)] <- 0
      data*100
    } else data
  }
  
  
  fixholes <- function(dat,Di,measure="inaccuracy",plot.model=FALSE) {
    
    cells <- function(d)
      apply(matrix(unlist(lapply(d,as.character)),ncol=dim(d)[2]),1,paste,collapse=" ")
    
    s <- Di$S
    dc=dat[is.infinite(dat[,Di$RT]) & as.logical(dat[,Di$SC[1]]),]
    de=dat[is.infinite(dat[,Di$RT]) & !as.logical(dat[,Di$SC[1]]),]
    cc <- cells(dc[,c(s,"rcell")]); ce <- cells(de[,c(s,"rcell")])
    # assumes only error responses can be non-unique!
    dups <- duplicated(ce) 
    deu <- de[!dups,]   
    de <- de[dups,]
    row.names(dc) <- cc; row.names(deu) <- ce[!dups]
    ce <- ce[dups]
    if (length(ce)>0) for (i in 1:length(ce)) {
      deu[ce[i],"qp"] <- deu[ce[i],"qp"] + de[i,"qp"]
      if (plot.model) 
        deu[ce[i],"qp.sim"] <- deu[ce[i],"qp.sim"] + de[i,"qp.sim"]      
    }
    out <- deu
    missing <- cc[!(cc %in% row.names(out))]
    if ( length(missing)!=0 ) {
      out <- rbind(out,dc[missing,])
      out[missing,"qp"] <- 1-dc[missing,"qp"]
      if (plot.model) 
        out[missing,"qp.sim"] <- 1-dc[missing,"qp.sim"]
    }
    if (measure=="accuracy") {
      out[,"qp"] <- 1 - out[,"qp"]
      if (plot.model) 
        out[,"qp.sim"] <- 1 - out[,"qp.sim"]
      out[,Di$SC[1]] <- T
    } else out[,Di$SC[1]] <- F
    out
  }
  
#   outc=fixholes(dat,Di,"accuracy",plot.model)
  fixstop <- function(dat,Di,outc,plot.model) {
    out <- outc
    out$qp <- 0
    out$qp.sim <- 0
    dc <- dat[is.infinite(dat[,Di$RT]),]
    fns <- do.call(paste,out[,c(Di$S,Di$F)])
    stops <- dc[dc[,Di$R]=="stop" & is.finite(dc[,"ssd"]),]
    stopfns <- do.call(paste,stops[,c(Di$S,Di$F)])
    is.in <- fns %in% stopfns
    out[is.in,"qp"] <- stops$qp
    if (plot.model) {
      out[is.in,"qp.sim"] <- stops$qp.sim
      no.stop <- do.call(paste,outc[!is.in,c(Di$S,"rcell")])
      indx <- do.call(paste,dc[,c(Di$S,"rcell")])
      for (i in no.stop) 
        out[i,"qp.sim"] <- 1-sum(dc[indx %in% i,"qp.sim"])
    }
    out <- out[out$SSD != "Inf",]
    out$SSD <- factor(as.character(out$SSD))
    out
  }
    
  
  model.bar <- function(rc,stat,measure,wsFs,probs=NA,pb=FALSE) {
    dname <- attr(rc$dat,"dname")
    model <- attr(rc$dat,"model")
    snams <- levels(rc$dat[[Di$S]])
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
      pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if ( measure=="accuracy" | measure=="correct" )  # target response to tabulate
      tr <- as.character(Di$D$R[as.logical(Di$D[[Di$SC[1]]])]) else
      tr <- as.character(Di$D$R[!as.logical(Di$D[[Di$SC[1]]])])
    cat("Calculating model errors for each subject ")
    for (s in 1:length(snams)) { # subjects
     if (!pb) load(paste(dname,"/",snams[s],"/sims/",model,"sims.RData",sep="")) else
        load(paste(dname,"/",snams[s],"/sims/",model,"=",model,".RData",sep=""))  
      if (s==1) av <- array(dim=c(length(sims),length(sims[[1]]),length(snams)))
      cat(".")
      if ( measure %in% c("accuracy","inaccuracy") ) { # need to implement "stop"
        for ( i in 1:length(tr) ) av[i,,s] <-
          tmp <- unlist(lapply(sims[[i]],function(x){100*mean(x$r==tr[[i]])}))
      } else { # RT
        for ( i in 1:length(tr) ) {
          if ( measure=="rt" ) 
            av[i,,s] <- unlist(lapply(sims[[i]],function(x){
              do.call(stat,list(x=x$rt,probs=probs,na.rm=TRUE))
            })) else
            av[i,,s] <- unlist(lapply(sims[[i]],function(x){
              do.call(stat,list(x=x$rt[x$r==tr[[i]]],probs=probs,na.rm=TRUE))
            }))
        }
      }
    }
    cat("\n")
    av <- apply(av,1:2,mean,na.rm=TRUE)
    Di <- attr(rc$dat,"D")[[1]]
    D <- Di$D[as.logical(Di$D[[Di$SC[1]]]),]
    apply(apply(av,2,function(x) {
      as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),1,quantile,probs=pci,na.rm=TRUE)
  }
    
  require(lattice)
  if (bar.type=="pbmodel") {
    bar.type="model"
    pb <- TRUE
  } else pb <- FALSE
  if (plot.subjects) bars <- NA
  dat <- rc[[1]]
  qp <- attr(dat,"qp")
  D <- attr(dat,"D")
  Di <- D[[1]]

  if (!any(is.na(exclude.subjects))) {
    dat <- dat[!(dat[,Di$S] %in% exclude.subjects),]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
  }  
  
  # Get factors
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
#   check(wsF,Di$F,data.fit.name)
  if (plot.subjects) {
    wsF <- c(wsF,Di$S)
    bars <- NA
  }
  if ( plot.model ) {
    wsF <- c(data.fit.name,wsF)
    plot.lines <- T
    wsFs <- wsF[-1]
  } else wsFs <- wsF
  
  if (!is.na(bars) & length(wsF)==1) plot.lines=FALSE
  bar.no.line=FALSE
  if ( !is.na(bars) & !plot.lines ) {
    plot.lines <- T
    wsF <- c(data.fit.name,wsF)
    if (!plot.model) bar.no.line=TRUE 
  }
  if (plot.subjects) dwsFs <- wsFs else dwsFs <- c(wsFs,Di$S)
  if (length(wsF)==1) plot.lines <- F 
  # plotting formula
  if (plot.lines) {
    plotform <- paste("y ~",wsF[2]) 
    if (length(wsF)>2) {
      if (length(wsF)==3) plotform <- paste(plotform,"|",wsF[3]) else
        plotform <- paste(plotform,"|",paste(wsF[-c(1,2)],collapse="+"))
    }
  } else {
    plotform <- paste("y ~",wsF[1])
    if (length(wsF)>1) {
      if (length(wsF)==2) plotform <- paste(plotform,"|",wsF[2]) else
        plotform <- paste(plotform,"|",paste(wsF[-1],collapse="+"))
    }
  }
  # Get measure
  measures <- c("rt","correct","error","accuracy","inaccuracy","stop")
  if (is.na(pmatch(measure[1],measures)))
    stop("Not a supported measure (rt, correct, error, accuracy, inaccuracy)")
  if (measure %in% c("accuracy","inaccuracy","stop")) {
    if (stat!="mean") stop("For accuracy & inaccuracy only mean can be used")
  } else if (!is.null(qp)) {
    if (measure=="rt")
      stop("rt measure not available for quantile data")
    if (stat != "quantile")
      stop("Only quantile stat can be applied to quantile data")
    if (!all(probs %in% qp))
      stop(paste("Quantile probs must be in set",paste(qp,collapse=","),"\n"))
  }
  dat <- switch(pmatch(measure[1],measures),
    {ylab <- "All RT (s):"; dat},
    {ylab <- "Correct RT (s):"; dat[as.logical(dat[,Di$SC[1]]),]},
    {ylab <- "Error RT (s):"; dat[!as.logical(dat[,Di$SC[1]]),]},
    {ylab <- "Accuracy (%):"; fixholes(dat,Di,measure,plot.model)},
    {ylab <- "Error Rate (%):"; fixholes(dat,Di,measure,plot.model)},
    {ylab <- "Stop (%):"; fixstop(dat,Di,
      outc=fixholes(dat,Di,"accuracy",plot.model),plot.model)}
  )  
  if (stat!="quantile") ylab <- paste(ylab,stat) else
    ylab <- paste(ylab,probs,stat)
  data <- getstat(dat,measure,stat,probs,dwsFs,Di)
  sdim <- c(1:length(dim(data)))[names(dimnames(data))==Di$S]
  snams <- dimnames(data)[[sdim]]
  bad <- apply(data,sdim,function(x){any(is.na(x))})
  if (all(bad))
    stop("All participants have missing data cells")
  if (exclude.incomplete && any(bad)) {
    cat(paste("Excluding participants with missing data cells",
              paste(snams[bad],collapse=","),"(",sum(!bad),"left )\n"))  
    dat <- dat[dat[,Di$S] %in% snams[!bad],]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstat(dat,measure,stat,probs,dwsFs,Di) 
  }
  if (!any(is.na(include.subjects))) {
    if (!all(include.subjects %in% snams))
      stop("Some specifed subjects not in set of subjects")
    cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
    dat <- dat[dat[,Di$S] %in% include.subjects,]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstat(dat,measure,stat,probs,dwsFs,Di) 
  }  
  if ( !is.na(bars) ) {
    if ( bar.type=="model" & !plot.model )
      stop("Set plot.model=TRUE to enable model error-bar plotting")
    if (!(any(bar.type %in% c("within","between","data","model")))) 
        stop("Parameter bar.type must be one of \"within\",\"between\", \"data\" or \"model\"")
    if ( bar.type %in% c("within","between") ) {
      ns <- apply(data,wsFs,function(x){sum(!is.na(x))})
      if (bar.type=="within") {
        m <- prod(unlist(lapply(dimnames(data)[wsFs],length)))
        se <- arr2df(apply(aperm(data,c(Di$S,wsFs)) - apply(data,Di$S,mean,na.rm=T),
                       wsFs,sd,na.rm=T)/sqrt(ns*(m-1)/m))
      } else se <- arr2df(apply(aperm(data,c(Di$S,wsFs)),wsFs,sd,na.rm=T)/sqrt(ns))
      if (is.numeric(bars)) se$y <- -se$y*qt((100-bars)/200,as.vector(ns))
    }
    if (bar.type=="data") {
      stop("Sorry data resampling option not implemented")
      if ( !is.null(qp) & is.null(rawdat) ) 
        stop("Provide raw.data for bootstrap standard errors for quantile fits")
      if ( !is.null(qp) ) 
        se <- boot.se(rawdat$dat,bars,wsFs) else se <- boot.se(rc$dat,bars)
    }
    if (bar.type=="model") LMU=model.bar(rc,stat,measure,wsFs,probs,pb)
  } 
  av <- arr2df(apply(data,wsFs,mean,na.rm=T))  
  if (length(wsFs)==1) {
    names(av)[1] <- wsFs
    av[[wsFs]] <- factor(as.character(av[[wsFs]]),levels=c(levels(dat[[wsFs]])))
  }
  if ( plot.model ) {
    ltys <- 1:2
    mdata <- getstat(dat,measure,stat,probs,wsFs,Di,F,T)   
    av.m <- arr2df(apply(mdata,wsFs,mean,na.rm=T))
    if (length(wsFs)==1) names(av.m)[1] <- wsFs    
    av=cbind(Type=factor(rep(c(data.name,model.name),each=dim(av)[1]),
                         levels=c(data.name,model.name)),rbind(av,av.m))
  } else if (bar.no.line) {
    av=cbind(Type=factor(rep(c(data.name),each=dim(av)[1])),av)
    ltys=1
  } else {
    if (length(wsFs)==1) names(av)[1] <- wsFs
    ltys <- rep(1,length(levels(av[[1]])))
  }
  if ( !is.na(bars) ) { # add to lines factor
    av[,1] <- factor(as.character(av[,1]),c(levels(av[,1]),c("lo","hi")))
    if ( plot.model ) lo <- av[av$Type==data.name,] else lo <- av
    lo[,1] <- "lo"; hi <- lo; hi[,1] <- "hi"
    if (bar.type %in% c("within","between")) {
      lo$y <- lo$y-se$y
      hi$y <- hi$y+se$y
    } else {
      lo$y <- LMU[1,]
      hi$y <- LMU[3,]
      av[av[,1]==model.name,"y"] <- LMU[2,]
    }
    av <- rbind(av,lo,hi)
  }
  if (plot.model) pchs <- c(16,1) else
    pchs <- 1:length(levels(av[[1]])[!(levels(av[[1]]) %in% c("lo","hi"))])
  cols <- "black"
  if (plot.lines) if (bar.no.line) key=NULL else
    key=list(text=list(levels(av[,wsF[1]])[!(levels(av[,wsF[1]]) %in% c("lo","hi"))]),
        lines=list(lty=ltys),points=list(pch=pchs,cex=cex),columns=2)
  if (is.na(xlab))
    if (plot.lines) xlab <- wsF[2] else
      xlab <- wsF[1]
  
  if (do.plot) {  
    if ( any(is.na(xlim)) & any(is.na(ylim)) ) {
      if ( !plot.lines ) 
        print(xyplot(formula(plotform),av,xlab=xlab,col=cols, #panel=panel.se,
                     type=c("p","l"),ylab=ylab,layout=layout,as.table=T,
                     strip = strip.custom(bg="white"))) else 
        print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
                     type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
                     layout=layout,as.table=T,xlab=xlab,panel=panel.se,
                   strip = strip.custom(bg="white")))
    } else if (!any(is.na(xlim)) & any(is.na(ylim))) {
      if (!plot.lines) 
        print(xyplot(formula(plotform),av,xlab=xlab,col=cols, # panel=panel.se,
                     type=c("p","l"),ylab=ylab,xlim=xlim,layout=layout,as.table=T,
                     strip = strip.custom(bg="white"))) else 
        print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
          type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
          layout=layout,xlim=xlim,as.table=T,xlab=xlab,panel=panel.se,
                     strip = strip.custom(bg="white")))    
    } else if (any(is.na(xlim)) & !any(is.na(ylim))) {
      if (!plot.lines) 
        print(xyplot(formula(plotform),av,xlab=xlab,col=cols, # panel=panel.se,
          type=c("p","l"),ylab=ylab,ylim=ylim,as.table=T,
                     strip = strip.custom(bg="white"))) else 
        print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
          type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
          layout=layout,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
                     strip = strip.custom(bg="white")))    
    } else {
      if (!plot.lines) 
        print(xyplot(formula(plotform),av,xlab=xlab,col=cols, # panel=panel.se,
          type=c("p","l"),ylab=ylab,xlim=xlim,ylim=ylim,as.table=T,
                     strip = strip.custom(bg="white"))) else 
        print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
          type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
          layout=layout,xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
                     strip = strip.custom(bg="white")))    
    }
  }
  if (save.dat) av
}

# measure = c("rt","correct","error")[2];factors=NA;probs=c(.1,.5,.9);xlim=NA;ylim=NA                                    # plot limits
# xlab=NA;plot.subjects=F;bars="se";bar.type=c("within","between","data","model")[1]
# include.subjects=NA;exclude.subjects=NA;exclude.incomplete=F;data.fit.name="Type"                               # data/model factor name
# plot.model=F;model.name="fit";data.name="data";save.dat=FALSE;cex=1;ylab.qs=TRUE;
# layout=NULL
# 
# rc=ejldt1pc;plot.model=T; bar.type=c("within","between","data","model")[4]
#  

plotqs.rc <- function(rc,
  # all/correct/error rt
  measure = c("rt","correct","error")[1],  
  factors=NA,                                         # subset of factor names
  probs=c(.1,.5,.9),                                  # quantile probability
  xlim=NA,ylim=NA,                                    # plot limits
  xlab=NA,
  plot.subjects=F,                                    # subjects as panels
  bars="se",                                          # if NOT subjects do "se" or % CI (if numeric)
  bar.type=c("within","between","data","model","pbmodel")[2],
      # Morey 2008, sd/sqrt(n), resample data (not fully implemented) or 
      # simulate model or parameteric bootstap model
  include.subjects=NA,                                # subjects to include in analysis NA => all
  exclude.subjects=NA,
  exclude.incomplete=F,                               # dump entire subject where any cell NA
  data.fit.name="Type",                               # data/model factor name
  plot.model=F,                                       # data and model as lines
  model.name="fit",
  data.name="data",
  save.dat=FALSE,
  cex=1,
  ylab.qs=TRUE,
  plot.key=TRUE,
  layout=NULL) {                                      # xyplot key layout parameter
  # If no factors values specified sets to all factors
  # First factor is x, subsequent factors are panels
  # If plot.subjects add subjects as last factor (i.e., panels)
  # If plot.model adds data vs. fit as first factor (i.e., lines) 
  # layout length=2 as c(0,n) => panels/pages, c(x,y) x columns, y rows, c(x,y,z) z=num pages
  
  panel.se <- function(x, y, groups, subscripts, ...) {
    if (any(groups %in% c("lo","hi"))) { # doing ses
      do.se <- T; is.datafit <- !(groups=="lo" | groups=="hi")      
      tmp <- unlist(lapply(strsplit(as.character(groups),".",fixed=T),function(x){x[1]}))
      if (any(tmp==model.name)) {
        is.lohi <- c(1:length(x))>length(x)/2
        is.data <- c(1:length(x))<=length(x)/4
      } else {
        is.lohi <- c(1:length(x))>length(x)/3
        is.data <- c(1:length(x))<=length(x)/3        
      }
      y.lohi <- y[is.lohi]
      if (bar.type!="model") 
        y.0 <- c(y[is.data],y[is.data]) else
        y.0 <- c(y[!is.data & !is.lohi],y[!is.data & !is.lohi])
      x.0 <- c(x[is.data],x[is.data]) 
      x <- x[is.datafit]; y <- y[is.datafit]; groups <- groups[is.datafit]
    } else do.se <- F
    panel.xyplot(x,y,groups=groups,subscripts=subscripts, ...)
    #      panel.text(4,8,paste(round(y.lohi,1),collapse="|"))
    if (do.se) for (i in 1:length(y.0)) if (bar.type!="model")
      panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1) else
        panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1,lty=2)
  }  
  
  check <- function(wsF,allF,data.fit.name) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
    if (any(wsF==data.fit.name))
      stop(paste("data.fit.name (",data.fit.name,,
                 "and factor name are the same, change the former"))
  }  
  
  # verbose=F;sim=T;qnam=data.fit.name
  
  getstats <- function(dat,measure,probs,wsF,Di,verbose=F,sim=F,qnam="qs") {
    
    getstat <- function(dat,measure,prob,wsF,Di,verbose=F,sim=F) {
      if (sim) {
        rtnam=paste(Di$RT,"sim",sep=".") 
        qpnam="qp.sim"
        if (!any(names(dat)==rtnam))
          stop("Simulation results have not been added to rc object")
      } else {
        rtnam <- Di$RT
        qpnam="qp"
      }
      if (verbose) cat(paste(" at p =",probs[1],"\n"))
      if (!is.null(attr(dat,"qp"))) {
        tmp <- dat[dat$qp==prob & is.finite(dat[[Di$RT]]),]
        data <- tapply(tmp[,rtnam],tmp[,wsF],function(x){mean(x)})
      } else data <- tapply(dat[,rtnam],dat[,wsF],quantile,probs=prob)
      data
    }
    
    if (sim) qpnam="qp.sim" else qpnam="qp"
    dat <- dat[is.finite(dat[,Di$RT]),]    
    tmp <- vector(mode="list",length=length(probs))
    for (i in 1:length(probs))
      tmp[[i]] <- getstat(dat,measure,probs[i],dwsFs,Di,verbose=verbose,sim=sim) 
    tmp.dn <- dimnames(tmp[[1]])
    dn <- vector(mode="list",length=length(tmp.dn))
    for (i in 1:length(tmp.dn)) dn[[i]] <- tmp.dn[[i]]
    dn[[length(tmp.dn)+1]] <- as.character(probs)
    names(dn) <- c(names(dimnames(tmp[[1]])),qnam)
    d <- c(dim(tmp[[1]]),length(probs))
    data <- array(unlist(tmp),dim=d,dimnames=dn)
    aperm(data,c(length(tmp.dn)+1,1:length(tmp.dn)))
  }
  
  
  model.bar <- function(rc,measure,wsFs,probs=NA,pb=FALSE) {
    dname <- attr(rc$dat,"dname")
    model <- attr(rc$dat,"model")
    snams <- levels(rc$dat[[Di$S]])
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
        pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if (measure=="accuracy" | measure=="correct")  # target response to tabulate
      tr <- as.character(Di$D$R[as.logical(Di$D[[Di$SC[1]]])]) else
      tr <- as.character(Di$D$R[!as.logical(Di$D[[Di$SC[1]]])])
    cat("Calculating model errors for each subject ")
    for (s in 1:length(snams)) { # subjects
      if (!pb) load(paste(dname,"/",snams[s],"/sims/",model,"sims.RData",sep="")) else
        load(paste(dname,"/",snams[s],"/sims/",model,"=",model,".RData",sep=""))  
      if (s==1) av <- array(0,dim=c(length(sims),length(probs),length(sims[[1]]),length(snams)))
      cat(".")
      for (i in 1:length(tr)) {
        av[i,,,s] <- unlist(lapply(sims[[i]],function(x){
          do.call(quantile,list(x=x$rt[x$r==tr[[i]]],probs=probs,na.rm=TRUE))
        }))
      }
    }
    cat("\n")
    av <- apply(av,1:3,mean,na.rm=TRUE)
    Di <- attr(rc$dat,"D")[[1]]
    D <- Di$D[as.logical(Di$D[[Di$SC[1]]]),]
    # probs, pci, conditions
    aperm(apply(apply(av,2:3,function(x){
        as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))
      }),1:2,quantile,probs=pci,na.rm=TRUE),c(1,3,2))
  }
  
  
  # main body
  if (bar.type=="pbmodel") {
    bar.type="model"
    pb <- TRUE
  } else pb <- FALSE
  if (plot.subjects) bars <- NA
  if (length(probs)<2)
    stop("Use plot.rc to plot a single quantile")
  probs <- sort(probs)
  require(lattice)
  dat <- rc[[1]]
  qp <- attr(dat,"qp")
  D <- attr(dat,"D")
  Di <- D[[1]]
  
  if (!any(is.na(exclude.subjects))) {
    dat <- dat[!(dat[,Di$S] %in% exclude.subjects),]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
  }  
  
  # Get factors
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
#   check(wsF,Di$F,data.fit.name)
  if (plot.subjects) {
    wsF <- c(wsF,Di$S)
    bars <- NA
  }
  wsF <- c(data.fit.name,wsF)
  plot.lines <- T
  wsFs <- wsF[-1]
  if (plot.subjects) dwsFs <- wsFs else dwsFs <- c(wsFs,Di$S)
  # plotting formula
  plotform <- paste("y ~",wsFs[1])
  if (length(wsFs)>1) if (length(wsFs)==2) plotform <- paste(plotform,"|",wsFs[2]) else
    plotform <- paste(plotform,"|",paste(wsFs[-1],collapse="+"))
  # Get measure
  measures <- c("rt","correct","error")
  if (is.na(pmatch(measure[1],measures)))
    stop("Not a supported measure (rt, correct, error)")
  if (!is.null(qp)) {
    if (measure=="rt")
      stop("rt measure not available for quantile data")
    if (!all(probs %in% qp))
      stop(paste("Quantile probs must be in set",paste(qp,collapse=","),"\n"))
  }
  dat <- switch(pmatch(measure[1],measures),
    {ylab <- "All RT (s):"; dat},
    {ylab <- "Correct RT (s):"; dat[as.logical(dat[,Di$SC[1]]),]},
    {ylab <- "Error RT (s):"; dat[!as.logical(dat[,Di$SC[1]]),]},
  )  
  if (!ylab.qs) ylab <- paste(ylab,"quantiles") else
    ylab <- paste(ylab,paste(probs,collapse=" ,"),"quantiles")
  data <- getstats(dat,measure,probs,dwsFs,Di,qnam=data.fit.name)
  sdim <- c(1:length(dim(data)))[names(dimnames(data))==Di$S]
  snams <- dimnames(data)[[sdim]]
  bad <- apply(data,sdim,function(x){any(is.na(x))})
  if (all(bad))
    stop("All participants have missing data cells")
  if (exclude.incomplete && any(bad)) {
    cat(paste("Excluding participants with missing data cells",
      paste(snams[bad],collapse=","),"(",sum(!bad),"left )\n"))  
    dat <- dat[dat[,Di$S] %in% snams[!bad],]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstats(dat,measure,probs,dwsFs,Di,qnam=data.fit.name) 
  }
  if (!any(is.na(include.subjects))) {
    if (!all(include.subjects %in% snams))
      stop("Some specifed subjects not in set of subjects")
    cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
    dat <- dat[dat[,Di$S] %in% include.subjects,]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstats(dat,measure,probs,dwsFs,Di,qnam=data.fit.name) 
  }
  if (plot.model) {
    av <- arr2df(apply(data,wsF,mean,na.rm=T))
    mdata <- getstats(dat,measure,probs,wsFs,Di,F,T,qnam=data.fit.name)
    av.m <- arr2df(apply(mdata,wsF,mean,na.rm=T))
    av <- rbind(av,av.m)
    av[,1] <- paste(rep(c(data.name,model.name),each=dim(av.m)[1]),av[,1],sep=".")
    av[,1] <- factor(av[,1],unique(av[,1]))
    ltys <- rep(1:2,each=length(levels(av.m[,1])))
    pchs <- rep(c(16,1),each=length(levels(av.m[,1])))
  } else {
    av <- arr2df(apply(data,wsF,mean,na.rm=T))
    ltys <- rep(1,length(levels(av[[1]])))
    pchs <- 1
  }
  if (!is.na(bars)) {
    if (bar.type=="model" & !plot.model)
      stop("Set plot.model=TRUE to enable model error-bar plotting")
    if (!(any(bar.type %in% c("within","between","data","model")))) 
      stop("Parameter bar.type must be one of \"within\",\"between\", \"data\" or \"model\"")
    facs <- c(data.fit.name,wsFs)    
    if ( bar.type %in% c("within","between") ) {
      ns <- apply(data,facs,function(x){sum(!is.na(x))})
      if (bar.type=="within") {
        m <- prod(unlist(lapply(dimnames(data)[facs],length)))
        se <- arr2df(apply(aperm(data,c(Di$S,facs)) - apply(data,Di$S,mean,na.rm=T),
                           facs,sd,na.rm=T)/sqrt(ns*(m-1)/m))
      } else se <- arr2df(apply(aperm(data,c(Di$S,facs)),facs,sd,na.rm=T)/sqrt(ns))
      if (is.numeric(bars)) 
        se$y <- -se$y*qt((100-bars)/200,as.vector(ns))
    } 
    if ( bar.type=="data" ) {
      stop("Sorry data resampling option not implemented")
    }
    if (bar.type=="model") LMU=model.bar(rc,measure,wsFs,probs)

    av[,1] <- factor(av[,1],c(levels(factor(av[,1])),c("lo","hi")))
    if (!plot.model) lo=av else {
      is.dat <- unlist(lapply(strsplit(
        as.character(av$Type),".",fixed=T),function(x){x[1]}))==data.name
      lo <- av[is.dat,]
    }
    lo[,1] <- "lo"; hi <- lo; hi[,1] <- "hi"
    if ( bar.type %in% c("within","between") ) {
      lo$y <- lo$y-se$y
      hi$y <- hi$y+se$y
    } else {
      lo$y <- as.vector(LMU[1,,])
      hi$y <- as.vector(LMU[3,,])
      av[!is.dat,"y"] <- as.vector(LMU[2,,])
    }
    av <- rbind(av,lo,hi)
    
  }
  cols <- "black"  
  if (plot.model) key=list(text=list( c(data.name,model.name) ),
    lines=list(lty=1:2),points=list(pch=c(16,1),cex=cex),columns=2)
  if (!plot.key) key <- NULL
  if (is.na(xlab))
    if (plot.lines) xlab <- wsF[2] else
      xlab <- wsF[1]  
  if (any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
        type=c("p","l"),ylab=ylab,layout=layout,as.table=T,col=cols,panel=panel.se,
                            strip = strip.custom(bg="white"))) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white")))
  } else if (!any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,panel=panel.se,
        type=c("p","l"),ylab=ylab,xlim=xlim,layout=layout,as.table=T,col=cols,
                            strip = strip.custom(bg="white"))) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,xlim=xlim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white")))    
  } else if (any(is.na(xlim)) & !any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
        type=c("p","l"),ylab=ylab,ylim=ylim,as.table=T,col=cols,panel=panel.se,
                            strip = strip.custom(bg="white"))) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white")))    
  } else {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
        type=c("p","l"),ylab=ylab,xlim=xlim,ylim=ylim,as.table=T,col=cols,panel=panel.se,
                            strip = strip.custom(bg="white"))) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white")))    
  }
  if (save.dat) av
}

# factors=NA; do.mean=F; prob=.5; xlim=NA;xlab=NA; plot.subjects=F
# bars="se"; bar.type=c("within","between","data","model")[1]
# include.subjects=NA; exclude.subjects=NA; exclude.incomplete=F
# data.fit.name="Type"; plot.model=F; model.name="fit"; data.name="data"
# correct.data.name="Correct data"; error.data.name="Error data"
# cex=1;error.col="black";layout=NULL

plotce.rc <- function(rc,
  factors=NA,                 # subset of factor names
  do.mean=F,                  # for ML fits only will do mean
  prob=.5,                    # quantile probability
  xlim=NA,ylim=NA,            # plot limits
  xlab=NA,
  plot.subjects=F,            # subjects as panels
  bars="se",                  # if NOT subjects do "se" or % CI (if numeric)
  bar.type=c("within","between","data","model","pbmoel")[2],
    # Morey 2008, sd/sqrt(n), resample data (not fully implemented) or 
    # simulate model or parameteric bootstrap
  include.subjects=NA,        # subjects to include in analysis NA => all
  exclude.subjects=NA,
  exclude.incomplete=F,       # dump entire subject where any cell NA
  data.fit.name="Type",       # data/model factor name
  plot.model=F,               # data and model as lines
  model.name="fit",
  data.name="data",
  correct.data.name="Correct data",
  error.data.name="Error data",
  cex=1,
  error.col="black",
  plot.key=TRUE,
  layout=NULL) {              # xyplot key layout parameter
  # If no factors values specified sets to all factors
  # First factor is x, subsequent factors are panels
  # If plot.subjects add subjects as last factor (i.e., panels)
  # If plot.model adds data vs. fit as first factor (i.e., lines) 
  # layout length=2 as c(0,n) => panels/pages, c(x,y) x columns, y rows, c(x,y,z) z=num pages
  
  panel.se <- function(x, y, groups, subscripts, ...) {
    if (any(groups %in% c("lo","hi"))) { # doing ses
      do.se <- T; is.datafit <- !(groups=="lo" | groups=="hi")      
      tmp <- unlist(lapply(strsplit(as.character(groups),".",fixed=T),function(x){x[1]}))
      if (any(tmp==model.name)) {
        is.lohi <- c(1:length(x))>length(x)/2
        is.data <- c(1:length(x))<=length(x)/4
      } else {
        is.lohi <- c(1:length(x))>length(x)/3
        is.data <- c(1:length(x))<=length(x)/3        
      }
      y.lohi <- y[is.lohi]
      if (bar.type!="model") 
        y.0 <- c(y[is.data],y[is.data]) else
        y.0 <- c(y[!is.data & !is.lohi],y[!is.data & !is.lohi])
      x.0 <- c(x[is.data],x[is.data]) 
      x <- x[is.datafit]; y <- y[is.datafit]; groups <- groups[is.datafit]
    } else do.se <- F
    panel.xyplot(x,y,groups=groups,subscripts=subscripts, ...)
    #      panel.text(4,8,paste(round(y.lohi,1),collapse="|"))
    if (do.se) for (i in 1:length(y.0)) if (bar.type!="model")
      panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1) else
        panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1,lty=2)
  }  

  check <- function(wsF,allF,data.fit.name) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
    if (any(wsF==data.fit.name))
      stop(paste("data.fit.name (",data.fit.name,
                 "and factor name are the same, change the former"))
  }  
  
  # verbose=F;sim=T;qnam=data.fit.name
  
  getstats <- function(dat,prob,wsF,Di,verbose=F,sim=F,qnam="qs",do.mean=F) {
    
    getstat <- function(dat,measure,prob,wsF,Di,verbose=F,sim=F,do.mean=F) {
      if (sim) {
        rtnam=paste(Di$RT,"sim",sep=".") 
        qpnam="qp.sim"
        if (!any(names(dat)==rtnam))
          stop("Simulation results have not been added to rc object")
      } else {
        rtnam <- Di$RT
        qpnam="qp"
      }
      if (verbose) cat(paste(" at p =",prob,"\n"))
      if (!is.null(attr(dat,"qp"))) {
        tmp <- dat[dat$qp==prob & is.finite(dat[[Di$RT]]),]
        data <- tapply(tmp[,rtnam],tmp[,wsF],function(x){mean(x)})
      } else if (do.mean) data <- tapply(dat[,rtnam],dat[,wsF],mean) else
        data <- tapply(dat[,rtnam],dat[,wsF],quantile,probs=prob)
      data
    }
    
    if (sim) qpnam="qp.sim" else qpnam="qp"
    dat <- dat[is.finite(dat[,Di$RT]),]    
    tmp <- vector(mode="list",length=2)
    tmp[[1]] <- getstat(dat[as.logical(dat[,Di$SC[1]]),],
      "correct",prob,dwsFs,Di,verbose=verbose,sim=sim,do.mean=do.mean) 
    tmp[[2]] <- getstat(dat[!as.logical(dat[,Di$SC[1]]),],
      "error",prob,dwsFs,Di,verbose=verbose,sim=sim,do.mean=do.mean) 
    tmp.dn <- dimnames(tmp[[1]])
    dn <- vector(mode="list",length=length(tmp.dn))
    for (i in 1:length(tmp.dn)) dn[[i]] <- tmp.dn[[i]]
    dn[[length(tmp.dn)+1]] <- c("Correct","Error")
    names(dn) <- c(names(dimnames(tmp[[1]])),qnam)
    d <- c(dim(tmp[[1]]),2)
    data <- array(unlist(tmp),dim=d,dimnames=dn)
    aperm(data,c(length(tmp.dn)+1,1:length(tmp.dn)))
  }
  
  model.bar <- function(rc,wsFs,do.mean=F,prob=NA,pb=FALSE) {
    dname <- attr(rc$dat,"dname")
    model <- attr(rc$dat,"model")
    snams <- levels(rc$dat[[Di$S]])
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
        pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
        pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if (do.mean) stat <- "mean" else stat <- "quantile"
    cat("Calculating model errors for each subject ")
    for (s in 1:length(snams)) { # subjects
      if (!pb) load(paste(dname,"/",snams[s],"/sims/",model,"sims.RData",sep="")) else
        load(paste(dname,"/",snams[s],"/sims/",model,"=",model,".RData",sep=""))  
      cat(".")
      if (s==1) {
        avc <- array(dim=c(length(sims),length(sims[[1]]),length(snams)))
        ave <- avc
      }
      tr <- as.character(Di$D$R[as.logical(Di$D[[Di$SC[1]]])]) # correct
      for (i in 1:length(tr)) avc[i,,s] <- unlist(lapply(sims[[i]],function(x){
          do.call(stat,list(x=x$rt[x$r==tr[[i]]],probs=prob,na.rm=TRUE))
      }))
      tr <- as.character(Di$D$R[!as.logical(Di$D[[Di$SC[1]]])]) # error
      for (i in 1:length(tr)) ave[i,,s] <- unlist(lapply(sims[[i]],function(x){
          do.call(stat,list(x=x$rt[x$r==tr[[i]]],probs=prob,na.rm=TRUE))
      }))
    }
    cat("\n")
    avc <- apply(avc,1:2,mean,na.rm=TRUE)
    ave <- apply(ave,1:2,mean,na.rm=TRUE)
    Di <- attr(rc$dat,"D")[[1]]
    D <- Di$D[as.logical(Di$D[[Di$SC[1]]]),]
    list(c=apply(apply(avc,2,function(x){
      as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),1,quantile,probs=pci,na.rm=TRUE),
      e=apply(apply(ave,2,function(x){
        as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),1,quantile,probs=pci,na.rm=TRUE)
    )
  }
  
  require(lattice)
  if (bar.type=="pbmodel") {
    bar.type="model"
    pb <- TRUE
  } else pb <- FALSE
  if (plot.subjects) bars <- NA
  cex <- cex*c(1.5,1.5,1,1)
  dat <- rc[[1]]
  qp <- attr(dat,"qp")
  if (do.mean & !is.null(qp))
    stop("Cant do mean on quantile data")
  D <- attr(dat,"D")
  Di <- D[[1]]
  
  if (!any(is.na(exclude.subjects))) {
    dat <- dat[!(dat[,Di$S] %in% exclude.subjects),]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
  }  
  
  # Get factors
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
#   check(wsF,Di$F,data.fit.name)
  if (plot.subjects) {
    wsF <- c(wsF,Di$S)
    bars=NA
  }
  wsF <- c(data.fit.name,wsF)
  plot.lines <- T
  wsFs <- wsF[-1]
  if (plot.subjects) dwsFs <- wsFs else dwsFs <- c(wsFs,Di$S)
  # plotting formula
  plotform <- paste("y ~",wsFs[1])
  if (length(wsFs)>1) if (length(wsFs)==2) plotform <- paste(plotform,"|",wsFs[2]) else
    plotform <- paste(plotform,"|",paste(wsFs[-1],collapse="+"))
  if (!is.null(qp)) {
    if (!all(prob %in% qp))
      stop(paste("Quantile prob must be in set",paste(qp,collapse=","),"\n"))
  }
  if (do.mean) ylab <- paste("RT (s): Mean") else
    ylab <- paste("RT (s):",prob,"quantile")
  data <- getstats(dat,prob,dwsFs,Di,qnam=data.fit.name,do.mean=do.mean)
  sdim <- c(1:length(dim(data)))[names(dimnames(data))==Di$S]
  snams <- dimnames(data)[[sdim]]
  bad <- apply(data,sdim,function(x){any(is.na(x))})
  if (exclude.incomplete & all(bad))
    stop("All participants have missing data cells")
  if (exclude.incomplete && any(bad)) {
    cat(paste("Excluding participants with missing data cells",
              paste(snams[bad],collapse=","),"(",sum(!bad),"left )\n"))  
    dat <- dat[dat[,Di$S] %in% snams[!bad],]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstats(dat,prob,dwsFs,Di,qnam=data.fit.name,do.mean=do.mean) 
  }
  if (!any(is.na(include.subjects))) {
    if (!all(include.subjects %in% snams))
      stop("Some specifed subjects not in set of subjects")
    cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
    dat <- dat[dat[,Di$S] %in% include.subjects,]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
    data <- getstats(dat,prob,dwsFs,Di,qnam=data.fit.name,do.mean=do.mean) 
  }
  if (plot.model) {
    av <- arr2df(apply(data,wsF,mean,na.rm=T))
    mdata <- getstats(dat,prob,wsFs,Di,F,T,qnam=data.fit.name,do.mean=do.mean)
    av.m <- arr2df(apply(mdata,wsF,mean,na.rm=T))
    av <- rbind(av,av.m)
    av[,1] <- factor(paste(rep(c(data.name,model.name),each=dim(av.m)[1]),av[,1],sep="."),
      levels=c(unique(paste(rep(c(data.name,model.name),each=dim(av.m)[1]),av[,1],sep="."))))
    ltys <- rep(1:2,each=length(levels(av.m[,1])))
    pchs=c(16,17,1,2)
   } else {
    av <- arr2df(apply(data,wsF,mean,na.rm=T))
    ltys <- rep(1,length(levels(av[[1]])))
    pchs <- c(16,17,1,2)
  }
  if (!is.na(bars)) {
    if (bar.type=="model" & !plot.model)
      stop("Set plot.model=TRUE to enable model error-bar plotting")
    if (!(any(bar.type %in% c("within","between","data","model")))) 
      stop("Parameter bar.type must be one of \"within\",\"between\", \"data\" or \"model\"")
    facs <- c(data.fit.name,wsFs)
    if (bar.type %in% c("within","between") ) {
      ns <- apply(data,facs,function(x){sum(!is.na(x))})
      if (bar.type=="within") {
        m <- prod(unlist(lapply(dimnames(data)[facs],length)))
        se <- arr2df(apply(aperm(data,c(Di$S,facs)) - apply(data,Di$S,mean,na.rm=T),
                           facs,sd,na.rm=T)/sqrt(ns*(m-1)/m))
      } else se <- arr2df(apply(aperm(data,c(Di$S,facs)),facs,sd,na.rm=T)/sqrt(ns))
      if (is.numeric(bars)) 
        se$y <- -se$y*qt((100-bars)/200,as.vector(ns))
    }     
    if (bar.type=="data") 
      stop("Sorry data resampling option not implemented")
    if (bar.type=="model") LMU=model.bar(rc,wsFs,do.mean,prob)    
    av[,1] <- factor(av[,1],c(levels(factor(av[,1])),c("lo","hi")))
    if (!plot.model) lo=av else {
      is.dat <- unlist(lapply(strsplit(
        as.character(av$Type),".",fixed=T),function(x){x[1]}))==data.name
      lo <- av[is.dat,]
    }
    lo[,1] <- "lo"; hi <- lo; hi[,1] <- "hi"

    if (bar.type %in% c("within","between")) {
      lo$y <- lo$y-se$y
      hi$y <- hi$y+se$y
    } else {
      odd <- even <- c(1:dim(lo)[1])
      odd <- odd[odd%%2==1]; even <- even[even%%2==0]
      lo$y[odd] <- LMU$c[1,]
      lo$y[even] <- LMU$e[1,]
      hi$y[odd] <- LMU$c[3,]
      hi$y[even] <- LMU$e[3,]
      av$y[!is.dat][odd] <- LMU$c[2,]
      av$y[!is.dat][even] <- LMU$e[2,]
    }
    av <- rbind(av,lo,hi) 
  }
  cols <- c("black",error.col) # correct error 
  keycols=c("black","black",error.col,error.col)
  if (plot.model) 
    key=list(text=list(
      c(correct.data.name,model.name,error.data.name,model.name),col=keycols),
      lines=list(lty=c(1:2,1:2),col=keycols),
      points=list(pch=c(16,1,17,2),cex=cex[c(1,3,2,4)],col=keycols),
      columns=4) else
    key=list(text=list(
      c("Correct","Error"),col=cols),
      lines=list(lty=c(1:1),col=cols),
      points=list(pch=c(16,17),cex=cex[c(1,2)],col=cols),
      columns=2)
  if (!plot.key) key <- NULL
  if (is.na(xlab))
    if (plot.lines) xlab <- wsF[2] else
      xlab <- wsF[1]  
  if (any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,key=key,panel=panel.se,
        type=c("p","l"),ylab=ylab,layout=layout,as.table=T,col=cols,pch=pchs,cex=cex,
             strip = strip.custom(bg="white")) else 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white"))
  } else if (!any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,key=key,panel=panel.se,cex=cex,
        type=c("p","l"),ylab=ylab,xlim=xlim,layout=layout,as.table=T,col=cols,pch=pchs,
             strip = strip.custom(bg="white")) else 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
      layout=layout,xlim=xlim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white"))    
  } else if (any(is.na(xlim)) & !any(is.na(ylim))) {
    if (!plot.model) 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,key=key,pch=pchs,cex=cex,
        type=c("p","l"),ylab=ylab,ylim=ylim,as.table=T,col=cols,panel=panel.se,
             strip = strip.custom(bg="white")) else 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white"))    
  } else {
    if (!plot.model) 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,key=key,panel=panel.se,cex=cex,
        type=c("p","l"),ylab=ylab,xlim=xlim,ylim=ylim,as.table=T,col=cols,pch=pchs,
             strip = strip.custom(bg="white")) else 
      xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
             strip = strip.custom(bg="white"))    
  }
}

# rc=ta;model=bestmod;parname="st";plot.subjects=T
# save.par=NA;do.plot=T;bars="SE"
# factors=NA; layout=NULL
# use.par=NA; xlab=NA; av.stat = "mean"
# trans=T;reduce=T;mhname=NA;ylim=NA; ylab=NA; save.par=F
# plot.lines=T;cex=1; exclude.subjects=NA
# xlim=NA; ylim=NA
# include.subjects=NA
# rename.levels=list()

plotpar.rc <- function(rc,
  model=NA,parname=1,trans=T,reduce=T,mhname=NA,
  factors=NA,             # subset of factor names
  xlim=NA,ylim=NA,        # plot limits
  plot.subjects=T,        # subjects as panels
  av.stat = "mean",
  include.subjects=NA,    # subjects to include in analysis NA => all
  exclude.subjects=NA,
  plot.lines=T,           # first factor as lines 
  layout=NULL,            # xyplot key layout parameter
  save.par=F,             # save plot data frame (av)
  use.par=NA,             # use as av instead of calculate 
  bars=NA,              # if NOT subjects do "se" or % CI (if numeric)
  xlab=NA,ylab=NA,
  do.plot=T,
  rename.levels=list(),
  reorder.levels=list(),
  cex=1) {
  # if no factors values specified sets to all factors
  # If plot.lines, plots first entry in factors as lines except if only one, 
  # then plots as x, otherwise second factor is x, further factors panels.
  # If not plot.lines first factor is x, subsequent factors are panels
  # if plot.subjects add subjects as last factor (i.e., panels)
  # if plot.model adds data vs. fit as first factor (i.e., lines) OVERRIDES plot.lines
  # layout length=2 as c(0,n) => panels/pages, c(x,y) x columns, y rows, c(x,y,z) 
  # z=num pages
  
  panel.se <- function(x, y, groups=NULL, subscripts, ...) {
    if (any(groups %in% c("lo","hi"))) { # doing ses
      do.se <- T; is.datafit <- !(groups=="lo" | groups=="hi")
      is.lohi <- c(1:length(x))>length(x)/3
      is.data <- c(1:length(x))<=length(x)/3        
      y.lohi <- y[is.lohi]; y.0 <- c(y[is.data],y[is.data])
      x.0 <- c(x[is.data],x[is.data]) 
      x <- x[is.datafit]; y <- y[is.datafit]; groups <- groups[is.datafit]
    } else do.se <- F
    panel.xyplot(x,y,groups=groups,subscripts=subscripts, ...)
    #      panel.text(4,8,paste(round(y.lohi,1),collapse="|"))
    if (do.se) for (i in 1:length(y.0)) {
      panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=.1)
    }    
  }  
  
  require(lattice)
  if (plot.subjects) bars <- NA
  av <- getpar.rc(rc,model=model,mhname=mhname,reduce=reduce,trans=trans)
  if (is.numeric(parname)) parname <- names(av)[parname]
  if (!(parname %in% names(av))) 
    stop("Parameter name not present in fit list")
  av <- av[[parname]] 
  if (length(rename.levels)) for (i in names(rename.levels)) {
    levels(av[[i]]) <- rename.levels[[i]]
    if (!is.null(reorder.levels[[i]]))
      av[[i]] <- factor(as.character(av[[i]]),levels=reorder.levels[[i]])
  }
  if ( dim(av)[2]==2 ) {
    names(av)[2] <- parname
    snam <- names(av)[1]
    if (!any(is.na(exclude.subjects))) {
      av <- av[!(av[[snam]] %in% exclude.subjects),]
      av[,snam] <- factor(as.character(av[,snam]))
    }
    if (!any(is.na(include.subjects))) {
      if (!all(include.subjects %in% snam))
        stop("Some specifed subjects not in set of subjects")
      cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
      av <- av[av[[snam]] %in% include.subjects,]
      av[,snam] <- factor(as.character(av[,snam]))
    }    
    if (do.plot) {
      if (plot.subjects) print(av) else 
        print(do.call(av.stat,list(x=av[,parname])))
    }
  } else {
    wsF <- names(av)
    snam <- wsF[1]    
    snams <- levels(factor(av[[snam]]))
    if ( !any(is.na(exclude.subjects)) ) {
      av <- av[!(av[[snam]] %in% exclude.subjects),]
      av[,snam] <- factor(as.character(av[,snam]))
    }
    if ( !any(is.na(include.subjects)) ) {
      if ( !all(include.subjects %in% snams) )
        stop("Some specifed subjects not in set of subjects")
      cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
      av <- av[av[[snam]] %in% include.subjects,]
      av[,snam] <- factor(as.character(av[,snam]))
    }
    av[[snam]] <- factor(av[[snam]])
    snams <- levels(av[[snam]])
    wsF <- wsF[-c(1,length(wsF))]
    if ( !all(is.na(factors)) ) {
      if ( !all(factors %in% wsF) )
        stop("Factor names not present in fit list")
      wsF <- factors
    }
    if (plot.subjects) {
      if (length(wsF)==1) plot.lines=F
      wsF <- c(wsF,snam)
    }
    if ( !is.na(bars) ) {
      if ( length(wsF)==1 ) m <- length(levels(av[,wsF])) else
        m <- prod(unlist(lapply(lapply(av[,wsF],levels),length)))
      ns <- tapply(av$y,av[,wsF],length)*(m-1)/m
      se <- av$s
      levels(se) <- tapply(av$y,av[,snam],mean)
      se <- av$y-as.numeric(as.character(se))
      se <- tapply(se,av[,wsF],sd)/sqrt(ns)
      if ( length(dim(se))==1 ) names(dimnames(se)) <- wsF
      se <- arr2df(se)
      if (is.numeric(bars)) 
        se$y <- -se$y*qt((100-bars)/200,as.vector(ns))
    }
    av <- tapply(av$y,av[,wsF],mean)
    names(dimnames(av)) <- wsF
    av <- arr2df(av)
    # plotting formula
    if ( length(wsF)==1 ) if ( is.na(bars) ) plot.lines <- F else { 
        wsF <- c("Type",wsF)
        av <- cbind.data.frame(Type=rep("data",dim(av)[1]),av)
        nokey <- T
    } else nokey <- F
    if ( plot.lines ) {
      plotform <- paste("y ~",wsF[2]) 
      if ( length(wsF)>2 ) {
        if (length(wsF)==3) plotform <- paste(plotform,"|",wsF[3]) else
          plotform <- paste(plotform,"|",paste(wsF[-c(1,2)],collapse="+"))
      }
    } else {
      plotform <- paste("y ~",wsF[1])
      if ( length(wsF)==2 ) plotform <- paste(plotform,"|",wsF[2]) else
        if ( length(wsF)>1 ) plotform <- paste(plotform,"|",paste(wsF[-1],collapse="+"))
    }
    if ( is.na(xlab) )
      if ( plot.lines ) xlab <- wsF[2] else
        xlab <- wsF[1]
    if ( is.na(ylab) ) ylab <- parname
    keylevs <- levels(av[,wsF[1]])
    if ( !all(is.na(use.par)) ) {
      if (!all(names(av)==names(use.par)))
        stop("Input parameters do not match plot specification")
      av <- use.par
      keylevs <- levels(av[,wsF[1]])
      if (!is.na(bars)) keylevs <- keylevs[!(keylevs %in% c("lo","hi"))]
    } else if ( !is.na(bars) ) { # add to lines factor
      av[,1] <- factor(as.character(av[,1]),c(levels(av[,1]),c("lo","hi")))
      lo <- av
      lo[,1] <- "lo"; hi <- lo; hi[,1] <- "hi"
      lo$y <- lo$y-se$y; hi$y <- hi$y+se$y
      av <- rbind(av,lo,hi)      
    }
    if ( length(dim(av))==1 ) {
      cat(paste(parname,"\n"))
      print(av) 
    } else {      
      ltys <- rep(1,length( keylevs ))
      if ( !plot.lines ) pchs <- 1 else
        pchs <- 1:length( keylevs )
      cols <- "black"
      if ( plot.lines ) if (nokey) key=NULL else key=list(text=list(keylevs),
        lines=list(lty=ltys),points=list(pch=pchs,cex=cex),columns=2)
      if ( do.plot ) {
        if ( any(is.na(xlim)) & any(is.na(ylim)) ) {
          if ( !plot.lines ) print(xyplot(formula(plotform),av,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            layout=layout,as.table=T,xlab=xlab,panel=panel.se,strip = strip.custom(bg="white"))) else 
           print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            layout=layout,as.table=T,xlab=xlab,panel=panel.se,strip = strip.custom(bg="white")))
        } else if (!any(is.na(xlim)) & any(is.na(ylim))) {
          if (!plot.lines) print(xyplot(formula(plotform),av,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            xlim=xlim,layout=layout,as.table=T,xlab=xlab,panel=panel.se,
            strip = strip.custom(bg="white"))) else 
          print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            layout=layout,xlim=xlim,as.table=T,xlab=xlab,panel=panel.se,
            strip = strip.custom(bg="white")))    
        } else if (any(is.na(xlim)) & !any(is.na(ylim))) {
          if (!plot.lines) print(xyplot(formula(plotform),av,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            ylim=ylim,as.table=T,xlab=xlab,layout=layout,panel=panel.se,
            strip = strip.custom(bg="white"))) else 
          print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            layout=layout,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
                       strip = strip.custom(bg="white")))    
        } else {
          if (!plot.lines) print(xyplot(formula(plotform),av,panel=panel.se,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,layout=layout,
                                        strip = strip.custom(bg="white"))) else 
          print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
            type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
            layout=layout,xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
                       strip = strip.custom(bg="white")))    
        }
      }
    }
  }
  if (save.par) av
}

# dname="di.LBA"
# mname="B~E*D & A~1 & v~D*C & sv~1 & ter~E*D & pc~1"
# ss=NA; name.plot=T
plot.oe <- function(dname,mname=1,ss=NA,name.plot=T,mark.outlier=T) {

  fitlist <- load(paste(dname,paste(strsplit(
    dname,".",fixed=T)[[1]][1],"RData",sep="."),sep="/"))
  fitlist <- get(fitlist)
  mnames <- names(fitlist[[2]]$MT[[1]])
  if (is.numeric(mname)) mname <- mnames[mname[1]]
  mname <- mname[1]
  if (!(mname %in% mnames))
    stop("Model not in hierarchy")
  snams <- names(fitlist[[2]]$M)
  if (any(is.na(ss))) ss <- snams else 
    if (is.numeric(ss)) ss <- snams[ss]
  if (!all(ss %in% snams))
    stop("Some subject names not in data file")
  dat <- get.simdat(dname,mname,ss)
  RTnam <- attr(dat,"D")[[1]]$RT
  RTsim <- paste(RTnam,"sim",sep=".")
  dat <- dat[is.finite(dat[[RTnam]]),]
  if (name.plot) tit <- paste(ss,collapse=" ") else tit=""
  outlier <- dat$lm/dat$lc < 1
  if (mark.outlier) tit = paste(tit,"(",sum(outlier),"outliers )") 
  plot(dat[[RTnam]],dat[[RTsim]],xlab="Observed",ylab="Predicted",main=tit)
  if (mark.outlier) 
    points(dat[[RTnam]][outlier],dat[[RTsim]][outlier],col="red")
  abline(0,1)
  fits <- fitlist[[2]]$fits
  dev <- sum(fits$stats[as.character(fits$stats.names[,"s"]) %in% ss &
    fits$stats.names[,"hmname"] %in% mname,"dev"])   
  text(max(dat[[RTnam]])*.9,min(dat[[RTsim]])*1.1,
    paste("D =",round(dev,2)))
}


plot.oes <- function(dname,mname=1,ss=NA) {
  fitlist <- load(paste(dname,paste(strsplit(
    dname,".",fixed=T)[[1]][1],"RData",sep="."),sep="/"))
  fitlist <- get(fitlist)
  mnames <- names(fitlist[[2]]$MT[[1]])
  if (is.numeric(mname)) mname <- mnames[mname[1]]
  mname <- mname[1]
  if (!(mname %in% mnames))
    stop("Model not in hierarchy")
  snams <- names(fitlist[[2]]$M)
  if (any(is.na(ss))) ss <- snams else 
    if (is.numeric(ss)) ss <- snams[ss]
  if (!all(ss %in% snams))
    stop("Some subject names not in data file")
  fits <- fitlist[[2]]$fits
  is.model <- fits$stats.names[,"hmname"] %in% mname
  dev <- fits$stats[is.model,"dev"]
  names(dev) <- fits$stats.names[is.model,"s"]
  snams <- names(sort(dev))
  cat("Deviance\n")
  print(round(sort(dev),1))
  for (i in snams) plot.oe(dname,mname,i,T)
}


# factors=NA;probs=c(.1,.5,.9);xlim=NA;xlim=NA;xlab="Accuracy (%)";plot.subjects=F;
# bars="se";bar.type=c("within","between","data","model","pbmodel")[2];len=.05                                            # bar cap length
# include.subjects=NA;exclude.subjects=NA;exclude.incomplete=F;data.fit.name="Type";plot.model=F
# model.name="fit"; data.name="data"; save.dat=FALSE; cex=1; ylab.qs=TRUE; layout=NULL
# 
# rc=lba1

plotqp.rc <- function(rc,
  factors=NA,                                         # subset of factor names
  probs=c(.1,.5,.9),                                  # quantile probability
  xlim=NA,ylim=NA,                                    # plot limits
  xlab="Accuracy (%)",
  plot.subjects=F,                                    # subjects as panels
  bars="se",                                          # if NOT subjects do "se" or % CI (if numeric)
  bar.type=c("within","between","data","model","pbmodel")[2],
      # Morey 2008, sd/sqrt(n), resample data (not fully implemented) or 
      # simulate model or parametric bootstrap
  len=.05,                                            # bar cap length
  include.subjects=NA,                                # subjects to include in analysis NA => all
  exclude.subjects=NA,
  exclude.incomplete=F,                               # dump entire subject where any cell NA
  data.fit.name="Type",                               # data/model factor name
  plot.model=F,                                       # data and model as lines
  model.name="fit",
  data.name="data",
  save.dat=FALSE,
  cex=1,
  ylab.qs=TRUE,
  layout=NULL) {                                      # xyplot key layout parameter
  # If no factors values specified sets to all factors
  # Probability of error/correct response for first factor is x, lines are quantiles,
  # remaining factors are panels
  # If plot.subjects add subjects as last factor (i.e., panels)
  # If plot.model adds data vs. fit as first factor (i.e., to lines) 
  # layout length=2 as c(0,n) => panels/pages, c(x,y) x columns, y rows, c(x,y,z) z=num pages
  
  panel.se <- function(x, y, groups, subscripts, ...) {
    if (any(groups %in% c("lo","hi"))) { # doing ses
      do.se <- T; is.datafit <- !(groups=="lo" | groups=="hi")      
      tmp <- unlist(lapply(strsplit(as.character(groups),".",fixed=T),function(x){x[1]}))
      if (any(tmp==model.name)) {
        is.lohi <- c(1:length(x))>length(x)/2
        is.data <- c(1:length(x))<=length(x)/4
      } else {
        is.lohi <- c(1:length(x))>length(x)/3
        is.data <- c(1:length(x))<=length(x)/3        
      }
      x.lohi <- x[is.lohi]
      y.lohi <- y[is.lohi]
      if (bar.type!="model") {
        y.0 <- c(y[is.data],y[is.data]) 
        x.0 <- c(x[is.data],x[is.data]) 
      } else {
        y.0 <- c(y[!is.data & !is.lohi],y[!is.data & !is.lohi])
        x.0 <- c(x[!is.data & !is.lohi],x[!is.data & !is.lohi]) 
      }  
      x <- x[is.datafit]; y <- y[is.datafit]; groups <- groups[is.datafit]
    } else do.se <- F
    panel.xyplot(x,y,groups=groups,subscripts=subscripts, ...)
    #      panel.text(4,8,paste(round(y.lohi,1),collapse="|"))
    if (do.se) {
      for (i in 1:length(y.0)) { 
        if (bar.type!="model") {
          panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=len) 
          panel.arrows(x.0[i],y.0[i],x.lohi[i],y.0[i],angle=90,length=len) 
        } else {
          panel.arrows(x.0[i],y.0[i],x.0[i],y.lohi[i],angle=90,length=len,lty=2)
          panel.arrows(x.0[i],y.0[i],x.lohi[i],y.0[i],angle=90,length=len,lty=2)
        }
      }
    }
  }  
  
  check <- function(wsF,allF,data.fit.name) {
    if (!all(wsF %in% allF))
      stop("Factors not in design")
    if (any(wsF==data.fit.name))
      stop(paste("data.fit.name (",data.fit.name,,
                 "and factor name are the same, change the former"))
  }  
  
  # verbose=F;sim=T;qnam=data.fit.name
  
  # verbose=F;sim=T; wsF=wsFs
  getstat <- function(dat,measure,stat,probs,wsF,Di,verbose=F,sim=F) {
    pc <-c("accuracy","inaccuracy","stop")
    if (!(measure %in% pc)) # && stat != "quantile") 
      dat <- dat[is.finite(dat[,Di$RT]),]
    if (sim) {
      rtnam=paste(Di$RT,"sim",sep=".") 
      qpnam="qp.sim"
      if (!any(names(dat)==rtnam))
        stop("Simulation results have not been added to rc object")
    } else {
      rtnam <- Di$RT
      qpnam="qp"
    }
    if (stat[1]=="quantile") {
      if (verbose) cat(paste(" at p =",probs[1],"\n"))
      if (!is.null(qp)) {
        tmp <- dat[dat$qp==probs[1] & is.finite(dat[[Di$RT]]),]
        if (any(is.na(tmp[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(tmp[,rtnam],tmp[,wsF],function(x){mean(x)})
      } else {
        if (any(is.na(dat[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(dat[,rtnam],dat[,wsF],stat,probs=probs[1])
      }
    } else {
      if (verbose) cat("\n")
      if (measure[1] %in% pc) {
        if (any(is.na(dat[,qpnam])))
          stop(paste("NAs present in statistic",qpnam))
        data <- tapply(dat[,qpnam],dat[,wsF],function(x){mean(x)}) 
        #        data[is.na(data)] <- 0    
      } else {
        if (any(is.na(dat[,rtnam])))
          stop(paste("NAs present in statistic",rtnam))
        data <- tapply(dat[,rtnam],dat[,wsF],stat)
      }
    }
    names(dimnames(data)) <- wsF
    if (measure[1] %in% pc) {
      data[is.na(data)] <- 0
      data*100
    } else data
  }
  
  
  fixholes <- function(dat,Di,measure="inaccuracy",plot.model=F) {
    
    cells <- function(d)
      apply(matrix(unlist(lapply(d,as.character)),ncol=dim(d)[2]),1,paste,collapse=" ")
    
    s <- Di$S
    #     f <- Di$F
    dc=dat[is.infinite(dat[,Di$RT]) & as.logical(dat[,Di$SC[1]]),]
    de=dat[is.infinite(dat[,Di$RT]) & !as.logical(dat[,Di$SC[1]]),]
    cc <- cells(dc[,c(s,"rcell")]); ce <- cells(de[,c(s,"rcell")])
    # assumes only error responses can be non-unique!
    dups <- duplicated(ce) 
    deu <- de[!dups,]   
    de <- de[dups,]
    row.names(dc) <- cc; row.names(deu) <- ce[!dups]
    ce <- ce[dups]
    if (length(ce)>0) for (i in 1:length(ce)) {
      deu[ce[i],"qp"] <- deu[ce[i],"qp"] + de[i,"qp"]
      if (plot.model) 
        deu[ce[i],"qp.sim"] <- deu[ce[i],"qp.sim"] + de[i,"qp.sim"]      
    }
    out <- deu
    missing <- cc[!(cc %in% row.names(out))]
    if ( length(missing)!=0 ) {
      out <- rbind(out,dc[missing,])
      out[missing,"qp"] <- 1-dc[missing,"qp"]
      if (plot.model) 
        out[missing,"qp.sim"] <- 1-dc[missing,"qp.sim"]
    }
    if (measure=="accuracy") {
      out[,"qp"] <- 1 - out[,"qp"]
      if (plot.model) 
        out[,"qp.sim"] <- 1 - out[,"qp.sim"]
      out[,Di$SC[1]] <- T
    } else out[,Di$SC[1]] <- F
    out
  }
   
  getstats <- function(dat,measure,probs,wsF,Di,verbose=F,sim=F,qnam="qs") {
    
    getstat <- function(dat,measure,prob,wsF,Di,verbose=F,sim=F) {
      if (sim) {
        rtnam=paste(Di$RT,"sim",sep=".") 
        qpnam="qp.sim"
        if (!any(names(dat)==rtnam))
          stop("Simulation results have not been added to rc object")
      } else {
        rtnam <- Di$RT
        qpnam="qp"
      }
      if (verbose) cat(paste(" at p =",probs[1],"\n"))
      if (!is.null(attr(dat,"qp"))) {
        tmp <- dat[dat$qp==prob & is.finite(dat[[Di$RT]]),]
        data <- tapply(tmp[,rtnam],tmp[,wsF],function(x){mean(x)})
      } else data <- tapply(dat[,rtnam],dat[,wsF],quantile,probs=prob)
      data
    }
    
    if (sim) qpnam="qp.sim" else qpnam="qp"
    dat <- dat[is.finite(dat[,Di$RT]),]    
    tmp <- vector(mode="list",length=length(probs))
    for (i in 1:length(probs))
      tmp[[i]] <- getstat(dat,measure,probs[i],dwsFs,Di,verbose=verbose,sim=sim) 
    tmp.dn <- dimnames(tmp[[1]])
    dn <- vector(mode="list",length=length(tmp.dn))
    for (i in 1:length(tmp.dn)) dn[[i]] <- tmp.dn[[i]]
    dn[[length(tmp.dn)+1]] <- as.character(probs)
    names(dn) <- c(names(dimnames(tmp[[1]])),qnam)
    d <- c(dim(tmp[[1]]),length(probs))
    data <- array(unlist(tmp),dim=d,dimnames=dn)
    aperm(data,c(length(tmp.dn)+1,1:length(tmp.dn)))
  }
  
  model.bar <- function(rc,wsFs,probs=NA,pb=FALSE) {
    dname <- attr(rc$dat,"dname")
    model <- attr(rc$dat,"model")
    snams <- levels(rc$dat[[Di$S]])
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
        pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if (!is.numeric(bars)) 
      pci <- c(pnorm(-1),.5,pnorm(1)) else 
        pci <- c((1-bars/100)/2,.5,1-(1-bars/100)/2)
    if (do.mean) stat <- "mean" else stat <- "quantile"
    cat("Calculating model errors for each subject ")
    for (s in 1:length(snams)) { # subjects
      if (!pb) load(paste(dname,"/",snams[s],"/sims/",model,"sims.RData",sep="")) else
        load(paste(dname,"/",snams[s],"/sims/",model,"=",model,".RData",sep=""))  
      cat(".")
      if (s==1) {
        avc <- array(dim=c(length(sims),length(sims[[1]]),length(snams)))
        ave <- avc
        av.qc <- array(0,dim=c(length(sims),length(probs),length(sims[[1]]),length(snams)))
        av.qe <- array(0,dim=c(length(sims),length(probs),length(sims[[1]]),length(snams)))
      }
      tr <- as.character(Di$D$R[as.logical(Di$D[[Di$SC[1]]])]) # correct
      for (i in 1:length(tr)) {        
        avc[i,,s] <- unlist(lapply(sims[[i]],function(x){100*mean(x$r==tr[[i]])}))
        av.qc[i,,,s] <- unlist(lapply(sims[[i]],function(x){
          do.call(quantile,list(x=x$rt[x$r==tr[[i]]],probs=probs,na.rm=TRUE))
        }))
      }
      tr <- as.character(Di$D$R[!as.logical(Di$D[[Di$SC[1]]])]) # error
      for (i in 1:length(tr)) {
        ave[i,,s] <- unlist(lapply(sims[[i]],function(x){100*mean(x$r==tr[[i]])}))
        av.qe[i,,,s] <- unlist(lapply(sims[[i]],function(x){
          do.call(quantile,list(x=x$rt[x$r==tr[[i]]],probs=probs,na.rm=TRUE))
        }))
      }
    }
    cat("\n")
    avc <- apply(avc,1:2,mean,na.rm=TRUE)
    ave <- apply(ave,1:2,mean,na.rm=TRUE)
    av.qc <- apply(av.qc,1:3,mean,na.rm=TRUE)
    av.qe <- apply(av.qe,1:3,mean,na.rm=TRUE)
    Di <- attr(rc$dat,"D")[[1]]
    D <- Di$D[as.logical(Di$D[[Di$SC[1]]]),]
    list(c=apply(apply(avc,2,function(x){
           as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),1,quantile,probs=pci,na.rm=TRUE),
         e=apply(apply(ave,2,function(x){
           as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),1,quantile,probs=pci,na.rm=TRUE),
        qc=aperm(apply(apply(av.qc,2:3,function(x){as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),
           1:2,quantile,probs=pci,na.rm=TRUE),c(1,3,2)),
        qe=aperm(apply(apply(av.qe,2:3,function(x){as.vector(tapply(x,D[,wsFs],mean,na.rm=TRUE))}),
           1:2,quantile,probs=pci,na.rm=TRUE),c(1,3,2))         
    )
  }
  
  # main body
  if (bar.type=="pbmodel") {
    bar.type="model"
    pb <- TRUE
  } else pb <- FALSE
  if (plot.subjects) bars <- NA
  probs <- sort(probs)
  require(lattice)
  dat <- rc[[1]]
  qp <- attr(dat,"qp")
  D <- attr(dat,"D")
  Di <- D[[1]]
  # Get factors
  if (any(is.na(factors))) wsF <- Di$F else wsF <- factors
  #   check(wsF,Di$F,data.fit.name)
  if ( plot.subjects ) {
    wsF <- c(wsF,Di$S)
    bars <- NA
  }
  wsF <- c(data.fit.name,wsF)
  plot.lines <- TRUE
  wsFs <- wsF[-1]
  if ( plot.subjects ) dwsFs <- wsFs else dwsFs <- c(wsFs,Di$S)
  # plotting formula
  plotform <- paste("y ~",wsFs[1])
  if ( length(wsFs)>1 ) if ( length(wsFs)==2 ) plotform <- paste(plotform,"|",wsFs[2]) else
    plotform <- paste(plotform,"|",paste(wsFs[-1],collapse="+"))
  # Get measure
  if (!is.null(qp) & !all(probs %in% qp))
      stop(paste("Quantile probs must be in set",paste(qp,collapse=","),"\n"))
  ylab <- "RT"
  cdat <- dat[as.logical(dat[,Di$SC[1]]),]
  edat <- dat[!as.logical(dat[,Di$SC[1]]),]
  datc <- fixholes(cdat,Di,"accuracy",plot.model)
  date <- fixholes(edat,Di,"inaccuracy",plot.model)
  if (!ylab.qs) ylab <- paste(ylab,"quantiles") else
    ylab <- paste(ylab,paste(probs,collapse=" ,"),"quantiles")
  snams <- levels(dat[[Di$S]])  
  if ( !any(is.na(exclude.subjects)) ) 
    snams <- snams[!(snams %in% exclude.subjects)]
  dat <- dat[dat[,Di$S] %in% snams,]
  if (!any(is.na(include.subjects))) {
    if (!all(include.subjects %in% snams))
      stop("Some specifed subjects not in set of subjects")
    cat(paste("Including only participants",paste(include.subjects,collapse=","),"\n"))  
    dat <- dat[dat[,Di$S] %in% include.subjects,]
    dat[,Di$S] <- factor(as.character(dat[,Di$S]))
  }
  pe <- getstat(date,"inaccuracy","mean",probs,dwsFs,Di)
  pc <- getstat(datc,"accuracy","mean",probs,dwsFs,Di)
  qe <- getstats(edat,"error",probs,dwsFs,Di,qnam=data.fit.name)
  qc <- getstats(cdat,"correct",probs,dwsFs,Di,qnam=data.fit.name)

#   badc <- apply(qc,2:length(dim(qc)),function(x){any(is.na(x))})
#   bade <- apply(qe,2:length(dim(qc)),function(x){any(is.na(x))})
#   if (all(bad))
#     stop("All participants have missing data cells")
#   if (exclude.incomplete && any(bad)) {
#     cat(paste("Excluding participants with missing data cells",
#               paste(snams[bad],collapse=","),"(",sum(!bad),"left )\n"))  
#     dat <- dat[dat[,Di$S] %in% snams[!bad],]
#     dat[,Di$S] <- factor(as.character(dat[,Di$S]))
#     data <- getstats(dat,measure,probs,dwsFs,Di,qnam=data.fit.name) 
#   }
#     data <- getstats(dat,measure,probs,dwsFs,Di,qnam=data.fit.name) 
#   }
  av <- arr2df(apply(qc,wsF,mean,na.rm=T))
  av[[wsFs[1]]] <- rep(arr2df(apply(pc,wsFs,mean,na.rm=T))$y,each=length(probs))
  tmp <- arr2df(apply(qe,wsF,mean,na.rm=T))
  tmp[[wsFs[1]]] <- rep(arr2df(apply(pe,wsFs,mean,na.rm=T))$y,each=length(probs))
  # First half correct quantiles, second error quantiles
  av <- rbind(av,tmp)
  xodr <- order(av[[wsFs[1]]])
#   av <- av[xodr,]
  if ( plot.model ) {
    pe.m <- getstat(date,"inaccuracy","mean",probs,dwsFs,Di,F,T)
    pc.m <- getstat(datc,"accuracy","mean",probs,dwsFs,Di,F,T)
    qe.m <- getstats(edat,"error",probs,dwsFs,Di,F,T,qnam=data.fit.name)
    qc.m <- getstats(cdat,"correct",probs,dwsFs,Di,F,T,qnam=data.fit.name)
    av.m <- arr2df(apply(qc.m,wsF,mean,na.rm=T))
    av.m[[wsFs[1]]] <- rep(arr2df(apply(pc.m,wsFs,mean,na.rm=T))$y,each=length(probs))
    tmp <- arr2df(apply(qe.m,wsF,mean,na.rm=T))
    tmp[[wsFs[1]]] <- rep(arr2df(apply(pe.m,wsFs,mean,na.rm=T))$y,each=length(probs))
    av.m <- rbind(av.m,tmp)
#     av.m <- av.m[xodr,]
    av <- rbind(av,av.m) 
    av[,1] <- paste(rep(c(data.name,model.name),each=dim(av.m)[1]),av[,1],sep=".")
    av[,1] <- factor(av[,1],unique(av[,1]))
    ltys <- rep(c(1:2),each=length(levels(av.m[,1])))
    pchs <- rep(c(16,1),each=length(levels(av.m[,1])))
  } else {
    ltys <- rep("blank",length(levels(av[[1]])))
    pchs <- 1
  }
  if ( !is.na(bars) ) {
    if ( bar.type=="model" & !plot.model )
      stop("Set plot.model=TRUE to enable model error-bar plotting")
    if ( !(any(bar.type %in% c("within","between","data","model"))) ) 
      stop("Parameter bar.type must be one of \"within\",\"between\", \"data\" or \"model\"")
    facs <- c(data.fit.name,wsFs)    
    if ( bar.type %in% c("within","between") ) {
      ns <- apply(qc.m,facs,function(x){sum(!is.na(x))})
      ns.p <- apply(qe.m,wsFs,function(x){sum(!is.na(x))})
      if (bar.type=="within") {
        m <- prod(unlist(lapply(dimnames(pe)[wsFs],length)))
        se.pc <- arr2df(apply(aperm(pe,c(Di$S,wsFs)) - apply(pe,Di$S,mean,na.rm=T),
                              wsFs,sd,na.rm=T)/sqrt(ns.p*(m-1)/m))
        m <- prod(unlist(lapply(dimnames(pc)[wsFs],length)))
        se.pe <- arr2df(apply(aperm(pc,c(Di$S,wsFs)) - apply(pc,Di$S,mean,na.rm=T),
                              wsFs,sd,na.rm=T)/sqrt(ns.p*(m-1)/m))
        m <- prod(unlist(lapply(dimnames(qe)[facs],length)))
        se.qe <- arr2df(apply(aperm(qe,c(Di$S,facs)) - apply(qe,Di$S,mean,na.rm=T),
                           facs,sd,na.rm=T)/sqrt(ns*(m-1)/m))
        m <- prod(unlist(lapply(dimnames(qc)[facs],length)))
        se.qc <- arr2df(apply(aperm(qc,c(Di$S,facs)) - apply(qc,Di$S,mean,na.rm=T),
                           facs,sd,na.rm=T)/sqrt(ns*(m-1)/m))
      } else {
        se.pe <- arr2df(apply(aperm(pe,c(Di$S,wsFs)),wsFs,sd,na.rm=T)/sqrt(ns.p))
        se.pc <- arr2df(apply(aperm(pc,c(Di$S,wsFs)),wsFs,sd,na.rm=T)/sqrt(ns.p))
        se.qe <- arr2df(apply(aperm(qe,c(Di$S,facs)),facs,sd,na.rm=T)/sqrt(ns))
        se.qc <- arr2df(apply(aperm(qc,c(Di$S,facs)),facs,sd,na.rm=T)/sqrt(ns))
      }
      if (is.numeric(bars)) { 
        se.pe$y <- -se.pe$y*qt((100-bars)/200,as.vector(ns.p))
        se.pc$y <- -se.pc$y*qt((100-bars)/200,as.vector(ns.p))
        se.qe$y <- -se.qe$y*qt((100-bars)/200,as.vector(ns))
        se.qc$y <- -se.qc$y*qt((100-bars)/200,as.vector(ns))
      }
      # error then correct
      se <- rbind(se.qc,se.qe)
      se.p <- rbind(se.pc,se.pe)
    } 
    if ( bar.type=="data" ) {
      stop("Sorry data resampling option not implemented")
    }
    if (bar.type=="model") LMU=model.bar(rc,wsFs,probs,pb)    
    av[,1] <- factor(av[,1],c(levels(factor(av[,1])),c("lo","hi")))
    if ( !plot.model ) lo=av else {
      is.dat <- unlist(lapply(strsplit(
        as.character(av$Type),".",fixed=T),function(x){x[1]}))==data.name
      lo <- av[is.dat,]
    }
    lo[,1] <- "lo"; hi <- lo; hi[,1] <- "hi"
    if (bar.type %in% c("within","between")) {
#       lo$y <- lo$y-se$y[xodr]
#       hi$y <- hi$y+se$y[xodr]
#       lo[[wsFs[1]]] <- lo[[wsFs[1]]]-rep(se.p$y,each=length(probs))[xodr]
#       hi[[wsFs[1]]] <- hi[[wsFs[1]]]+rep(se.p$y,each=length(probs))[xodr]
      lo$y <- lo$y-se$y
      hi$y <- hi$y+se$y
      lo[[wsFs[1]]] <- lo[[wsFs[1]]]-rep(se.p$y,each=length(probs))
      hi[[wsFs[1]]] <- hi[[wsFs[1]]]+rep(se.p$y,each=length(probs))
    } else {
      lo$y <- c(LMU$qc[1,,],LMU$qe[1,,])
      hi$y <- c(LMU$qc[3,,],LMU$qe[3,,]) 
      lo[[wsFs[1]]] <- c(rep(LMU$c[1,],each=length(probs)),
                         rep(LMU$e[1,],each=length(probs)))
      hi[[wsFs[1]]] <- c(rep(LMU$c[3,],each=length(probs)),
                         rep(LMU$e[3,],each=length(probs)))
      av[!is.dat,"y"] <- c(LMU$qc[2,,],LMU$qe[2,,])
      av[!is.dat,wsFs[1]] <- c(rep(LMU$c[2,],each=length(probs)),
                               rep(LMU$e[2,],each=length(probs)))
    }
    if (!plot.model) av <- av[xodr] else
      av <- av[c(xodr,length(xodr)+xodr),] 
    lo <- lo[xodr,]
    hi <- hi[xodr,]
    av <- rbind(av,lo,hi)    
  } else {
    if (!plot.model) av <- av[xodr] else
      av <- av[c(xodr,length(xodr)+xodr),]
  }
      
  cols <- "black"  
  if (plot.model) key=list(text=list( c(data.name,model.name) ),
        lines=list(lty=1:2),points=list(pch=c(16,1),cex=cex),columns=2)
  if (is.na(xlab))
    if (plot.lines) xlab <- wsF[2] else
      xlab <- wsF[1]  

  if (any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
      type=c("p","l"),ylab=ylab,layout=layout,as.table=T,col=cols,panel=panel.se,
      strip = strip.custom(bg="white"),lty=ltys)) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
          type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
          layout=layout,as.table=T,xlab=xlab,panel=panel.se,
          strip = strip.custom(bg="white")))
  } else if (!any(is.na(xlim)) & any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,panel=panel.se,
        type=c("p","l"),ylab=ylab,xlim=xlim,layout=layout,as.table=T,col=cols,
        strip = strip.custom(bg="white"),lty=ltys)) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,xlim=xlim,as.table=T,xlab=xlab,panel=panel.se,
        strip = strip.custom(bg="white")))    
  } else if (any(is.na(xlim)) & !any(is.na(ylim))) {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
        type=c("p","l"),ylab=ylab,ylim=ylim,as.table=T,col=cols,panel=panel.se,
        strip = strip.custom(bg="white"),lty=ltys)) else 
      print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
        type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
        layout=layout,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
        strip = strip.custom(bg="white")))    
  } else {
    if (!plot.model) print(xyplot(formula(plotform),av,groups=av[,wsF[1]],xlab=xlab,
      type=c("p","l"),ylab=ylab,xlim=xlim,ylim=ylim,as.table=T,col=cols,panel=panel.se,
      strip = strip.custom(bg="white"),lty=ltys)) else 
    print(xyplot(formula(plotform),av,groups=av[,wsF[1]],key=key,
      type=c("p","l"),ylab=ylab,pch=pchs,lty=ltys,col=cols,cex=cex,
      layout=layout,xlim=xlim,ylim=ylim,as.table=T,xlab=xlab,panel=panel.se,
      strip = strip.custom(bg="white")))    
  }
  if (save.dat) av
}

