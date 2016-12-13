################################################################################
# LBA no transforms
################################################################################

  make_cnams <- function(pnams) {
    cnams=character(0)
    if (all(c("ab","bb","cb") %in% pnams)) {
      cnams <- c(cnams,c("ab","bb","cb"))
    } else cnams <- c(cnams,"b") 
    if (all(c("aa","ba","ca") %in% pnams)) {
      cnams <- c(cnams,c("aa","ba","ca"))
    } else cnams <- c(cnams,"a") 
    if (all(c("aV","bV","cV") %in% pnams)) {
      cnams <- c(cnams,c("aV","bV","cV"))
    } else cnams <- c(cnams,"v") 
    if (all(c("aSV","bSV","cSV") %in% pnams)) {
      cnams <- c(cnams,c("aSV","bSV","cSV"))
    } else cnams <- c(cnams,"sv") 
    if (all(c("aT","bT","cT") %in% pnams)) {
      cnams <- c(cnams,c("aT","bT","cT"))
    } else cnams <- c(cnams,"ter")
    cnams
  }  

  make_cols <- function(nacc,pnams,b,A,vC,sv,ter) {
    if (all(c("ab","bb","cb") %in% pnams)) {
      cc.cols <- cbind(rep(b,nacc),rep(1e-5,nacc),rep(1e-5,nacc))
    } else cc.cols <- matrix(rep(b,nacc),ncol=1) 
    if (all(c("aa","ba","ca") %in% pnams)) {
      cc.cols <- cbind(cc.cols,rep(A/b,nacc),rep(1e-5,nacc),rep(1e-5,nacc))  
    } else cc.cols <- cbind(cc.cols,rep(A/b,nacc)) 
    if (all(c("aV","bV","cV") %in% pnams)) {
      cc.cols <- cbind(cc.cols,rep(vC,nacc),rep(1e-5,nacc),rep(1e-5,nacc))  
    } else cc.cols <- cbind(cc.cols,rep(vC,nacc)) 
    if (all(c("aSV","bSV","cSV") %in% pnams)) {
      cc.cols <- cbind(cc.cols,rep(sv,nacc),rep(1e-5,nacc),rep(1e-5,nacc))  
    } else cc.cols <- cbind(cc.cols,rep(sv,nacc)) 
    if (all(c("aT","bT","cT") %in% pnams)) {
      cc.cols <- cbind(cc.cols,rep(ter,nacc),rep(1e-5,nacc),rep(1e-5,nacc))
    } else cc.cols <- cbind(cc.cols,rep(ter,nacc))
    cc.cols
  }

# trial by trial expansion of parlist
parlistt <- function(parlist,cv) {
  nt <- length(cv$T)
  parlist.t <- vector(mode="list")
  if (all(c("aT","bT","cT") %in% names(parlist))) {
    parlist.t$ter <- t(matrix(parlist$aT[1]+parlist$bT[1]*
        fun(rep(cv$T,each=length(parlist$aT)),parlist$cT[1]),ncol=nt))  
  } else parlist.t$ter <- matrix(rep(parlist$ter,each=nt),nrow=nt)
  if ( all(c("ab","bb","cb") %in% names(parlist)) ) {
    parlist.t$b <- t(matrix(parlist$ab + parlist$bb*
        fun(rep(cv$T,each=length(parlist$cb)),parlist$cb),ncol=nt))  
  } else parlist.t$b <- matrix(rep(parlist$b,each=nt),nrow=nt)
  if (all(c("aa","ba","ca") %in% names(parlist))) {
    parlist.t$A <- parlist.t$b * t(matrix(parlist$aa + parlist$ba*
        fun(rep(cv$T,each=length(parlist$ca)),parlist$ca),ncol=nt))         
  } else parlist.t$A <- matrix(rep(parlist$b*parlist$a,each=nt),nrow=nt) 
  if (all(c("aV","bV","cV") %in% names(parlist))) {
    parlist.t$v <- t(matrix(parlist$aV + parlist$bV*
        fun(rep(cv$T,each=length(parlist$cV)),parlist$cV),ncol=nt))      
  } else parlist.t$v <- matrix(rep(parlist$v,each=nt),nrow=nt)
  if (all(c("aSV","bSV","cSV") %in% names(parlist))) {
    parlist.t$sv <- t(matrix(parlist$aSV + parlist$bSV*
        fun(rep(cv$T,each=length(parlist$cSV)),parlist$cSV),ncol=nt))      
  } else parlist.t$sv <- matrix(rep(parlist$sv,each=nt),nrow=nt)
  parlist.t
}



start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,unitv=F,st0=.05)
}

M2P <- function(p) {
  # p is a data frame (made by start.model) with b, A, v, sv, ter 
  p
}


P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p
}

