################################################################################
# Original LBA B = b-A, v(Correct)=1-v(Error), sv free parameterization
################################################################################

start.model <- function(dati,pnams,Bounds) {
  sfun(dati,pnams,pc=.01,aprior=.5,tertune=.9,Atune=.5,unitv=T,st0=.05)
}

M2P <- function(p,D) {
  # p is a data frame (made by start.model) with b, A, v, sv, ter 
  p$B <- p$b-p$A
  p$vc[as.logical(D$D[[D$SC[1]]])] <- p$v[as.logical(D$D[[D$SC[1]]])]
  p$vc[!as.logical(D$D[[D$SC[1]]])] <- p$v[as.logical(D$D[[D$SC[1]]])]

  p$ve[!as.logical(D$D[[D$SC[1]]])] <- p$v[!as.logical(D$D[[D$SC[1]]])]  
  p$ve[as.logical(D$D[[D$SC[1]]])] <- p$v[!as.logical(D$D[[D$SC[1]]])]  

  p$svc[as.logical(D$D[[D$SC[1]]])] <- p$sv[as.logical(D$D[[D$SC[1]]])]
  p$svc[!as.logical(D$D[[D$SC[1]]])] <- p$sv[as.logical(D$D[[D$SC[1]]])]

  p$sve[!as.logical(D$D[[D$SC[1]]])] <- p$sv[!as.logical(D$D[[D$SC[1]]])]  
  p$sve[as.logical(D$D[[D$SC[1]]])] <- p$sv[!as.logical(D$D[[D$SC[1]]])]   
  
  p
}

P2Mfit <- function(p,pmap,D) {
  p
}

P2Mnatural <- function(p,pmap,D) {
  # changes to parameters used by density, CDF and ML etc.
  p$b = p$A + p$B
  p$v[as.logical(D$D[[D$SC[1]]])] <- p$vc[as.logical(D$D[[D$SC[1]]])]
  p$v[!as.logical(D$D[[D$SC[1]]])] <- p$ve[!as.logical(D$D[[D$SC[1]]])]  
  p$sv[as.logical(D$D[[D$SC[1]]])] <- p$svc[as.logical(D$D[[D$SC[1]]])]
  p$sv[!as.logical(D$D[[D$SC[1]]])] <- p$sve[!as.logical(D$D[[D$SC[1]]])]  
  p
}


