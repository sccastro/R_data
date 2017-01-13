rm(list=ls())
source("rc/utils.R")
source("rc/data.analysis.R")

myData = read.csv("Spencer.csv")
myData$X = NULL
myData$s = factor(myData$s)
levels(myData$S)
myData$C <- as.numeric(myData$S) == as.numeric(myData$R)

myData <- myData[!is.na(myData$RT),]

##########################################################################

s4.df <- score.rc(myData,S="s",R="R",RT="RT",SC="C",F=c("S","F"),
                  autoscore=list(C=list(sfac="S",rfac="R")))

s4=make.rc(s4.df,correct.name="C",minrt=.25,maxrt=6) 

save(s4,file="spen.RData")

###########################################################

rm(list=ls())
load("spen.RData")
source("rc/utils.R")
source("rc/mt.R")
source("rc/AA-2n.R")
source("rc/MP-BNU2-B2v.R") 
source("rc/DF-BNU2p.R")

model.rc(dname="spen",hmodelname="LBA",fitlist=s4,
         Latents=list(lR=c("r1","r2")),Cpars=c("v","sv"),showmodels=F,
         Formulas=list(B="~F",A="~1",v="~S*C",sv="~S*C",ter="~1",pc="~1"),
         stopFormulas=list(B="~F",A="~1",v="~S*C",sv="~S*C",ter="~1",pc="~1"),
         Bounds=list(B=c(0,Inf),A=c(0,Inf),v=c(-Inf,Inf),sv=c(0,Inf),ter=c(0,Inf),pc=c(0,1)),
         Constants=list(sv=c(I=0),pc=c(I=qlogis(.00001))),
         saturated=T
)

# use terminal, change directory to "Spencer", then type in:
# nohup R CMD BATCH srim.R &

######################

rm(list=ls())
source("rc/utils.R")
source("rc/data.analysis.R")

load("spen.LBA/spen.RData")
lba <- spen
rm(spen)

lba.p <- summary.rc(lba,showbest=1,type="aic",output="pars")
summary.rc(lba,showbest=1,type="aic",top.summary=FALSE)

aiclba <- "B~F & A~1 & v~S*C & sv~S*C & ter~1 & pc~1"

lba.aicp <- summary.rc(lba,showbest=1,type="aic",output="pars",topnam=aiclba)

################################## Plots
bestlba <- aiclba

lba$dat <- get.simdat(dname="spen.LBA",mname=aiclba)

btype="within"; bars=95
tmp <- lba; mnam <- "LBA"

facs=c("F", "S")

require(ggplot2)

plot.rc(tmp,measure="accuracy",plot.subjects=F,factors=facs, xlab="Load",
        plot.model=T,model.name=mnam,bar.type=btype,bars=bars)

plotqs.rc(tmp,measure="correct",probs=c(.1,.5,.9),plot.subjects=F,factors=facs, xlab="Load", 
          plot.model=T,model.name=mnam,bar.type=btype,bars=bars)

plotqs.rc(tmp,measure="error",probs=c(.1,.5,.9),plot.subjects=F,factors=facs, xlab="Load", 
          plot.model=T,model.name=mnam,bar.type=btype,bars=bars)


# Threshold ---------------------------------------------------------------
plotpar.rc(lba,bestlba,"B",plot.subjects=F,bars="se", xlab="Load", factors=c("F"), ylab="Threshold (B)")

round(tapply(lba.aicp$B$y, lba.aicp$B$F, mean), 2)

# Drift ---------------------------------------------------------------
plotpar.rc(lba,bestlba,"v",plot.subjects=F,bars="se",
           xlab="Match Factor",factors=c("S", "C"), ylab="Mean Drift Rate (v)")

round(tapply(lba.aicp$v$y, lba.aicp$v$S, mean), 2)
round(tapply(lba.aicp$v$y, lba.aicp$v$C, mean), 2)
