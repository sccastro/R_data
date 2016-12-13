
rm(list=ls())
source("~/rc/utils.R")
source("~/rc/data.analysis.R")

# This is the dataset. It has 74063 rows total, 52 participants, 
# and I’ve put the columns into DMC format with
# s: subjects, 
# S: stimulus by s1 and s2 (S or H)
# F: accuracy and speed instruction by A and S
# R: response by r1 and r2
# RT: reaction time in seconds
# congruent: 1 for congruent and 0 for incongruent


dats  <- read.csv("FlankerLBA.csv")

dats$s <- factor(dats$s)

# # Remove two cases where S=="", convert levels to s1="S" and s2="H", reorder to
# # alphabetic H, S
# dats[dats$S=="",]
# #       s S F R RT congruent OverallAcc
# # 7811  2   A    0         1          0
# # 34212 5   A    0         1          0
dats <- dats[dats$S!="",]
dats$S <- factor(as.character(dats$S),levels=c("s2","s1"),labels=c("H","S"))

# # Some participants have 50 practice trials: remove and label Accuracy and Speed
# table(dats$F=="Practice",dats$s)
# #           1    2    3    4    5    6    7    8    9   10   11   12   13
# #   FALSE 2300 1910 1000 1200 1200 1200 1200 1100 1200 1100 1200 1100 1000
# #   TRUE    50   50   50   50   50   50   50   50   50   50    0   50   50
# #        
# #           14   15   16   17   18   19   20   21   22   23   24   25   26
# #   FALSE 1200  900 1200 1100 1200 1000 1200 1100 1200 1200 1100 1200 1100
# #   TRUE    50    0   50   50   50    0   50   50   50    0   50   50   50
# #        
# #           27   28   29   30   31   32   33   34   35   36   37   38   39
# #   FALSE 1200 1200 1200 1100 1100 1000 1200 1200 1200 1100 1200 1200 1200
# #   TRUE     0    0    0    0   50   50   50    0   50   50   50    0   50
# #        
# #           40   41   42   43   44   45   46   47   48   49   50   51   52
# #   FALSE 1100 1100 1300 1200 1100 1200 1100 1000 1100 1200 1200 1200 1100
# #   TRUE     0   50    0   50    0   50   50   50    0   50   50   50   50
# #        
# #           53   54   55   56   57   58   59   60   61
# #   FALSE 1200 1100 1100 1200 1200 1200 1000 1200 1200
# #   TRUE     0   50   50   50    0   50   50   50   50
dats <- dats[dats$F!="Practice",]
dats$F <- factor(as.character(dats$F),labels=c("Accuracy","Speed"))

#  great variety of responses, treat as NA all but r1 and r2,rename r1=S and r2=H
table(dats$R)
#                 {?}       {.} {CONTROL}        r1        r2         s         x 
#         4         1       102         2     35968     35710         1        22 
dats$R[!(dats$R %in% c("r1","r2"))] <- NA
dats$R <- factor(as.character(dats$R),levels=c("r2","r1"),labels=c("H","S"))
 
# # At worst 4% missed responses  
# round(sort(100*tapply(is.na(dats$R),dats$s,mean)))
# #   3  4  5  8  9 10 11 12 13 14 15 16 17 18 19 21 24 25 26 27 28 29 30 31 32 33 34 35 
# #  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
# # 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 54 55 56 57 59 60 61  6 20 23 53 58 
# #  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
# # 36  2  1  7 22 
# #  0  0  1  4  4
 
# Check accuracy scoring
dats$OverallAcc <- as.logical(dats$OverallAcc)
# tmp <- dats[!is.na(dats$R),]
# all(tmp$OverallAcc==(tmp$S==tmp$R))
# # [1] TRUE
# 
# # Looks like 25% incongruent design so expect big flanker effect!
# table(dats$congruent)
# #     0     1 
# # 17905 53905 
# # 17905/(17905+53905)
# # [1] 0.2493385

# Rename levels and reorder so alphabetic
dats$congruent <- factor(dats$congruent,levels=1:0,
  labels=c("Congruent","Incongruent"))

# Rename the way I like, note RS = relevant (central) stimulus
names(dats)[c(2,3,6,7)] <- c("RS","AS","CI","C")

# Figure out IS = irrelevant stimulus = flanker
dats$IS <- dats$RS # Congruent cases
# Incongruent case
is.H <- dats$IS=="S" & dats$CI=="Incongruent"
is.S <- dats$IS=="H" & dats$CI=="Incongruent"
dats$IS[is.H] <- "H"
dats$IS[is.S] <- "S"

# Remove non-responses
dat <- dats[!is.na(dats$R),]

# # Nice accuracy data
# tmp <- tapply(dat$C,dat[c("CI","s")],mean)
# round(100*tmp)
# #              s
# # CI             1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
# #   Congruent   94 98 97 98 99 93 94 95 98 98 96 88 97 99 96 94 96 92 99 99 96 98 99 98
# #   Incongruent 66 89 88 87 94 79 73 77 92 91 91 79 86 91 88 81 78 81 97 96 88 94 95 88
# #              s
# # CI            25 26 27 28 29 30 31 32 33 34 35 36  37 38 39 40 41 42 43 44 45 46 47
# #   Congruent   99 97 95 93 97 99 97 97 92 96 96 88 100 94 98 92 96 96 95 98 95 92 97
# #   Incongruent 96 96 82 78 90 99 82 81 85 90 85 84  87 83 92 79 89 72 82 93 88 88 89
# #              s
# # CI            48 49 50 51 52 53 54 55 56 57 58 59 60 61
# #   Congruent   97 99 99 97 98 98 97 93 95 99 97 99 96 97
# #   Incongruent 88 92 90 89 88 92 88 85 76 96 94 92 88 91

# # Interference effect
# -round(sort(100*apply(tmp,2,diff)))
# #  1 42  7 56  8 17 32 31 28  6 40 37 43 27 16 18 13 38 35  4 24 52 54  3 12 50  2 51 
# # 29 24 21 19 18 18 16 15 15 14 13 13 13 13 12 11 11 11 11 10 10  9  9  9  9  9  9  9 
# # 21 48 15 60 33 47 55 10 14 29 59 49 41 61 39 45 53 34 44  9 11  5 22 36 23 25 46 20 
# #  8  8  8  8  8  8  8  7  7  7  7  7  7  7  6  6  6  6  6  6  5  5  4  4  4  4  3  3 
# # 58 57 19 26 30 
# #  3  3  2  2  1

# # Nice accuracy data
# tmp <- tapply(dat$RT,dat[c("CI","s")],mean)
# round(1000*tmp)
# #              s
# # CI              1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
# #   Congruent   427 395 426 367 432 407 387 368 364 432 440 357 403 450 400 439 327 391
# #   Incongruent 483 439 471 414 465 473 456 418 405 486 473 404 462 534 424 486 369 420
# #              s
# # CI             19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
# #   Congruent   483 448 447 484 457 395 359 492 363 349 450 416 409 380 375 409 378 518
# #   Incongruent 554 538 459 546 495 450 393 578 390 380 527 508 453 418 437 462 405 560
# #              s
# # CI             37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
# #   Congruent   363 404 459 375 438 363 410 420 419 379 418 395 426 406 457 380 410 444
# #   Incongruent 395 449 506 437 489 403 442 467 445 418 466 442 458 442 520 436 475 502
# #              s
# # CI             55  56  57  58  59  60  61
# #   Congruent   397 406 481 449 436 381 442
# #   Incongruent 448 445 561 485 502 417 487

# # Interference effect
# round(sort(1000*apply(tmp,2,diff)))
# # 21 15 45 27 35 18 28 37 49 43  5 11 25 60 58 50 23 32 56 46 42  9 36 17  2 38 31 61 
# # 11 23 26 27 28 29 31 32 32 33 33 33 34 36 36 36 38 38 39 40 40 41 42 43 44 45 45 45 
# #  3 39 44 12 16  4 48 47  8 41 55 10 34 24 52  1 54 13 51 40 22 33 59  6 53  7 19 29 
# # 45 47 47 47 47 47 47 47 51 51 52 53 54 56 56 56 58 59 62 62 63 63 65 66 66 69 72 76 
# # 57 14 26 20 30 
# # 80 84 87 90 92 
 
# # Nice RT distribution
# hist(dat$RT,breaks="fd")
# 
# # Lower bound of .15s seems reasonable
# round(head(sort(dats$RT),300),3)
# #   [1] 0.000 0.000 0.000 0.000 0.001 0.001 0.001 0.001 0.002 0.002 0.003 0.003 0.005
# #  [14] 0.005 0.013 0.013 0.014 0.016 0.016 0.017 0.017 0.018 0.018 0.020 0.020 0.022
# #  [27] 0.026 0.028 0.028 0.029 0.029 0.029 0.029 0.030 0.030 0.031 0.036 0.036 0.036
# #  [40] 0.037 0.038 0.042 0.044 0.044 0.044 0.045 0.045 0.046 0.047 0.048 0.048 0.048
# #  [53] 0.049 0.051 0.051 0.053 0.053 0.054 0.054 0.055 0.056 0.056 0.058 0.058 0.058
# #  [66] 0.060 0.060 0.061 0.061 0.062 0.063 0.063 0.064 0.064 0.065 0.066 0.067 0.067
# #  [79] 0.068 0.068 0.068 0.068 0.069 0.070 0.072 0.072 0.072 0.073 0.073 0.074 0.074
# #  [92] 0.076 0.076 0.076 0.076 0.078 0.079 0.080 0.081 0.082 0.082 0.083 0.083 0.083
# # [105] 0.084 0.084 0.085 0.085 0.086 0.087 0.087 0.088 0.088 0.089 0.089 0.089 0.090
# # [118] 0.091 0.092 0.092 0.093 0.093 0.093 0.093 0.093 0.095 0.095 0.095 0.096 0.097
# # [131] 0.097 0.097 0.098 0.099 0.100 0.100 0.100 0.101 0.101 0.101 0.101 0.101 0.102
# # [144] 0.102 0.102 0.103 0.103 0.103 0.105 0.106 0.106 0.107 0.108 0.109 0.109 0.109
# # [157] 0.110 0.110 0.110 0.110 0.110 0.111 0.111 0.111 0.112 0.113 0.113 0.113 0.114
# # [170] 0.115 0.115 0.116 0.116 0.117 0.119 0.119 0.119 0.120 0.120 0.120 0.120 0.121
# # [183] 0.121 0.122 0.123 0.124 0.124 0.124 0.124 0.124 0.124 0.125 0.125 0.125 0.125
# # [196] 0.126 0.126 0.127 0.127 0.127 0.127 0.128 0.129 0.129 0.130 0.130 0.131 0.131
# # [209] 0.131 0.131 0.132 0.132 0.132 0.132 0.132 0.133 0.133 0.134 0.134 0.134 0.134
# # [222] 0.135 0.135 0.135 0.135 0.135 0.135 0.137 0.137 0.139 0.139 0.140 0.140 0.142
# # [235] 0.142 0.142 0.143 0.143 0.143 0.143 0.144 0.144 0.144 0.145 0.145 0.145 0.146
# # [248] 0.146 0.147 0.147 0.148 0.148 0.148 0.149 0.149 0.149 0.149 0.150 0.150 0.150
# # [261] 0.150 0.152 0.152 0.153 0.153 0.154 0.154 0.155 0.156 0.157 0.157 0.158 0.158
# # [274] 0.159 0.159 0.159 0.160 0.160 0.160 0.160 0.160 0.160 0.161 0.161 0.161 0.161
# # [287] 0.161 0.162 0.162 0.163 0.164 0.165 0.165 0.165 0.165 0.166 0.167 0.167 0.167
# # [300] 0.168

# # Upper bound of 1s seesm reasonable
# round(tail(sort(dats$RT),200),3)
# #   [1] 0.954 0.954 0.954 0.955 0.955 0.956 0.956 0.957 0.958 0.958 0.958 0.959 0.959
# #  [14] 0.961 0.964 0.964 0.964 0.966 0.967 0.970 0.970 0.970 0.971 0.971 0.972 0.972
# #  [27] 0.973 0.973 0.975 0.977 0.979 0.980 0.983 0.986 0.986 0.987 0.988 0.988 0.989
# #  [40] 0.991 0.991 0.992 0.993 0.995 0.996 0.997 0.997 0.998 0.998 1.000 1.000 1.003
# #  [53] 1.005 1.007 1.008 1.009 1.009 1.009 1.009 1.011 1.016 1.016 1.017 1.018 1.018
# #  [66] 1.024 1.024 1.026 1.026 1.028 1.028 1.029 1.029 1.029 1.030 1.030 1.033 1.034
# #  [79] 1.035 1.037 1.038 1.038 1.040 1.043 1.043 1.046 1.047 1.049 1.054 1.054 1.055
# #  [92] 1.056 1.057 1.060 1.068 1.068 1.069 1.071 1.073 1.079 1.081 1.082 1.082 1.083
# # [105] 1.083 1.088 1.088 1.092 1.094 1.094 1.111 1.111 1.116 1.116 1.119 1.120 1.122
# # [118] 1.122 1.123 1.139 1.140 1.140 1.141 1.151 1.155 1.159 1.162 1.164 1.168 1.168
# # [131] 1.170 1.175 1.181 1.186 1.186 1.192 1.192 1.195 1.196 1.196 1.198 1.199 1.207
# # [144] 1.209 1.214 1.216 1.219 1.224 1.229 1.238 1.248 1.250 1.260 1.267 1.281 1.284
# # [157] 1.285 1.286 1.289 1.298 1.300 1.319 1.322 1.324 1.343 1.343 1.345 1.349 1.360
# # [170] 1.361 1.376 1.389 1.409 1.432 1.449 1.459 1.471 1.473 1.478 1.488 1.492 1.493
# # [183] 1.510 1.521 1.534 1.553 1.555 1.578 1.597 1.637 1.681 1.744 1.775 1.836 1.863
# # [196] 1.887 1.895 1.908 1.932 1.938

# Format ready for RC 
source("~/rc/utils.R")
source("~/rc/data.analysis.R")
flanker.df <- score.rc(dat,S="s",R="R",RT="RT",F=c("RS","IS"),
  CV=c("CI"),SC="C")
# Spreading 71481 of 71678 RTs that are ties given preceision 0.001 .
#    839 have ties out of 1036 unique values
# 
# Added the following manifest design
#   IS RS R rcell
# 1  S  S S     1
# 2  S  S H     1
# 3  S  H S     2
# 4  S  H H     2
# 5  H  S S     3
# 6  H  S H     3
# 7  H  H S     4
# 8  H  H H     4

flanker=make.rc(flanker.df,correct.name="C",minrt=.15,maxrt=1) 
# Overall % RTs censored per subject
#  40  46  56  36  33   8  41  22  21   5   1  26  54  20  15  12  45  24  38  19  59 
# 7.1 4.4 3.4 3.3 2.6 1.5 1.1 1.0 0.8 0.8 0.7 0.6 0.6 0.6 0.6 0.5 0.5 0.5 0.4 0.4 0.4 
#  23  16  13  58   3  44  53  18  27  29  39  49  51   2  17   6  14  60   4   7   9 
# 0.3 0.3 0.3 0.3 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.1 0.0 0.0 0.0 
#  10  11  25  28  30  31  32  34  35  37  42  43  47  48  50  52  55  57  61 
# 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# Overall percent censored: 0.6 
# 
# Min/Median/Max number of correct responses per cell
# [1]  86 292 818
# 
# Min/Median/Max number of error responses per cell
# [1]   0  15 157

save(flanker,file="spencer.RData")

#######################################################################
rm(list=ls())
load("spencer.RData")
source("~/rc/utils.R")
source("~/rc/mt.R")
source("~/rc/AA-2n.R")
source("~/rc/MP-BNU2-B2v.R")           
source("~/rc/DF-BNU2.R") 
model.rc(dname="flanker",hmodelname="LBA",fitlist=flanker,
         Latents=list(lR=c("H","S")),Cpars=c("v","sv"),showmodels=F,
         Formulas=list(B="~lR*IS",A="~1",v="~RS*CI*C",sv="~C",ter="~1",pc="~1"),
         Bounds=list(B=c(0,Inf),A=c(0,Inf),v=c(-Inf,Inf),sv=c(0,Inf),ter=c(0,Inf),pc=c(0,1)),
         stopFormulas=list(v="~C"),
         Constants=list(sv=c(I=0),pc=c(I=qlogis(.00001))),
         saturated=T
         )

# TOP model has 15 parameters
# 
# Creating start point maps (slow for large model trees)
# 
# Creating model trees for subjects:
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 
# 
# Model 1 
# B~lR*IS & A~1 & v~RS*CI*C & sv~C & ter~1 & pc~1 
#    Start points from: 
#     B~IS & A~1 & v~RS*CI*C & sv~C & ter~1 & pc~1 
#     B~lR & A~1 & v~RS*CI*C & sv~C & ter~1 & pc~1 
#     B~lR*IS & A~1 & v~RS*C & sv~C & ter~1 & pc~1 
#     B~lR*IS & A~1 & v~CI*C & sv~C & ter~1 & pc~1 
#     B~lR*IS & A~1 & v~RS*CI*C & sv~1 & ter~1 & pc~1 
# 
# ...
# 
# Model 32 
# B~1 & A~1 & v~C & sv~1 & ter~1 & pc~1 
# No nested models
# 
# NUMBER OF FITS REQUIRED PER SUBJECT = 81 
# NUMBER OF FITS REQUIRED IN TOTAL = 4941 


### FOLLOWING NOT UPDATED YET

###############################################################

rm(list=ls())
source("~/rc/utils.R")
source("~/rc/data.analysis.R")

load(paste("spencer.LBA","spencer.RData",sep="/"))
lba <- spencer

##########  MODEL SELECTION

lba.top <- summary.rc(lba,showS=F,showbest=1,type="aic",show.TOP=F,output="pars",excludeS=excludeS)
summary.rc(lba,showS=F,showbest=1,type="bic",show.TOP=F,top.summary=F,excludeS=excludeS)

bestlba=toplba; lba.pars <- lba.top
bestlba=biclba; lba.pars <- lba.bic
bestlba=biclba1; lba.pars <- lba.bic1


# save.pars <- function(ps,type="LBA_BIC") {
#   pn <- names(ps)
#   for (i in pn[-length(pn)]) {
#     write.table(ps[[i]],sep="\t",quote=FALSE,row.names=FALSE,
#     file=paste(type,i,sep="."))
#   }
# }
# 
# save.pars(lba.top,"LBA_TOP")
# save.pars(ddm.top,"DDM_TOP")
# save.pars(lba.bic,"LBA_BIC")
# save.pars(ddm.bic,"DDM_BIC")
# # save.pars(lba.bic1,"LBA_BIC_noFter")


##### FITS

lba$dat <- get.simdat(dname="spencer.LBA",mname=bestlba)

# save.fits <- function(lba,type="LBA_TOP") {
#   rt <- lba$dat[is.finite(lba$dat$RT),]
#   p <- lba$dat[!is.finite(lba$dat$RT),]
#   write.table(rt[,-c(1,2,10:13,15)],sep="\t",quote=FALSE,row.names=FALSE,
#     file=paste(type,"RT",sep="."))
#   write.table(p[,-c(1,2,9,11:14)],sep="\t",quote=FALSE,row.names=FALSE,
#     file=paste(type,"qp",sep="."))
# }
# 
# save.fits(lba,"LBA_TOP")
# save.fits(lba,"LBA_BIC")

tmp=lba # LBA

plot.rc(tmp,measure="inaccuracy",plot.subjects=FALSE,bar.type="between",
        plot.model=T,model.name="LBA",bars=95,exclude.subjects=excludeS); #,factors=facs
plotqs.rc(tmp,measure="correct",probs=c(.1,.5,.9),plot.subjects=FALSE,bar.type="between",
          plot.model=T,model.name="LBA",bars=95,exclude.subjects=excludeS) # ,factors=facs
plotqs.rc(tmp,measure="error",probs=c(.1,.5,.9),plot.subjects=FALSE,bar.type="between",
          plot.model=T,model.name="LBA",bars=95,exclude.subjects=excludeS) # factors=facs,
plotce.rc(lba,plot.subjects=F,plot.model=T,model.name="LBA",bar.type="between",bars=95) # factors=facs,


plot.rc(tmp,measure="inaccuracy",plot.subjects=T,layout=c(6,3),
         plot.model=T,model.name="LBA",exclude.subjects=excludeS) # ,factors=facs
plotqs.rc(tmp,measure="correct",probs=c(.1,.5,.9),plot.subjects=T,layout=c(6,3),
          plot.model=T,model.name="LBA",exclude.subjects=excludeS) # ,factors=facs
plotqs.rc(tmp,measure="error",probs=c(.1,.5,.9),plot.subjects=T,layout=c(6,2),
          plot.model=T,model.name="LBA",exclude.subjects=excludeS) # ,factors=facs

######### 

plotpar.rc(lba,model=bestlba,parname="B",plot.subjects=F,bars="SE")
wsAnova(lba.pars$B)
mneffects(lba.pars$B,list("L"))

plotpar.rc(lba,model=bestlba,parname="ter",plot.subjects=F,bars="SE")
wsAnova(lba.pars$ter)
mneffects(lba.pars$ter,list("L"))


plotpar.rc(lba,model=bestlba,parname="v",plot.subjects=F,bars="SE",
  factors=c("C","S","P"))
wsAnova(lba.pars$v)


plotpar.rc(lba,model=bestlba,parname="sv",plot.subjects=F,bars="SE")
wsAnova(lba.pars$sv)


plotpar.rc(lba,model=bestlba,parname="A",plot.subjects=F,bars="SE")

