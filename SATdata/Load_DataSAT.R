#######################Spencer's Data Manipulation Functions######################################

rm(list = ls())
require(doBy)
require(plyr)
require(tidyverse)
require(lme4)
require(ggplot2)
require(TeachingDemos)
library(splitstackshape)
library(tcltk)
# My plot config ----------------------------------------------------------

my.axis.font<-theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), 
                    plot.title=element_text(size=18,face="bold"),
                    legend.title=element_text(size=14),legend.text=element_text(size=14))



# Functions ---------------------------------------------------------------

# Get Value

getValue <- function(x, data) {
  tmp <- data %>%
    filter(between(x, tlow, tup))
  return(tmp$targcol)
}

#Combine DRT files
ExtractDRT <- function(topdir) {
  filepath <- getwd()
  setwd(filepath)
  
  dirs <- list.dirs(path = '4_3_2017', full.names = TRUE) #get directories
  files <- list.files(path=dirs, pattern="*.txt", full.names=T, recursive=F) #get files
  subids1 <- regmatches(files, regexpr("3\\d{3}_", files)) #get subject ids
  subids <- substr(subids1, 2,3) #get just the id number
  condition1 <- as.character(lapply(strsplit(files, "_"), "[",4)) #get conditions
  condition <- substr(condition1, 4,4) #get just the one number
  
  
  pb <- tkProgressBar(title = "progress bar", min = 0,
                      max = length(files), width = 300)
  
  for (i in 1:length(files)) {
    drts <- readLines(files[i])
    matches <- regexpr("^DRT.*$", drts)
    drts <- regmatches(drts, matches)
    drts <- drts[drts != ""]
    dat <- cSplit(data.frame(drts),"drts", ",")
    dat <- dplyr::select(dat, time = drts_3, rt = drts_4, s2 = drts_5, s1 = drts_6, R = drts_7, response = drts_8)
    infile <- dat
    infile$subid <- as.factor(subids[i])
    infile$condition <- as.factor(condition[i])
    
    Sys.sleep(0.1)
    setTkProgressBar(pb, i, label=paste( round(i/length(files)*100, 0),
                                         "% done"))
    
    if(!exists("drt.data")) {
      drt.data <- infile
    }
    else {
      drt.data <- rbind(drt.data,infile)
    }
  }
  close(pb)
  return(drt.data)
}

#Combine Steering data files

ExtractSteering <- function(topdir) {
  filepath <- getwd()
  setwd(filepath)
  
  dirs <- list.dirs(path = topdir, full.names = TRUE) #get directories
  files <- list.files(path=dirs, pattern="*.txt", full.names=T, recursive=F) #get files
  subids1 <- regmatches(files, regexpr("3\\d{3}_", files)) #get subject ids
  subids <- substr(subids1, 2,3) #get just the id number
  condition1 <- as.character(lapply(strsplit(files, "_"), "[",4)) #get conditions
  condition <- substr(condition1, 4,4) #get just the one number
  
  pb <- tkProgressBar(title = "progress bar", min = 0,
                      max = length(files), width = 300)
  
  for (i in 1:length(files)) {
    drts <- readLines(files[i])
    matches <- regexpr("^Follow.*$", drts)
    drts <- regmatches(drts, matches)
    dat <- cSplit(data.frame(drts),"drts", ",")
    dat <- dplyr::select(dat, ballpos = drts_3, cursorpos = drts_4)
    infile <- dat
    infile$subid <- as.factor(subids[i])
    infile$condition <- as.factor(condition[i])
    
    Sys.sleep(0.1)
    setTkProgressBar(pb, i, label=paste( round(i/length(files)*100, 0),
                                         "% done"))
    
    if(!exists("drt.data")) {
      drt.data <- infile
    }
    else {
      drt.data <- rbind(drt.data,infile)
    }
  }
  close(pb)
  return(drt.data)
}

#Combine .csv files

CombineData <- function(dir) {
  
  filepath <- getwd()
  setwd(filepath)
  
  dirs <- list.dirs(path = dir, full.names = TRUE) #get directories
  files <- list.files(path=dirs, pattern="*.csv", full.names=T, recursive=F)
  subids <- regmatches(files, regexpr("1\\d{3}", files))
  condition1 <- as.character(lapply(strsplit(files, "_"), "[",4))
  condition <- substr(condition1, 4,4)
  
  
  for (i in 1:length(files)) {
    infile <- read.csv(files[i])
    infile$subids <- as.factor(subids[i])
    infile$condition <- as.factor(condition[i])
    
    if(!exists("combined.data")) {
      combined.data <- infile
    }
    else {
      combined.data <- rbind(combined.data,infile)
    }
  }
  return(combined.data)
}

# Clean up dataset --------------------------------------------------------

CleanData <- function(combined.data) {
  clean.data <- combined.data %>%
    select(subid = subids, condition, isPractice, pressCount,
           r1g0 = LED.Red..1.Red.0.Green., rt = Response.Time)
  
  levels(clean.data$condition) <- list(PreTesting = ("PreTesting"), 
                                       Testing = ("Testing"),
                                       Posttesting = c("Posttesting","Postsetting","Posttetting"))
  clean.data$isPractice <- factor(clean.data$isPractice)
  levels(clean.data$isPractice) <- c("no","yes")
  clean.data$r1g0 <- factor(clean.data$r1g0)
  levels(clean.data$r1g0) <- c("g", "r")
  clean.data$subid <- factor(clean.data$subid)
  
  
  
  return(clean.data) 
}


# Get Normalized Standard Errors for One Subject --------------------------

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# Error Bars for Repeated Measures ----------------------------------------

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}



# Normalize Data Within Groups --------------------------------------------

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

# Calculate D' -------------------------------------------------------

dprime = function(data) {
  yes         = subset(data, resp=="Y")
  no          = subset(data, resp=="N")
  hit         = subset(data, resp=="Y" & acc == 1)
  falsealarm  = subset(data, resp=="Y" & acc == 0)
  
  Hrate = xtabs(~subject, data=hit)/xtabs(~subject, data=yes)
  Frate = xtabs(~subject, data=falsealarm)/xtabs(~subject, data=no)
  dprime_score = qnorm(Hrate) - qnorm(Frate)
  
  return(dprime_score)
}


# Calculate A' ------------------------------------------------------------



# Standard Error ----------------------------------------------------------

SE <- function(dependent_variable) {
  sd(dependent_variable)/sqrt(length(dependent_variable))
}


# Confidence Intervals ----------------------------------------------------


lower <- function(dependent_variable) {
  mean(dependent_variable)-2*SE(dependent_variable)
}
upper <- function(dependent_variable) {
  mean(dependent_variable)+2*SE(dependent_variable)
}


# Hadley Counter Function -------------------------------------------------
init.counter <- function(){
  x <- 0
  function(){
    x <<- x + 1
    x
  }
}  #source: hadley wickham

counter <- init.counter()

# Trial Counter Function --------------------------------------------------



trial <- rep(1:5, 4)
block <- 0
blocklist <- 0

for (i in seq_along(df$trial)){
  if (df$trial[i]==1){
    block = block + 1}else
 if (df$trial!=1){
   block = block}
   blocklist<- c(blocklist, block)
}

blocklist <- blocklist[-1]
df$block <- blocklist
