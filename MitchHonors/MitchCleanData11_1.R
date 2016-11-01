# Clear/setup workspace ---------------------------------------------------

rm(list=ls()); par(mfrow = c(1,1))
getwd()
library(dplyr)
#Other users of this script will have to set their own working directory
setwd("C:/Users/u1010240/Desktop/Studies/R_data/MitchHonors")


# Get all files and apply functions to them -------------------------------

#Get filenames
files <- list.files(path="data", pattern="*.csv", full.names=T, recursive=FALSE)
#Get Subject IDs and conditions
subids <- regmatches(files, regexpr("\\d[a-z]?", files))
condition <- regmatches(subids, regexpr("[gmpw]", subids))

#List of files with all data
datalist = lapply(files, function(x){read.csv(file=x,header=T)})

#merge all of the files together
combined <- do.call(rbind, datalist)
summary(combined)



# Clean Data --------------------------------------------------------------


