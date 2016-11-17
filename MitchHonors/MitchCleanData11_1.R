# Clear/setup workspace ---------------------------------------------------
my.axis.font<-theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), 
                    plot.title=element_text(size=18,face="bold"),
                    legend.title=element_text(size=14),legend.text=element_text(size=14))

rm(list=ls()); par(mfrow = c(1,1))

library(dplyr)
library(ggplot2)
#Other users of this script will have to set their own working directory
setwd("C:/Users/u1010240/Desktop/Studies/R_data/MitchHonors")


# Get all files and apply functions to them -------------------------------

#Get filenames
files <- list.files(path="data", pattern="*.csv", full.names=T, recursive=FALSE)

#Get Subject IDs and conditions
# subids <- regmatches(files, regexpr("\\d[a-z]?", files))
# condition <- regmatches(subids, regexpr("[gmpw]", subids))

#List of files with all data
datalist = lapply(files, function(x){read.csv(file=x,header=T)})

#merge all of the files together
combined <- do.call(rbind, datalist)


#rename the sub column subid (matches a function sub())
combined <- rename(combined, subid = sub)
# Clean Data --------------------------------------------------------------

#create a new column that finds the letter in the subid column and creates the condition.
combined$condition <- factor(ifelse(grepl("p", combined$subid, ignore.case = T), "practice", 
                  ifelse(grepl("w", combined$subid, ignore.case = T), "white", 
                         ifelse(grepl("g", combined$subid, ignore.case = T), "green","mixed"))))

#Select the columns we want, clean out NA's
clean <- combined %>%
  select(subid,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(!desc(subid))

#Get the letters out of subid
clean$subid <- factor(gsub('[a-z]+', '', clean$subid))



# Descriptive Stats -------------------------------------------------------

#Percentage correct by group
acc.mean <- clean %>%
  group_by(condition) %>%
  summarise(perc_corr = mean(corr)*100)

#RT


accplot <- ggplot(acc.mean, aes(x = as.factor(condition), y=perc_corr, colour = perc_corr)) +
  geom_point(size=5) + ylim(90,100) + theme_minimal() +  my.axis.font +
  xlab("Condition")
accplot

rtplot <- ggplot(clean, aes(x = as.factor(condition), y=rt1, colour = condition)) + 
  geom_boxplot() + ylim(0,1300) + theme_minimal() + my.axis.font + xlab("Condition") 
# + guides(colour=FALSE) +
# scale_x_discrete(breaks = 1:2, labels = c("Green","White"))
#can add target color --(colour = targcol, group=targcol)
rtplot
