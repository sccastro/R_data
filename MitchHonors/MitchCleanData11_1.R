# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("~/Documents/R_data/Nature_DRT/Load_Data.R")

library(dplyr)
library(ggplot2)

my.axis.font<-theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), 
                    plot.title=element_text(size=18,face="bold"),
                    legend.title=element_text(size=14),legend.text=element_text(size=14))



#Other users of this script will have to set their own working directory
setwd("~/Documents/R_data/MitchHonors")


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
combined <- dplyr::rename(combined, subid = sub)
# Clean Data --------------------------------------------------------------

#create a new column that finds the letter in the subid column and creates the condition.
combined$condition <- factor(ifelse(grepl("p", combined$subid, ignore.case = T), "practice", 
                  ifelse(grepl("w", combined$subid, ignore.case = T), "white", 
                         ifelse(grepl("g", combined$subid, ignore.case = T), "green","mixed"))))

#Select the columns we want, clean out NA's
clean <- combined %>%
  select(subid,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(!desc(subid)) %>%
  filter(rt1 > 300 & rt1 < 6000, rt2 > 150 & rt2 < 3000)

cleanac <- combined %>%
  select(subid,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(!desc(subid))
#Get the letters out of subid
cleanac$subid <- factor(gsub('[a-z]+', '', cleanac$subid))
cleanac$condition <- factor(cleanac$condition, levels = c("practice","white","green","mixed"))



#Get the letters out of subid
clean$subid <- factor(gsub('[a-z]+', '', clean$subid))

clean$condition <- factor(clean$condition, levels = c("practice","white","green","mixed"))


# Descriptive Stats -------------------------------------------------------
#Colors for plots
#Create a custom color scale
myColors <- c("red","grey","green","purple")
names(myColors) <- levels(clean$condition)
colScale <- scale_colour_manual(name = "conditions",values = myColors)


#Percentage correct by group
acc.mean <- clean %>%
  group_by(condition) %>%
  summarise(perc_corr = mean(corr)*100)

datac <- summarySEwithin(cleanac, measurevar ="corr", withinvars = c("condition"), idvar = "subid")


accplot <- ggplot(datac, aes(x = as.factor(condition), y=corr)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=corr-ci, ymax=corr+ci)) +
  xlab("Condition")
accplot



#RT
rthist1 <- ggplot(data=clean, aes(x=rt1, colour = as.factor(condition))) + 
  stat_density(position="identity", geom="line", aes(colour = as.factor(condition))) + 
  ggtitle("Reaction Time \n by Condition \n in Selecting Target") + theme_minimal() +
  colScale
rthist1

rthist2 <- ggplot(data=clean, aes(x=rt2, colour = as.factor(condition))) + 
  stat_density(position="identity", geom="line", aes(colour = as.factor(condition))) + 
  colScale +
  ggtitle("Reaction Time \n by Condition \n in Selecting Target") + theme_minimal()
rthist2

#Box and whisker plots

rtplot1 <- ggplot(clean, aes(x = as.factor(condition), y=rt1, colour = condition)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 1) + ylim(0,6000) + 
  stat_boxplot(geom ='errorbar', width = 0.1)+
  geom_boxplot(width = 0.2)+
  theme_minimal() + my.axis.font + xlab("Condition") +
  colScale + theme(legend.position="none")
rtplot1

rtplot2 <- ggplot(clean, aes(x = as.factor(condition), y=rt2, colour = condition)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,3000) +
  stat_boxplot(geom ='errorbar', width = 0.1)+
  geom_boxplot(width = 0.2)+
  theme_minimal() + my.axis.font + xlab("Condition") +
  theme(legend.position="none") + colScale
rtplot2



