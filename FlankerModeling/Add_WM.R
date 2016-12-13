rm(list = ls())
require(dplyr)


# Read Data ---------------------------------------------------------------
demos <-read.csv("~/Documents/R_data/FlankerModeling/PhaseTwoParticipantTrackingSheet.csv")
df <- read.csv("~/Documents/R_data/FlankerModeling/FlankerLBA.csv")


# Get subids and filter bad data ------------------------------------------
demosclean <- demos %>%
  select(s = X, gdata = Good.Data.,WMC = H.L.Span) %>%
  filter(gdata == "yes") #original file has notes as to which participants to drop

gsubids <- demosclean$s #Make a list of subids
length(gsubids) #56 good participants

dfclean <- df %>%
  arrange(s) %>%
  filter(s %in% gsubids) #Get only the good participants


df.full <- merge(dfclean,demosclean, by = "s",all = TRUE) #Merge datasets

l1 <- c("s","S","F", "R", "congruent","WMC") #Make sure factors are factors list
l2 <- c("RT", "OverallAcc") #Make sure numbers are numbers list

df.fc <- df.full %>%
  select(-gdata)%>% #Get rid of good data column
  mutate_each_(funs(factor), l1) %>% mutate_each_(funs(as.numeric), l2) %>% #factors and numbers are right.
  filter(S != "", F != "Practice", R == "r2" | R == "r1") #Get rid of extra levels with weird names/practice

#drop unused levels
df.fc$R <- factor(df.fc$R)
df.fc$F <- factor(df.fc$F)
df.fc$S <- factor(df.fc$S)

write.csv(x=df.fc, file = "FlankerDataClean.csv")

print(unique(subset(df.fc, WMC == "L")[,1]))

print(unique(subset(df.fc, WMC != "L")[,1]))
