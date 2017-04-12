# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_Data.R")
require(tidyverse)
require(gridExtra)
require(lubridate)
require(car)
require(lme4)
require(stringr)
library(forcats)

#Other users of this script will have to set their own working directory
setwd("~/Documents/R_data/MitchHonors/MitchChoice")


# 1. Get all files and apply functions to them -------------------------------
## 1.a Visual Data, or Vdata (This is on the touch screen)
vdirs <- list.dirs(path = "data", full.names = TRUE) #get directories 
vfiles <- list.files(path=vdirs, pattern="*.csv", full.names=T, recursive=FALSE) #Get filenames
vdatalist = lapply(vfiles, function(x){read.csv(file=x,header=T)}) #List of files with all data
vis.combined <- do.call(rbind, vdatalist) #merge all of the files together
vis.combined <- dplyr::rename(vis.combined, subid = sub) #rename the sub column subid (matches a function sub())

subids1 <- regmatches(vis.combined$subid, regexpr("2\\d{3}", vis.combined$subid)) #get subject ids
subids <- substr(subids1, 2,3) #get just the id number
condition1 <- regmatches(vis.combined$subid, regexpr("2\\d{3}", vis.combined$subid)) #get conditions
condition <- substr(condition1, 4,4) #get just the one number

vis.combined.choice <- vis.combined %>%
  mutate(subid = subids, condition = condition)
  


## 1.b Get all DRT files and apply functions to them
if(!exists("drt.combined.choice")) {
  drt.combined.choice <- ExtractDRT('DRTdata') #Only run this once or it will append to drt.data
}

##1.c Get all Steering data files and apply functions to them
if(!exists("steering.combined.choice")) {
  steering.combined.choice <- ExtractSteering('DRTdata') #Only run this once or it will append to drt.data
}


##1.d Save data frames for Visual and DRT
save(drt.combined.choice, vis.combined.choice, steering.combined.choice, file = "alldatachoice.Rdata")

# 2. clean.choice Data --------------------------------------------------------------
#####################################################################################
#                               NEW START HERE                                    #
#####################################################################################
rm(list=ls())
source("Load_Data.R")
load(file = "alldatachoice.Rdata")
factors <- c("subid", "condition", "numtargs")


## 2a. Visual Data

vis.combined.choice$subid <- factor(vis.combined.choice$subid)


clean.choice <- vis.combined.choice %>%
  filter(!str_detect(subid, 'pmix'), subid != "01", subid != "02", subid != "03") %>%
  dplyr::select(subid,time,condition,targcol,numtargs,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid) %>%
  mutate_each_(funs(factor), factors) %>%
  mutate(condition = fct_recode(condition, "practice" = "1", 
                                "green" = "2", 
                                "white" = "3", 
                                "mixed" = "4")) %>%
  filter((targcol == "green" & condition == "green" | condition == "mixed") |
           targcol == "white" & condition == "white" | condition == "mixed") %>%
  filter(rt1 > 150 & rt1 < 8000, rt2 > 150 & rt2 < 3000)

hist(clean.choice$rt1, breaks = 100)
hist(clean.choice$rt2, breaks = 100)


clean.choiceac <- vis.combined.choice %>%
  filter(!str_detect(subid, 'pmix')) %>%
  dplyr::select(subid,time,condition,targcol,numtargs,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid) %>%
  mutate_each_(funs(factor), factors) %>%
  mutate(condition = fct_recode(condition, "practice" = "1", 
                                "green" = "2", 
                                "white" = "3", 
                                "mixed" = "4",
                                "bad" = "5")) %>%
  filter((targcol == "green" & condition == "green" | condition == "mixed") |
           targcol == "white" & condition == "white" | condition == "mixed") %>%
  arrange(subid)

##2.b clean.choice drt.combined data 

### Get steering deviation

steering.combined.choice$cursorpos <- as.numeric(as.character(steering.combined.choice$cursorpos))

steer.dev <- steering.combined.choice %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos)) %>%
  dplyr::group_by(subid, condition) %>%
  dplyr::summarise(mean.dev = mean(deviation))


tapply(steer.dev$mean.dev,as.factor(steer.dev$condition),mean)
### combine steering with drt

drt.combined.choice <- dplyr::left_join(drt.combined.choice, steer.dev, by = c("subid","condition"))

drt.combined.choice$response <- recode_factor(drt.combined.choice$response, 
                                       wrong ="-1", correct ="1", miss = "0",
                                       .default = "NA") #refactor for getting percentages

drt.combined.choice$response <-as.numeric(as.character(drt.combined.choice$response))
drt.combined.choice <- drt.combined.choice[complete.cases(drt.combined.choice),]


drt.clean.choice <- drt.combined.choice %>%
  dplyr::select(subid = subid,time,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 & rt < 3000)

drt.ac <- drt.combined.choice %>%
  dplyr::select(subid = subid,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 | rt < 1, rt < 3000)

save(drt.combined.choice, vis.combined.choice, drt.clean.choice, 
     drt.ac, clean.choice, clean.choiceac, steer.dev, file = "cleanchoiceData.Rdata")


# Visual Analysis ---------------------------------------------------------

###Colors for plots
####Create a custom color scale
myColors <- c("forestgreen","grey50")
names(myColors) <- levels(clean.choice$targcol)
colScale <- scale_colour_manual(name = "targcol",values = myColors)


### Visual Search Plots

#### Accuracy
#### Percentage correct by group
acc.mean <- clean.choiceac %>%
  group_by(condition, targcol) %>%
  summarise(perc_corr = mean(corr)*100, n=n())


datac <- summarySEwithin(clean.choiceac, measurevar ="corr", withinvars = c("condition","targcol"), idvar = "subid")

accplot <- ggplot(datac, aes(x = condition, y=corr, fill = targcol)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=corr-ci, ymax=corr+ci)) +
  xlab("Condition") + coord_cartesian(ylim = c(.9, 1)) +
  scale_fill_manual(values=myColors) + ylab("Percent Correct (%)") +
  ggtitle("Accuracy of the \n Visual Search Task \n by Condition")
accplot

#### Reaction Time
rthist1 <- ggplot(data=clean.choice, aes(x=rt1, colour = targcol, linetype = condition)) + 
  stat_density(position="identity", geom="line", aes(group=interaction(targcol,condition))) +
  scale_linetype_manual(values = c("longdash","dotted","solid")) +
  theme_minimal() + colScale + theme(legend.position = "left") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
rthist1

rthist2 <- ggplot(data=clean.choice, aes(x=rt1, colour = as.factor(condition))) + 
  stat_density(position="identity", geom="line", aes(colour = as.factor(condition))) + 
  scale_colour_manual(name = "targcol",values = c("forestgreen","grey50","purple")) +
  theme_minimal() + ylab("") + scale_y_continuous(position = "right") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
rthist2

grid.arrange(rthist1, rthist2, bottom = "Reaction Time", 
             top = "Reaction Time by Condition in Selecting Target", 
             layout_matrix = matrix(c(1,2), ncol =2))


##### Violin plots of Reaction Time

rtplot1 <- ggplot(clean.choice, aes(x = condition, y=rt1, colour = targcol)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 1) + 
  stat_boxplot(geom ='errorbar', width = 0.1,position = position_dodge(.9))+
  geom_boxplot(width = 0.2, position = position_dodge(.9))+
  theme_minimal() + my.axis.font + xlab("") +
  colScale + theme(legend.position="none") + ylab("Reaction Time (ms)")
rtplot1

rtplot2 <- ggplot(clean.choice, aes(x = condition, y=rt1, color = condition)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_y_continuous(position = "right") +
  scale_colour_manual(name = "targcol",values = c("forestgreen","grey","purple"))+
  stat_boxplot(geom ='errorbar', width = 0.1)+
  geom_boxplot(width = 0.2)+
  theme_minimal() + my.axis.font  + xlab("") +
  theme(legend.position="none") + ylab("")
rtplot2

grid.arrange(rtplot1, rtplot2, bottom = "Condition", 
             top = "Reaction Time by Condition in Selecting Target", 
             layout_matrix = matrix(c(1,2), ncol =2))

##### ANOVA plots of Reaction Time
clean.choice.id <- clean.choice %>%
  group_by(condition,targcol) %>%
  summarise(mean.rt = mean(rt1))

clean.choice.id$new_c <- c("green.white","green.white","mixed","mixed")

main.plot <- ggplot(clean.choice.id, 
                    aes(x = new_c, y=mean.rt, colour = targcol, 
                        group= targcol, 
                        label = c("condition green", 
                                  "condition white", 
                                  "condition mixed", 
                                  "condition mixed"))) +
  scale_x_discrete(breaks = 0:1, labels = c("Single","Mixed")) +
  geom_line() + geom_point(size=5) + theme_minimal() +
  ggtitle("Reaction Time by Condition and Target Color \n for the DRT") + 
  scale_color_manual(name = "targcol",values = c("forestgreen","grey")) +
  geom_label(stat = "identity", nudge_y = -100) + ylab("Reaction Time") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ylim(1300, 3000) + xlab("Single or Mixed")

main.plot


# Stats -------------------------------------------------------------------

clean.choice.na <- lmer(rt1~condition*targcol+(1|subid),data=clean.choice)
Anova(clean.choice.na)


##2.b clean.choice drt.combined.choice data 

### Get steering deviation


steering.combined.choice$cursorpos <- as.numeric(as.character(steering.combined.choice$cursorpos))

steer.dev <- steering.combined.choice %>%
  filter(subid != "01", subid != "02", subid != "03") %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos)) %>%
  dplyr::group_by(subid, condition) %>%
  dplyr::summarise(mean.dev = mean(deviation))


tapply(steer.dev$mean.dev,as.factor(steer.dev$condition),mean)
### combine steering with drt

drt.combined.choice <- dplyr::left_join(drt.combined.choice, steer.dev, by = c("subid","condition"))

drt.combined.choice$response <- recode_factor(drt.combined.choice$response, 
                                       wrong ="-1", correct ="1", miss = "0",
                                       .default = "NA") #refactor for getting percentages

drt.combined.choice$response <-as.numeric(as.character(drt.combined.choice$response))
drt.combined.choice <- drt.combined.choice[complete.cases(drt.combined.choice),]


drt.clean.choice <- drt.combined.choice %>%
  dplyr::select(subid = subid,time,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  filter(subid != "01", subid != "02", subid != "03") %>%
  arrange(subid) %>%
  filter(rt > 150 & rt < 3000)

drt.ac <- drt.combined.choice %>%
  dplyr::select(subid = subid,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  filter(subid != "01", subid != "02", subid != "03") %>%
  arrange(subid) %>%
  filter(rt > 150 | rt < 1, rt < 3000)


## 3.b DRT data ------------------------------------------------------------
myColors2 <- c("red","grey","green","purple")
names(myColors2) <- levels(drt.clean.choice$condition)
colScale2 <- scale_colour_manual(name = "conditions",values = myColors2)

### RT Plots

rthist3 <- ggplot(data=drt.clean.choice, aes(x=rt, colour = as.factor(condition))) +
  stat_density(position = "identity", geom="line", aes(colour = as.factor(condition))) +
  ggtitle("Reaction Time \n by Condition \n in Selecting Target") + theme_minimal() +
  colScale2 + xlim(0,1000)
rthist3

####Miss Rate
####Percentage correct by group
drt.missac <- drt.ac %>%
  na.omit() %>%
  filter(response != -1) 

drt.dat.missac <- summarySEwithin(drt.missac, measurevar ="response", withinvars = c("condition"), idvar = "subid")

drt.missplot <- ggplot(drt.dat.missac, aes(x = as.factor(condition), y=response)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=response-ci, ymax=response+ci)) +
  xlab("Condition") + coord_cartesian(ylim = c(.5,1))
drt.missplot

miss.acc <- glmer(response~condition+(1|subid),data=drt.missac,family=binomial("probit"))
Anova(miss.acc)
#   Response: response
#             Chisq Df Pr(>Chisq)    
# condition 344.39  3  < 2.2e-16 ***

round(100*tapply(drt.missac$response,drt.missac[,"condition"],mean),1)
# Simple
#   1    2    3    4 
# 98.7 91.8 92.8 92.9 
# Choice
# 1    2    3    4 
# 86.4 57.5 62.1 59.2 

#### Incorrect Answers

drt.wrongac <- drt.ac %>%
  na.omit() %>%
  filter(response != 0) 

drt.wrongac$response <- factor(drt.wrongac$response)
drt.wrongac$response <- recode_factor(drt.wrongac$response,
                                      "-1" = "0", "1" ="1", "0" = "-1", 
                                      .default = "NA")
drt.wrongac$response <-as.numeric(as.character(drt.wrongac$response))

drt.dat.wrongac <- summarySEwithin(drt.wrongac, measurevar ="response", withinvars = c("condition"), idvar = "subid")

drt.accplot <- ggplot(drt.dat.wrongac, aes(x = as.factor(condition), y=response)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=response-ci, ymax=response+ci)) +
  xlab("Condition")
drt.accplot

wrong.acc <- glmer(response~condition+(1|subid),data=drt.wrongac,family=binomial("probit"))
Anova(wrong.acc)

# Simple
# Response: response
#             Chisq Df Pr(>Chisq)    
# condition 151.04  3  < 2.2e-16 ***
# Choice
# Response: response
# Chisq Df Pr(>Chisq)
# condition 1.1164  3     0.7731

# 3.c Steering Deviation --------------------------------------------------



steer.dat <- summarySEwithin(steer.dev, measurevar ="mean.dev", withinvars = c("condition"), idvar = "subid")

steer.plot <- ggplot(steer.dat, aes(x = as.factor(condition), y=mean.dev)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=mean.dev-ci, ymax=mean.dev+ci)) +
  xlab("Condition")
steer.plot


# 4. Hypothesis and Error Plots - Full Model ----------------------------------------------
library(heplots)
library(candisc)

clean.choice.mean <- clean.choice %>%
  group_by(subid,condition) %>%
  summarise(mean.rt = mean(rt1))

acc.mean <- clean.choiceac %>%
  group_by(subid, condition) %>%
  summarise(perc_corr = mean(corr)*100)

vis.mean <- dplyr::left_join(acc.mean,clean.choice.mean,by = c("subid","condition"))

drt.clean.choice$condition <- recode_factor(drt.clean.choice$condition, 
                                     "1" = "practice","2" = "green","3" ="white", "4" = "mixed",
                                     .default = "NA")
drt.ac$condition <- recode_factor(drt.ac$condition, 
                                  "1" = "practice","2" = "green","3" ="white", "4" = "mixed",
                                  .default = "NA")

drt.mean <- drt.clean.choice %>%
  group_by(subid,condition) %>%
  summarise(drt.mean = mean(rt), mean.dev = mean(mean.dev))

drt.meanac <- drt.ac %>%
  group_by(subid,condition) %>%
  summarise(drt.acc = mean(response)*100)



drt.meanall <- left_join(drt.mean,drt.meanac,by = c("subid","condition"))
drt.meanall$subid <- as.numeric(as.character(drt.meanall$subid))
drt.meanall <- drt.meanall %>%
  arrange(subid)
drt.meanall$subid <- factor(drt.meanall$subid,levels = 1:38)

all.combined <- left_join(drt.meanall, vis.mean, by = c("subid","condition"))
##DRT data and steering deviation

all.clean.choice <- all.combined %>%
  filter(condition != "practice", subid != 7)

all.clean.choice$condition <- factor(all.clean.choice$condition)

drt.model <- lm(cbind(drt.mean,mean.dev,drt.acc,perc_corr,mean.rt) ~ condition, data=all.clean.choice)
Anova(drt.model)

