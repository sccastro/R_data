# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_Data.R")
require(tidyverse)
require(gridExtra)
require(lubridate)
require(car)


#Other users of this script will have to set their own working directory
setwd("~/Documents/R_data/MitchHonors")


# 1. Get all files and apply functions to them -------------------------------
## 1.a Visual Data, or Vdata (This is on the touch screen)
vdirs <- list.dirs(path = "data", full.names = TRUE) #get directories 
vfiles <- list.files(path=vdirs, pattern="*.csv", full.names=T, recursive=FALSE) #Get filenames
vdatalist = lapply(vfiles, function(x){read.csv(file=x,header=T)}) #List of files with all data
vis.combined <- do.call(rbind, vdatalist) #merge all of the files together
vis.combined <- dplyr::rename(vis.combined, subid = sub) #rename the sub column subid (matches a function sub())



## 1.b Get all DRT files and apply functions to them
if(!exists("drt.combined")) {
  drt.combined <- ExtractDRT('DRTdata') #Only run this once or it will append to drt.data
}

##1.c Get all Steering data files and apply functions to them
if(!exists("steering.combined")) {
  steering.combined <- ExtractSteering('DRTdata') #Only run this once or it will append to drt.data
}


##1.d Save data frames for Visual and DRT
save(drt.combined, vis.combined, steering.combined, file = "alldata.Rdata")

# 2. Clean Data --------------------------------------------------------------
#####################################################################################
#                               MITCH START HERE                                    #
#####################################################################################
rm(list=ls())
source("Load_Data.R")
load(file = "alldata.Rdata")
#2.a Visual Data, or Vdata cleaned first
###create a new column that finds the letter in the subid column and creates the condition.
vis.combined$condition <- factor(ifelse(grepl("p", vis.combined$subid, ignore.case = T), "practice", 
                  ifelse(grepl("w", vis.combined$subid, ignore.case = T), "white", 
                         ifelse(grepl("g", vis.combined$subid, ignore.case = T), "green","mixed"))))

vis.combined$subid <- factor(gsub('[a-z]+', '', vis.combined$subid), levels = 1:38)
vis.combined$condition <- factor(vis.combined$condition, levels = c("practice","green","white","mixed"))

vis.combined %>%
  dplyr::select(subid,time,numtargs,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid)

  
  # mutate(UTC = paste(year(UTC), month(UTC), mday(UTC), hour(UTC), sep="-"))

###Select the columns we want, clean out NA's
clean <- vis.combined %>%
  dplyr::select(subid,time,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter((targcol == "green" & condition == "green" | condition == "mixed") |
           targcol == "white" & condition == "white" | condition == "mixed") %>%
  filter(rt1 > 150 & rt1 < 6000, rt2 > 150 & rt2 < 3000)
  

hist(clean$time, breaks = 100)

cleanac <- vis.combined %>%
  dplyr::select(subid,time,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  filter((targcol == "green" & condition == "green" | condition == "mixed") |
           targcol == "white" & condition == "white" | condition == "mixed") %>%
  arrange(subid)


##2.b Clean drt.combined data 

### Get steering deviation

steering.combined$cursorpos <- as.numeric(as.character(steering.combined$cursorpos))

steer.dev <- steering.combined %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos)) %>%
  dplyr::group_by(subid, condition) %>%
  dplyr::summarise(mean.dev = mean(deviation))


tapply(steer.dev$mean.dev,as.factor(steer.dev$condition),mean)
### combine steering with drt

drt.combined <- dplyr::left_join(drt.combined, steer.dev, by = c("subid","condition"))

drt.combined$response <- recode_factor(drt.combined$response, 
                                       wrong ="-1", correct ="1", miss = "0",
                                       .default = "NA") #refactor for getting percentages

drt.combined$response <-as.numeric(as.character(drt.combined$response))
drt.combined <- drt.combined[complete.cases(drt.combined),]


drt.clean <- drt.combined %>%
  dplyr::select(subid = subid,time,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 & rt < 3000)

drt.ac <- drt.combined %>%
  dplyr::select(subid = subid,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 | rt < 1, rt < 3000)

save(drt.combined, vis.combined, drt.clean, drt.ac, clean, cleanac, steer.dev, file = "CleanData.Rdata")

#####################################################################################
#                                   New Start                                       #
#####################################################################################
# 3. Descriptive Stats -------------------------------------------------------
## 3.a Visual Plots first
rm(list = ls())
source("Load_Data.R")
load("CleanData.Rdata")
###Colors for plots
####Create a custom color scale
myColors <- c("forestgreen","grey50")
names(myColors) <- levels(clean$targcol)
colScale <- scale_colour_manual(name = "targcol",values = myColors)


### Visual Search Plots

#### Accuracy
#### Percentage correct by group
acc.mean <- cleanac %>%
  group_by(condition, targcol) %>%
  summarise(perc_corr = mean(corr)*100, n=n())


datac <- summarySEwithin(cleanac, measurevar ="corr", withinvars = c("condition","targcol"), idvar = "subid")

datac$new_c <- c("green.white","mixed","mixed","green.white")

accplot <- ggplot(datac, aes(x = new_c, y=corr, 
                             color = targcol, group = targcol, 
                             label = c("condition green",
                                       "condition mixed",
                                       "condition mixed",
                                       "condition white"))) +
  scale_x_discrete("Single or Mixed", breaks = 0:1, 
                   labels = c("Single","Mixed")) +
  geom_point(stat = "identity") + theme_minimal() + geom_line() +
  scale_y_continuous(position = "right") +
  geom_errorbar(width = .25, aes(ymin=corr-ci, ymax=corr+ci)) +
  geom_label_repel(stat = "identity", nudge_y = -.02, segment.alpha = .01) +
  coord_cartesian(ylim = c(.9, 1)) +
  scale_color_manual(values=myColors) + ylab("Percent Correct (%)") +
  ggtitle("Accuracy of the \n Visual Search Task \n by Condition") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5))

accplot

#### Reaction Time
rthist1 <- ggplot(data=clean, aes(x=rt1, colour = targcol, linetype = condition)) + 
  stat_density(position="identity", geom="line", aes(group=interaction(targcol,condition))) +
  scale_linetype_manual(values = c("longdash","dotted","solid")) +
  theme_minimal() + colScale + theme(legend.position = "left") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

rthist2 <- ggplot(data=clean, aes(x=rt1, colour = as.factor(condition))) + 
  stat_density(position="identity", geom="line", aes(colour = as.factor(condition))) + 
  scale_colour_manual(name = "targcol",values = c("forestgreen","grey50","purple")) +
  theme_minimal() + ylab("") + scale_y_continuous(position = "right") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

grid.arrange(rthist1, rthist2, bottom = "Reaction Time", 
             top = "Reaction Time by Condition in Selecting Target", 
             layout_matrix = matrix(c(1,2), ncol =2))


##### Violin plots of Reaction Time

# rtplot1 <- ggplot(clean, aes(x = condition, y=rt1, colour = targcol)) + 
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 1) + 
#   stat_boxplot(geom ='errorbar', width = 0.1,position = position_dodge(.9))+
#   geom_boxplot(width = 0.2, position = position_dodge(.9))+
#   theme_minimal() + my.axis.font + xlab("") +
#   colScale + theme(legend.position="none") + ylab("Reaction Time (ms)")
# 
# 
# rtplot2 <- ggplot(clean, aes(x = condition, y=rt1, color = condition)) + 
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
#   scale_y_continuous(position = "right") +
#   scale_colour_manual(name = "targcol",values = c("forestgreen","grey","purple"))+
#   stat_boxplot(geom ='errorbar', width = 0.1)+
#   geom_boxplot(width = 0.2)+
#   theme_minimal() + my.axis.font  + xlab("") +
#   theme(legend.position="none") + ylab("")
# 
# 
# grid.arrange(rtplot1, rtplot2, bottom = "Condition", 
#              top = "Reaction Time by Condition in Selecting Target", 
#              layout_matrix = matrix(c(1,2), ncol =2))

##### ANOVA plots of Reaction Time
clean.id <- clean %>%
  group_by(condition,targcol) %>%
  summarise(mean.rt = mean(rt1))

datrt <- summarySEwithin(clean, measurevar ="rt1", withinvars = c("condition","targcol"), idvar = "subid")

datrt$new_c <- c("green.white","mixed","mixed","green.white")

main.plot <- ggplot(datrt, 
                    aes(x = new_c, y=rt1, colour = targcol, 
                        group= targcol, 
                        label = c("condition green", 
                                  "condition white", 
                                  "condition mixed", 
                                  "condition mixed"))) +
  scale_x_discrete(breaks = 0:1, labels = c("Single","Mixed")) +
  geom_line() + geom_point(size=1) + theme_minimal() +
  geom_errorbar(width = .25, aes(ymin = rt1 - ci, ymax = rt1 + ci)) +
  ggtitle("Reaction Time by Condition and \n Target Color for the DRT") + 
  scale_color_manual(name = "targcol",values = c("forestgreen","grey")) +
  geom_label(stat = "identity", nudge_y = -150) + ylab("Reaction Time") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  ylim(1300, 3000) + xlab("Single or Mixed")

main.plot


grid.arrange(main.plot, accplot, bottom = "RT and Accuracy", 
             layout_matrix = matrix(c(1,2), ncol =2))


## 3.b DRT data ------------------------------------------------------------
myColors2 <- c("red","grey","green","purple")
names(myColors2) <- levels(drt.clean$condition)
colScale2 <- scale_colour_manual(name = "conditions",values = myColors2)

### RT Plots

rthist3 <- ggplot(data=drt.clean, aes(x=rt, colour = as.factor(condition))) +
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
  xlab("Condition") + coord_cartesian(ylim = c(.90,1))
drt.missplot

miss.acc <- glmer(response~condition+(1|subid),data=drt.missac,family=binomial("probit"))
Anova(miss.acc)
#   Response: response
#             Chisq Df Pr(>Chisq)    
# condition 344.39  3  < 2.2e-16 ***

round(100*tapply(drt.missac$response,drt.missac[,"condition"],mean),1)
#   1    2    3    4 
# 98.7 91.8 92.8 92.9 

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
# Response: response
#             Chisq Df Pr(>Chisq)    
# condition 151.04  3  < 2.2e-16 ***


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

clean.mean <- clean %>%
  group_by(subid,condition) %>%
  summarise(mean.rt = mean(rt1))

acc.mean <- cleanac %>%
  group_by(subid, condition) %>%
  summarise(perc_corr = mean(corr)*100)

vis.mean <- dplyr::left_join(acc.mean,clean.mean,by = c("subid","condition"))

drt.clean$condition <- recode_factor(drt.clean$condition, 
                                     "1" = "practice","2" = "green","3" ="white", "4" = "mixed",
                                     .default = "NA")
drt.ac$condition <- recode_factor(drt.ac$condition, 
                                  "1" = "practice","2" = "green","3" ="white", "4" = "mixed",
                                     .default = "NA")

drt.mean <- drt.clean %>%
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

all.clean <- all.combined %>%
  filter(condition != "practice", subid != 7)

all.clean$condition <- factor(all.clean$condition)

drt.model <- lm(cbind(drt.mean,mean.dev,drt.acc,perc_corr,mean.rt) ~ condition, data=all.clean)
Anova(drt.model)

#                drt.mean mean.dev   drt.acc  perc_corr   mean.rt
# (Intercept)    372.0016 3.902780  94.97846 99.5439189 1539.1465
# conditionmixed 223.7075 2.478629 -11.21800  0.1689189  467.8533
# conditionwhite 203.5509 1.795724 -13.76713 -0.1916307  901.5393


#Spencer's Model
base.mod <- lm(cbind(drt.mean,mean.dev,drt.acc,perc_corr,mean.rt) ~ condition, data=all.clean)
Anova(base.mod, test="Roy")

base.mod2 <- update(base.mod, . ~ .^2)
Anova(base.mod2, test="Roy")

summary(base.mod2)
#Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)          3.99481    0.18888  21.151  < 2e-16 ***
#   condition2           1.67905    0.22141   7.584 3.48e-14 ***
#   condition3           2.24291    0.22758   9.855  < 2e-16 ***
#   condition4           3.64677    0.22254  16.387  < 2e-16 ***
#   response             0.03071    0.18888   0.163  0.87083    
#   condition2:response -0.01566    0.22141  -0.071  0.94363    
#   condition3:response -0.10849    0.22758  -0.477  0.63357    
#   condition4:response -0.66033    0.22254  -2.967  0.00301 ** 


pairs(base.mod2)
heplot3d(base.mod2, wire = FALSE)


base.can <- candiscList(base.mod2)
# extract canonical R^2s
unlist(lapply(base.can, function(x) x$canrsq))

# condition1 condition2 
# 0.5902517  0.1639893 

op <- par(xpd=TRUE)
heplot(base.can, term="condition", scale=4, fill=TRUE, var.col="black", var.lwd=2)
par(op)


plot(base.can, term="condition")


# Fix Times ---------------------------------------------------------------
drt.combined.time <- drt.combined %>%
  mutate(subid = as.numeric(as.character(subid))) %>%
  mutate(UTC = as.POSIXct(UTC/1000, origin="1970-01-01", tz="UTC")) %>%
  mutate(trial_start =
           as.POSIXct(ifelse(subid != lag(subid),
                                         UTC - seconds(5 + rt/1000),
                                         UTC - seconds(rt/1000)),
                                  origin="1970-01-01", tz="UTC")) %>%
  mutate(stim_onset = UTC - seconds(rt/1000)) %>%
  select(trial_start,stim_onset, UTC, rt, subid, condition, response) %>%
  mutate(condition = recode_factor(condition, "1" = "practice", "2" = "green", "3" = "white", "4" = "mixed")) %>%
  filter(condition == "mixed")
  # mutate(windows = difftime(UTC, lag(UTC))) %>%
  # mutate(cycles = time - lag(time))


start_times <- drt.combined.time %>%
  group_by(subid) %>%
  summarise(starts = stim_onset[1])




vis.time <- vis.combined %>%
    filter(condition == "mixed") %>%
    arrange(subid) %>%
    group_by(subid) %>%
    mutate(tup = cumsum(rt1 + rt2)) %>%
    mutate(tlow = lag(tup, default = 0)) %>%
    mutate(block = cumsum(trial==1)) %>%
    select(subid, block, trial, time, targcol:rt2, tlow, tup) %>%
    ungroup(subid)


    # mutate(tlow = as.POSIXct("2016-10-28 15:18:37.063")) %>%
    # mutate(tup = as.POSIXct(tlow + seconds(tup/1000)))
    
    
  
  

# drt.time <- drt.clean %>%
#   filter(condition == "mixed") %>%
#   select(subid, rt, R, mean.dev) %>%
#   group_by(subid) %>%
#   mutate(pos = cumsum(rt + 5000)) %>%
#   ungroup(subid)

# drt.time$subid <- as.numeric(as.character(drt.time$subid))
# drt.time <- drt.time %>%
#   arrange(subid)
# drt.time$subid <- factor(drt.time$subid,levels = 1:38)






drt.time %>%
  group_by(subid) %>%
  summarise(n(), min(pos), max(pos))

vis.time %>%
  group_by(subid) %>%
  summarise(n(), min(tup), max(tup))




#group by subject id for both datasets
#If drt.time$pos == between vis.time$tlow & vis.time%tup, then bind vis.time to that row in drt.time

ggplot(all.test, aes(x = targcol, fill = targcol)) +
  scale_fill_manual(name = "targcol",values = c("forestgreen","grey")) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = 1,
              aes(y = vis.rt, col = "vis.rt")) + 
  stat_boxplot(geom ='errorbar', width = 0.1,position = position_dodge(.9),
               aes(y = vis.rt, col = "vis.rt")) +
  geom_violin(aes(y = drt, col = "drt")) +
  stat_boxplot(geom ='errorbar', width = 0.1,position = position_dodge(.9),
               aes(y = drt, col = "drt")) + theme_minimal()


# Crazy data.table technique ----------------------------------------------
  
bins <- unique(vis.time$subid)
tab1 <- data.table(drt.time)
setkey(tab1,"subid","pos")


tab2 <- data.table(vis.time)

tab2[,pos:=tlow]
setkey(tab2,"subid","pos")

x<-tab2[tab1, roll=TRUE, nomatch=0]

tab2[,pos:=tup]
setkey(tab2,"subid","pos")
y<-tab2[tab1, roll=-Inf, nomatch=0]

setkey(x,"subid","pos","tlow")
setkey(y,"subid","pos","tlow")
inBin<-x[y,nomatch=0]
inBin[, between:=TRUE]

setkey(tab1,"subid","pos")
setkey(inBin,"subid","pos")

result<-inBin[,list(subid,pos,between)][tab1]
result[is.na(between), between:=FALSE]

