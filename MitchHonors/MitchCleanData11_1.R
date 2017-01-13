# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_Data.R")
require(tidyr)

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
vis.combined$condition <- factor(vis.combined$condition, levels = c("practice","white","green","mixed"))

###Select the columns we want, clean out NA's
clean <- vis.combined %>%
  select(subid,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt1 > 150 & rt1 < 6000, rt2 > 150 & rt2 < 3000)

hist(clean$rt1, breaks = 100)

cleanac <- vis.combined %>%
  select(subid,condition,targcol,corr,rt1,rt2) %>%
  na.omit() %>%
  arrange(subid)


##2.b Clean drt.combined data 

### Get steering deviation

steering.combined$cursorpos <- as.numeric(as.character(steering.combined$cursorpos))

steer.dev <- steering.combined %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos)) %>%
  group_by(subid, condition) %>%
  summarise(mean.dev = mean(deviation), n=n())

### combine steering with drt

drt.combined <- left_join(drt.combined, steer.dev, by = c("subid","condition"))

drt.combined$response <- recode_factor(drt.combined$response, 
                                       wrong ="-1", correct ="1", miss = "0",
                                       .default = "NA") #refactor for getting percentages

drt.combined$response <-as.numeric(as.character(drt.combined$response))
drt.combined <- drt.combined[complete.cases(drt.combined),]


drt.clean <- drt.combined %>%
  select(subid = subid,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 & rt < 3000)

drt.ac <- drt.combined %>%
  select(subid = subid,condition,rt,s1,R,response,mean.dev) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 | rt < 1, rt < 3000)



# 3. Descriptive Stats -------------------------------------------------------
## 3.a Visual Plots first

###Colors for plots
####Create a custom color scale
myColors <- c("red","grey","green","purple")
names(myColors) <- levels(clean$condition)
colScale <- scale_colour_manual(name = "conditions",values = myColors)


### Visual Search Plots

#### Accuracy
#### Percentage correct by group
acc.mean <- cleanac %>%
  group_by(condition) %>%
  summarise(perc_corr = mean(corr)*100)


datac <- summarySEwithin(cleanac, measurevar ="corr", withinvars = c("condition"), idvar = "subid")

accplot <- ggplot(datac, aes(x = as.factor(condition), y=corr)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25,
                               aes(ymin=corr-ci, ymax=corr+ci)) +
  xlab("Condition")
accplot

#### Reaction Time
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

##### Violin plots of Reaction Time

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
  xlab("Condition")
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


# 4. Hypothesis and Error Plots - Full Model ----------------------------------------------
library(heplots)
library(candisc)

clean.mean <- clean %>%
  group_by(subid,condition) %>%
  summarise(mean.rt = mean(rt1))

acc.mean <- cleanac %>%
  group_by(subid, condition) %>%
  summarise(perc_corr = mean(corr)*100)

vis.mean <- left_join(acc.mean,clean.mean,by = c("subid","condition"))

drt.clean$condition <- recode_factor(drt.clean$condition, 
                                     "1" = "green","2" = "white","3" ="mixed", "4" = "practice",
                                     .default = "NA")
drt.ac$condition <- recode_factor(drt.ac$condition, 
                                     "1" = "green","2" = "white","3" ="mixed", "4" = "practice",
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
