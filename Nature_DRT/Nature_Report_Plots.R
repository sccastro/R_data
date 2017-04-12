rm(list = ls())
source("Load_Data.R")


combined.data <- CombineData()

clean.data <- CleanData(combined.data)


# Descriptives ------------------------------------------------------------
str(clean.data)

summaryBy(rt ~ subid + condition + r1g0, data = clean.data)


clean.data %>%
  group_by(subid, condition, r1g0) %>%
  dplyr::summarise(total.count = n())



# Calculate DV's -----------------------------------------------------
#Create accuracy
clean <- clean.data %>%
  filter(isPractice == "no") %>%
  group_by(subid,condition,r1g0) %>%
  dplyr::mutate(falsealarm = ifelse(pressCount == 1 & r1g0 == "r", 1,0)) %>%
  dplyr::mutate(miss = ifelse(pressCount == 0 & r1g0 == "g", 1,0)) %>%
  dplyr::mutate(totalerrors = miss + falsealarm) %>%
  dplyr::mutate(correct = ifelse(totalerrors == 0, 1,0))

cleanrt <- clean %>%
  filter(rt > 100 & rt < 2000)

#summary of Signal, Noise, Hits and False Alarms
aprime <- clean %>%
  select(subid,condition,r1g0,correct,falsealarm)%>%
  group_by(subid,condition,r1g0) %>%
  dplyr::summarise(hitrate = mean(correct), farate = mean(falsealarm))

hitrate1 <- aprime[c(TRUE,FALSE),]$hitrate
farate1 <- aprime[aprime$r1g0 == "r",]
farate2 <- farate1$farate 

aprime1 <- aprime %>%
  group_by(subid,condition) %>%
  summarise(hitrate = mean(hitrate),farate = mean(farate))

aprime1$hitrate <- hitrate1
aprime1$farate <- farate2

aprime2 <- Getaprime(aprime1)
aprimeclean <- aprime2 %>%
  filter(aprimescore > .5)

# Get some Error Bars -----------------------------------------------------
datac <- summarySEwithin(clean, measurevar ="correct", withinvars = c("condition","r1g0"), idvar = "subid")
datac

datac1 <- summarySEwithin(clean, measurevar ="correct", withinvars = c("condition"), idvar = "subid")
datac1

dataca <- summarySE(aprimeclean, measurevar = "aprimescore", groupvars = "condition", na.rm = TRUE)
dataca

datart <- summarySEwithin(clean, measurevar ="rt", withinvars = c("condition","correct"), idvar = "subid")
datart


# Make Plots --------------------------------------------------------------

# #Accuracy
# acc.mean <- clean %>%
#   group_by(condition, r1g0) %>%
#   dplyr::summarise(perc_fa = mean(falsealarm)*100, perc_miss = mean(miss)*100, perc_error = mean(totalerrors)*100)

#Count Errors
# real_errors <-clean %>%
#   group_by(condition, totalerrors = as.factor(totalerrors)) %>%
#   dplyr::summarise(total.count = n()) %>%
#   mutate(freq = total.count/sum(total.count)) %>%
#   filter(totalerrors == "0")

# errorplot <- ggplot(real_errors, aes(x = as.factor(condition), y=freq, fill = totalerrors)) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = totalerrors)) +
#   ylim(0,1) + theme_minimal() +  my.axis.font +
#   xlab("Condition")
# errorplot

# accplot1 <- ggplot(acc.mean, aes(x = as.factor(condition), y=perc_error, fill = as.factor(perc_miss))) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = as.factor(perc_miss))) + ylim(0,100) + theme_minimal() +  my.axis.font +
#   xlab("Condition")
# accplot1

accplot2 <- ggplot(datac, aes(x = condition, y=correct, fill = r1g0)) +
  geom_bar(stat = "identity", position = "dodge", inherit.aes = TRUE, colour = "black") + 
  theme_minimal() +  my.axis.font + xlab("Condition") + coord_cartesian(ylim=c(.7,1)) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=correct-ci, ymax=correct+ci)) + 
  ggtitle("Percent \n Correct Responses in a \n Go-NoGo DRT Task \n") + 
  scale_fill_brewer(type = "qual", palette = 2,direction = 1)
accplot2

accplot3 <- ggplot(datac1, aes(x = condition, y=correct)) +
  geom_bar(stat = "identity", position = "dodge", inherit.aes = TRUE, colour = "black") + 
  theme_minimal() +  my.axis.font + xlab("Condition") + coord_cartesian(ylim=c(.7,1)) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=correct-ci, ymax=correct+ci)) + 
  ggtitle("Percent \n Correct Responses in a \n Go-NoGo DRT Task \n")
accplot3

#Aprime Score Plots
# accplot4 <- ggplot(aprime2, aes(x = condition, y=aprimescore, color = subid)) + 
#   geom_point(position = "jitter") + theme_minimal() + ylim(.8,1.1) + 
#   ggtitle("A\' \n by Condition \n and Subject") + my.axis.font
# accplot4

accplot5 <- ggplot(dataca, aes(x = condition, y=aprimescore)) + 
  geom_bar(stat = "identity") + theme_minimal() +
  my.axis.font + ggtitle("A\' \n by Condition") +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=aprimescore-ci, ymax=aprimescore+ci))
accplot5



#Reaction Time
rt.mean <- clean %>%
  group_by(condition, r1g0) %>%
  filter(r1g0 == "g" | falsealarm == 1) %>%
  dplyr::summarise(rt_mean = mean(rt))

#Count Rt's in group
clean %>%
  group_by(condition, r1g0) %>%
  filter(r1g0 == "g" | falsealarm == 1) %>%
  dplyr::summarise(total.count = n())


# rtplot1 <- ggplot(rt.mean, aes(x = as.factor(condition), y=rt_mean, fill = as.factor(r1g0))) +
#   geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(r1g0))) + ylim(0,1200) + theme_minimal() +  my.axis.font +
#   xlab("Condition") + scale_fill_brewer(type = "qual", palette = 2,direction = 1)
# rtplot1

rtplot2 <- ggplot(datart, aes(x = as.factor(condition), y=rt, fill = as.factor(correct))) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(correct))) + 
  ylim(0,510) + theme_minimal() +  my.axis.font +
  xlab("Condition") + scale_fill_brewer(type = "qual", palette = 2,direction = 1) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=rt-ci, ymax=rt+ci)) + 
  ggtitle("Reaction Time \n by Condition in a \n Go-NoGo DRT Task \n")
rtplot2

rtplot3 <- ggplot(datart, aes(x = as.factor(correct), y=rt, fill = as.factor(condition))) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(condition))) + 
  ylim(0,510) + theme_minimal() +  my.axis.font +
  xlab("Correct") + scale_fill_brewer(type = "qual", palette = 2,direction = 1) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=rt-ci, ymax=rt+ci)) + 
  ggtitle("Reaction Time \n by Condition in a \n Go-NoGo DRT Task \n")
rtplot3

# Simple Statistics -------------------------------------------------------
require(lmerTest)

#Errors first
clean.mean <- clean %>%
  group_by(subid,condition, r1g0) %>%
  dplyr::summarise(perc_correct = mean(correct)*100) %>%
  filter(perc_correct > 50)

summaryBy(perc_correct~condition,data=as.data.frame(clean.mean),FUN=function(x) {any(is.na(x))})

acc.lm <- with(as.data.frame(clean.mean), lmer(perc_correct ~ condition + (1|subid)),na.omit()) 

summary(acc.lm)

# Linear mixed model fit by REML 
# t-tests use  Satterthwaite approximations to degrees of freedom ['lmerMod']
# Formula: perc_correct ~ condition + (1 | subid)
# 
# REML criterion at convergence: 1412.5
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.5690  0.1834  0.2061  0.3819  1.5859 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subid    (Intercept) 107.1    10.35   
# Residual             589.6    24.28   
# Number of obs: 153, groups:  subid, 26
# 
# Fixed effects:
#   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)           90.6184     4.0006  82.4800  22.651   <2e-16 ***
#   conditionTesting      -0.2709     4.8447 126.6200  -0.056    0.955    
# conditionPosttesting  -0.2539     4.8191 126.3600  -0.053    0.958    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#   (Intr) cndtnT
# condtnTstng -0.613       
# cndtnPsttst -0.617  0.509

with(as.data.frame(clean.mean), lmer(perc_error ~ condition  + (1|subid)),na.omit()) %>% anova()

# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
#            Sum Sq Mean Sq NumDF  DenDF   F.value Pr(>F)
# condition 2.3095  1.1548     2 126.14 0.0019587  0.998


#Reaction Time
summaryBy(rt~condition + correct ,data=as.data.frame(cleanrt))
summaryBy(rt~condition + totalerrors,data=as.data.frame(clean),FUN=function(x) {any(is.na(x))})

rt.lm <- with(as.data.frame(cleanrt), lmer(rt ~ condition * correct  + (1|subid)),na.omit()) 

summary(rt.lm)

# Linear mixed model fit by REML 
# t-tests use  Satterthwaite approximations to degrees of freedom ['lmerMod']
# Formula: rt ~ condition * totalerrors + (1 | subid)

# REML criterion at convergence: 226314.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.5646 -0.3820 -0.2787 -0.1864 17.4555 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subid    (Intercept)  2754     52.47  
# Residual             97326    311.97  
# Number of obs: 15798, groups:  subid, 26
# 
# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                        120.836     11.327    33.000  10.668 3.38e-12 ***
#   conditionTesting                    -6.719      6.579 15563.000  -1.021    0.307    
# conditionPosttesting               -27.076      6.403 15569.000  -4.228 2.37e-05 ***
#   totalerrors                        -14.350     18.521 10430.000  -0.775    0.438    
# conditionTesting:totalerrors        14.840     24.162 15381.000   0.614    0.539    
# conditionPosttesting:totalerrors   282.702     27.359 13505.000  10.333  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#   (Intr) cndtnT cndtnP ttlrrr cndtT:
#   condtnTstng -0.298                            
# cndtnPsttst -0.311  0.530                     
# totalerrors -0.146  0.220  0.259              
# cndtnTstng:  0.093 -0.298 -0.167 -0.620       
# cndtnPstts:  0.096 -0.135 -0.271 -0.691  0.417

with(as.data.frame(clean), lmer(rt ~ condition * totalerrors  + (1|subid)),na.omit()) %>% anova()

# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
# Sum Sq Mean Sq NumDF DenDF F.value    Pr(>F)    
#   condition              1942171  971086     2 15780   9.978 4.672e-05 ***
#   totalerrors            5228641 5228641     1 12032  53.723 2.456e-13 ***
#   condition:totalerrors 11995654 5997827     2 14617  61.626 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# D' ----------------------------------------------------------------------


#Pre D's
dprimepre <- clean %>%
  filter(condition == "PreTesting") %>%
  dplyr::mutate(resp = ifelse(pressCount == 1, "Y","N")) %>%
  dplyr::rename(acc = correct) %>%
  select(subject=subid,condition,r1g0,resp,acc)

dprimepre$condition <- droplevels(dprimepre$condition)

pre <- data.frame(dprime(dprimepre)) %>%
  mutate(condition = "pre")

#Test D's
dprimetest <- clean %>%
  filter(condition == "Testing") %>%
  dplyr::mutate(resp = ifelse(pressCount == 1, "Y","N")) %>%
  dplyr::rename(acc = correct) %>%
  select(subject=subid,condition,r1g0,resp,acc)

dprimetest$condition <- droplevels(dprimetest$condition)

test <- data.frame(dprime(dprimetest)) %>%
  mutate(condition = "test")

#Post D's

dprimepost <- clean %>%
  filter(condition == "Posttesting") %>%
  dplyr::mutate(resp = ifelse(pressCount == 1, "Y","N")) %>%
  dplyr::rename(acc = correct) %>%
  select(subject=subid,condition,r1g0,resp,acc)

dprimepost$condition <- droplevels(dprimepost$condition)

post <- data.frame(dprime(dprimepost)) %>%
  mutate(condition = "post")

#gather them
dprimes <- bind_rows(pre, bind_rows(test,post)) %>%
  arrange(subject)

dprimes$condition <- factor(dprimes$condition, levels = c("pre","test","post"))


#plot them
dprimeplot <- ggplot(dprimes, aes(x=condition, y = Freq, color = as.factor(condition))) + 
  geom_boxplot(na.rm = TRUE) + ylim(-1,5)
dprimeplot + theme_minimal() + my.axis.font

#D prime seems pretty useless to me.

