rm(list = ls())
source("Load_Data.R")


combined.data <- CombineData()

clean.data <- CleanData(combined.data)


# Descriptives ------------------------------------------------------------
summaryBy(rt ~ subid + condition + r1g0, data = clean.data)


clean.data %>%
  group_by(subid, condition, r1g0) %>%
  dplyr::summarise(total.count = n())



# Calculate DV's -----------------------------------------------------
#Create accuracy
clean <- clean.data %>%
  filter(isPractice == "no", pressCount == 1 | pressCount == 0) %>%
  group_by(subid,condition,r1g0) %>%
  dplyr::mutate(falsealarm = ifelse(pressCount == 1 & r1g0 == "r", 1,0)) %>%
  dplyr::mutate(miss = ifelse(pressCount == 0 & r1g0 == "g", 1,0)) %>%
  dplyr::mutate(totalerrors = miss + falsealarm) %>%
  dplyr::mutate(correct = ifelse(totalerrors == 0, 1,0))

#summary of Signal, Noise, Hits and False Alarms
aprime <- clean %>%
  select(subid,condition,r1g0,pressCount,falsealarm)%>%
  group_by(subid,condition,r1g0) %>%
  dplyr::summarise(total.count = n())



clean.mean <- clean %>%
  group_by(subid,condition, r1g0) %>%
  dplyr::summarise(perc_error = mean(totalerrors)*100)

clean.rt <- clean %>%
  select(subid,condition,r1g0,rt,falsealarm) %>%
  filter(rt >= 150 & rt <= 2000)

# Outliers ----------------------------------------------------------------
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r") # Load the function


boxplot.with.outlier.label(clean.rt$rt~clean.rt$condition*clean.rt$r1g0,
                           label_name = clean.rt$subid)


notwanted <- c(2,3,5,8,14,15,20,22,23)

subjectsremoved <- clean %>%
  filter(!as.integer(subid) %in% notwanted)

levels(subjectsremoved$subid) <- droplevels(subjectsremoved$subid)


#rthist
medlines <- clean.rt %>%
  group_by(condition) %>%
  dplyr::summarise(rt.med = median(rt))

ggplot(data.frame(clean.rt), aes(rt, fill = condition)) + 
  geom_density(alpha=.8) + xlim(0,2000) + 
  geom_vline(data = as.data.frame(medlines), aes(xintercept=medlines$rt.med, color=condition),
             linetype="dashed", size=1)



# Get some Error Bars -----------------------------------------------------
datac <- summarySEwithin(clean, measurevar ="correct", withinvars = c("condition","r1g0"), idvar = "subid")
datac

<<<<<<< HEAD
datart <- summarySEwithin(clean.rt, measurevar ="rt", withinvars = c("condition","falsealarm"), idvar = "subid")
=======
datac1 <- summarySEwithin(clean, measurevar ="correct", withinvars = c("condition"), idvar = "subid")
datac1

datart <- summarySEwithin(clean, measurevar ="rt", withinvars = c("condition","correct"), idvar = "subid")
>>>>>>> origin/master
datart


# Make Plots --------------------------------------------------------------

<<<<<<< HEAD
acc.mean <- clean %>%
  group_by(condition, r1g0) %>%
  dplyr::summarise(perc_fa = mean(falsealarm)*100, perc_miss = mean(miss)*100, perc_error = mean(totalerrors)*100)

accplot1 <- ggplot(acc.mean, aes(x = as.factor(condition), y=perc_error, fill = as.factor(perc_miss))) +
  geom_bar(stat = "identity", position = "stack", aes(fill = as.factor(perc_miss))) + ylim(0,100) + theme_minimal() +  my.axis.font +
  xlab("Condition")
accplot1

accplot2 <- ggplot(datac, aes(x = condition, y=1-totalerrors, fill = r1g0)) +
=======
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
>>>>>>> origin/master
  geom_bar(stat = "identity", position = "dodge", inherit.aes = TRUE, colour = "black") + 
  theme_minimal() +  my.axis.font + xlab("Condition") + coord_cartesian(ylim=c(.7,1)) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=correct-ci, ymax=correct+ci)) + 
  ggtitle("Percent \n Correct Responses in a \n Go-NoGo DRT Task \n")
accplot3

datac2 <- datac %>% 
  group_by(condition) %>%
  dplyr::summarise(totalerrors = mean(totalerrors),ci = mean(ci))
  

<<<<<<< HEAD
accplot3 <- ggplot(datac2, aes(x = as.factor(condition), y=1-totalerrors)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=1-totalerrors-ci, ymax=1-totalerrors+ci)) +
  xlab("Condition") + coord_cartesian(ylim=c(.7,1))
accplot3
=======
#Reaction Time
rt.mean <- clean %>%
  group_by(condition, r1g0) %>%
  filter(r1g0 == "g" | falsealarm == 1) %>%
  dplyr::summarise(rt_mean = mean(rt))
>>>>>>> origin/master

#Count Rt's in group
clean %>%
  group_by(condition, r1g0) %>%
  filter(r1g0 == "g" | falsealarm == 1) %>%
  dplyr::summarise(total.count = n())

# Reaction Time Plots -----------------------------------------------------

rtplot1 <- ggplot(datart, aes(x = as.factor(falsealarm), 
                              fill = condition,y=rt_norm)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font() + geom_errorbar(position=position_dodge(.9), width = .25, 
                               aes(ymin=rt_norm-ci, ymax=rt_norm+ci)) +
  scale_x_discrete("False Alarm", labels = c("incorrect","correct"))
rtplot1


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
<<<<<<< HEAD
=======
require(lmerTest)

#Errors first
clean.mean <- clean %>%
  group_by(subid,condition, r1g0) %>%
  dplyr::summarise(perc_correct = mean(correct)*100) %>%
  filter(perc_correct > 50)
>>>>>>> origin/master

summaryBy(perc_correct~condition,data=as.data.frame(clean.mean),FUN=function(x) {any(is.na(x))})

<<<<<<< HEAD
acc.lm <- with(as.data.frame(clean.mean), lmer(perc_error ~ condition * r1g0 + (1|subid)),na.omit()) 

summary(acc.lm)

# Linear mixed model fit by REML ['lmerMod']
# Formula: perc_error ~ condition * r1g0 + (1 | subid)
=======
acc.lm <- with(as.data.frame(clean.mean), lmer(perc_correct ~ condition  + (1|subid)),na.omit()) 

summary(acc.lm)

# Linear mixed model fit by REML 
# t-tests use  Satterthwaite approximations to degrees of freedom ['lmerMod']
# Formula: perc_correct ~ condition + (1 | subid)
>>>>>>> origin/master
# 
# REML criterion at convergence: 1385.8
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
<<<<<<< HEAD
# -2.0552 -0.4722 -0.2102  0.0404  3.4156 
=======
# -3.5690  0.1834  0.2061  0.3819  1.5859 
>>>>>>> origin/master
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subid    (Intercept) 111.6    10.56   
# Residual             558.6    23.63   
# Number of obs: 153, groups:  subid, 26
# 
# Fixed effects:
<<<<<<< HEAD
#                             Estimate Std. Error t value
# (Intercept)                  10.144      5.170   1.962
# conditionTesting              5.171      6.699   0.772
# conditionPosttesting          8.269      6.627   1.248
# r1g0r                        -1.502      6.685  -0.225
# conditionTesting:r1g0r       -9.593      9.413  -1.019
# conditionPosttesting:r1g0r  -16.053      9.362  -1.715
# 
# Correlation of Fixed Effects:
#   (Intr) cndtnT cndtnP r1g0r  cnT:10
# condtnTstng -0.648                            
# cndtnPsttst -0.655  0.506                     
# r1g0r       -0.646  0.499  0.504              
# cndtnTst:10  0.459 -0.710 -0.358 -0.710       
# cndtnPst:10  0.462 -0.356 -0.706 -0.714  0.507
=======
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
>>>>>>> origin/master

with(as.data.frame(clean.mean), lmer(perc_error ~ condition  + (1|subid)),na.omit()) %>% anova()

# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
#            Sum Sq Mean Sq NumDF  DenDF   F.value Pr(>F)
# condition 2.3095  1.1548     2 126.14 0.0019587  0.998


#Reaction Time

summaryBy(rt~condition + totalerrors,data=as.data.frame(clean),FUN=function(x) {any(is.na(x))})

rt.lm <- with(as.data.frame(clean), lmer(rt ~ condition * totalerrors  + (1|subid)),na.omit()) 

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


# Full Model --------------------------------------------------------------

full.lm <- with(as.data.frame(), lmer(perc_error ~ condition * r1g0 + (1|subid)),na.omit())





