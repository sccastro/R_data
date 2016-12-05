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
  dplyr::mutate(totalerrors = miss + falsealarm)

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
datac <- summarySEwithin(clean, measurevar ="totalerrors", withinvars = c("condition","r1g0"), idvar = "subid")
datac

datart <- summarySEwithin(clean.rt, measurevar ="rt", withinvars = c("condition","falsealarm"), idvar = "subid")
datart


# Make Plots --------------------------------------------------------------

acc.mean <- clean %>%
  group_by(condition, r1g0) %>%
  dplyr::summarise(perc_fa = mean(falsealarm)*100, perc_miss = mean(miss)*100, perc_error = mean(totalerrors)*100)

accplot1 <- ggplot(acc.mean, aes(x = as.factor(condition), y=perc_error, fill = as.factor(perc_miss))) +
  geom_bar(stat = "identity", position = "stack", aes(fill = as.factor(perc_miss))) + ylim(0,100) + theme_minimal() +  my.axis.font +
  xlab("Condition")
accplot1

accplot2 <- ggplot(datac, aes(x = condition, y=1-totalerrors, fill = r1g0)) +
  geom_bar(stat = "identity", position = "dodge", inherit.aes = TRUE, colour = "black") + 
  theme_minimal() +  my.axis.font + xlab("Condition") + coord_cartesian(ylim=c(.7,1)) +
  geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=1-totalerrors-ci, ymax=1-totalerrors+ci)) + 
  ggtitle("Percent \n Correct Responses in a \n Go-NoGo DRT Task \n")
accplot2

datac2 <- datac %>% 
  group_by(condition) %>%
  dplyr::summarise(totalerrors = mean(totalerrors),ci = mean(ci))
  

accplot3 <- ggplot(datac2, aes(x = as.factor(condition), y=1-totalerrors)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font + geom_errorbar(position=position_dodge(.9), width = .25, 
                aes(ymin=1-totalerrors-ci, ymax=1-totalerrors+ci)) +
  xlab("Condition") + coord_cartesian(ylim=c(.7,1))
accplot3


# Reaction Time Plots -----------------------------------------------------

rtplot1 <- ggplot(datart, aes(x = as.factor(falsealarm), 
                              fill = condition,y=rt_norm)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
  my.axis.font() + geom_errorbar(position=position_dodge(.9), width = .25, 
                               aes(ymin=rt_norm-ci, ymax=rt_norm+ci)) +
  scale_x_discrete("False Alarm", labels = c("incorrect","correct"))
rtplot1


# Simple Statistics -------------------------------------------------------

summaryBy(perc_error~condition,data=as.data.frame(clean.mean),FUN=function(x) {any(is.na(x))})

acc.lm <- with(as.data.frame(clean.mean), lmer(perc_error ~ condition * r1g0 + (1|subid)),na.omit()) 

summary(acc.lm)

# Linear mixed model fit by REML ['lmerMod']
# Formula: perc_error ~ condition * r1g0 + (1 | subid)
# 
# REML criterion at convergence: 1385.8
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.0552 -0.4722 -0.2102  0.0404  3.4156 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subid    (Intercept) 111.6    10.56   
# Residual             558.6    23.63   
# Number of obs: 153, groups:  subid, 26
# 
# Fixed effects:
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

with(as.data.frame(clean.mean), lmer(perc_error ~ condition  + (1|subid)),na.omit()) %>% anova()

# Analysis of Variance Table
#            Df Sum Sq Mean Sq F value
# condition  2 2.3095  1.1548   0.002






ggplot(data = combined.data, aes(x = UniqueID, y = Response.Time)) + 
  geom_point() +
  geom_smooth(method = "lm")


# Full Model --------------------------------------------------------------

full.lm <- with(as.data.frame(), lmer(perc_error ~ condition * r1g0 + (1|subid)),na.omit())





