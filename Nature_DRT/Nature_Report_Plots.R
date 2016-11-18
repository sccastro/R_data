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
  dplyr::mutate(totalerrors = miss + falsealarm)


# Get some Error Bars -----------------------------------------------------
datac <- summarySEwithin(clean, measurevar ="totalerrors", withinvars = c("condition","r1g0"), idvar = "subid")
datac

datart <- summarySEwithin(clean, measurevar ="rt", withinvars = c("condition","totalerrors"), idvar = "subid")
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





# Simple Statistics -------------------------------------------------------
clean.mean <- clean %>%
  group_by(subid,condition, r1g0) %>%
  dplyr::summarise(perc_error = mean(totalerrors)*100)

summaryBy(perc_error~condition,data=as.data.frame(clean.mean),FUN=function(x) {any(is.na(x))})

acc.lm <- with(as.data.frame(clean.mean), lmer(perc_error ~ condition  + (1|subid)),na.omit()) 

summary(acc.lm)

# Linear mixed model fit by REML ['lmerMod']
# Formula: perc_error ~ condition + (1 | subid)
# 
# REML criterion at convergence: 1412.5
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.5859 -0.3819 -0.2061 -0.1834  3.5690 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subid    (Intercept) 107.1    10.35   
# Residual             589.6    24.28   
# Number of obs: 153, groups:  subid, 26
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            9.3816     4.0006   2.345
# conditionTesting       0.2709     4.8447   0.056
# conditionPosttesting   0.2539     4.8191   0.053
# 
# Correlation of Fixed Effects:
#             (Intr) cndtnT
# condtnTstng -0.613       
# cndtnPsttst -0.617  0.509

with(as.data.frame(clean.mean), lmer(perc_error ~ condition  + (1|subid)),na.omit()) %>% anova()

# Analysis of Variance Table
#            Df Sum Sq Mean Sq F value
# condition  2 2.3095  1.1548   0.002






ggplot(data = combined.data, aes(x = UniqueID, y = Response.Time)) + 
  geom_point() +
  geom_smooth(method = "lm")

