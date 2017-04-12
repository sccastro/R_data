# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_DataSAT.R")
require(tidyverse)
require(gridExtra)
require(lubridate)
require(car)
require(lme4)
require(stringr)
library(forcats)


# Get Data ----------------------------------------------------------------


## 1.b Get all DRT files and apply functions to them
if(!exists("drt.combined")) {
  drt.combined <- ExtractDRT('4_3_2017') #Only run this once or it will append to drt.data
}


##1.c Get all Steering data files and apply functions to them
if(!exists("steering.combined")) {
  steering.combined <- ExtractSteering('4_3_2017') #Only run this once or it will append to drt.data
}


### Get steering deviation

steering.combined$cursorpos <- as.numeric(as.character(steering.combined$cursorpos))

steering.combined$condition <- recode_factor(steering.combined$condition, "2" = "DRT", "1" = "Equal", "3" = "Steering")

steer.dev <- steering.combined %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos))
  # dplyr::group_by(subid, condition) %>%
  # dplyr::summarise(mean.dev = mean(deviation))

mean.steer.dev <- steering.combined %>%
  na.omit() %>%
  mutate(deviation = abs(ballpos - cursorpos)) %>%
  dplyr::group_by(subid,condition) %>%
  dplyr::summarise(mean.dev = mean(deviation))


tapply(steer.dev$deviation,as.factor(steer.dev$condition),mean)
### combine steering with drt

datdev <- summarySEwithin(steer.dev, measurevar ="deviation", withinvars = c("condition"), idvar = "subid")

RMSplot <- ggplot(datdev, aes(x = condition, y=deviation, group = 1)) +
  geom_point(stat = "identity") + theme_minimal() + geom_line() +
  my.axis.font + geom_errorbar(width = .25,
                               aes(ymin=deviation-ci, ymax=deviation+ci)) +
  xlab("Condition") + coord_cartesian(ylim = c(0,7.5)) +
  ylab("RMS Steering Error") + 
  scale_y_continuous(position = "right")




# DRT ---------------------------------------------------------------------

##Combine steering and RT
drt.combined$condition <- recode_factor(drt.combined$condition, "2" = "DRT", "1" = "Equal", "3" = "Steering")
drt.combined <- dplyr::left_join(as_tibble(drt.combined), mean.steer.dev, by = c("subid","condition"))

drt.combined$response <- recode_factor(drt.combined$response, miss = "0", hit = "1") #refactor for getting percentages

drt.combined$response <-as.numeric(as.character(drt.combined$response))
drt.combined <- drt.combined[complete.cases(drt.combined),]


##Reaction Time

drt.clean <- drt.combined %>%
  dplyr::select(subid,time,condition,rt,s1,R,response) %>%
  na.omit() %>%
  arrange(subid) %>%
  filter(rt > 150 & rt < 3000)

tapply(drt.clean$rt,as.factor(drt.clean$condition),mean)

datrt <- summarySEwithin(drt.clean, measurevar ="rt", withinvars = c("condition"), idvar = "subid")

rtplot <- ggplot(datrt, aes(x = condition, y=rt, group = 1)) +
  geom_point(stat = "identity") + theme_minimal() + geom_line() +
  my.axis.font + geom_errorbar(width = .25,
                               aes(ymin=rt-ci, ymax=rt+ci)) +
  xlab("Condition") + coord_cartesian(ylim = c(200,400)) +
  scale_fill_manual(values=myColors) + ylab("Reaction Time (ms)")


##Accuracy
drt.ac <- drt.combined %>%
  dplyr::select(condition,rt,s1,R,response) %>%
  na.omit() %>%
  filter(rt > 150 | rt < 1, rt < 3000)

tapply(as.numeric(as.character(drt.ac$response)), as.factor(drt.ac$condition), mean)



# summary -----------------------------------------------------------------

grid.arrange(rtplot,RMSplot, bottom = "Summary", 
             top = "Reaction Time and RMS Steering Error \n by Condition", 
             layout_matrix = matrix(c(1,2), ncol =2))




# By Subject --------------------------------------------------------------

RMSplotsub <- ggplot(steer.dev, aes(x = condition, y=deviation, color = subid, group = subid)) +
  stat_summary(fun.y="mean", geom="point") + 
  stat_summary(fun.y="mean", geom="line") + theme_minimal() +
  my.axis.font  +
  xlab("Condition") + coord_cartesian(ylim = c(0,7.5)) +
  ylab("RMS Steering Error") + 
  scale_y_continuous(position = "right")


rtplotsub <- ggplot(drt.clean, aes(x = condition, y=rt, group = subid, color = subid)) +
  stat_summary(fun.y="mean", geom="point") + 
  stat_summary(fun.y="mean", geom="line") + theme_minimal() +
  xlab("Condition") + coord_cartesian(ylim = c(200,400)) + ylab("Reaction Time (ms)") +
  my.axis.font +  guides(color=FALSE)


grid.arrange(rtplotsub,RMSplotsub, bottom = "Summary", 
             top = "Reaction Time and RMS Steering Error \n by Condition", 
             layout_matrix = matrix(c(1,2), ncol =2))


# Hit Rate by Reaction Time -----------------------------------------------
rt_hit <- drt.combined %>%
  na.omit() %>%
  dplyr::filter(rt > 150 | rt < 1, rt < 3000) %>%
  summarise(mean.rt = ifelse(rt == -1, ))


rtbyhit <- ggplot(drt.clean, aes(x = mean.hit, y=rt, group = condition, color = condition)) +
  stat_summary(fun.y="mean", geom="point") + 
  stat_summary(fun.y="mean", geom="line") + theme_minimal() +
  xlab("Condition") + coord_cartesian(ylim = c(200,400)) + ylab("Reaction Time (ms)") +
  my.axis.font +  guides(color=FALSE)
