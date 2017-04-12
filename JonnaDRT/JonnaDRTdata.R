rm(list = ls()); par(mfrow = c(1,1))
source("Load_Data.R")

require(gridExtra)

# Extract and Clean Data --------------------------------------------------------------

if(!exists("fulldata")) {
fulldata <- ExtractJDRT("data")
}


l1 <- c("time", "rt")
l2 <- c("date")

newdata <- fulldata %>%
  mutate_each_(funs(as.character), l1) %>%
  mutate_each_(funs(as.numeric), l1) %>%
  mutate_each_(funs(as.Date.factor), l2) %>%
  mutate(rt = round(rt/1000,0)) %>%
  mutate(late = ifelse(2500 > rt & rt > 1000, 0,1), miss = ifelse(rt >= 2500,0,1)) %>%
  mutate(condition = ifelse(condition == "easy1" | condition == "easy2", "easy","hard")) %>%
  na.omit() %>%
  filter(code != "X",  code != "V") %>%
  mutate(code = droplevels(code))

ggplot(as.data.frame(newdata), aes(x = rt, linetype = condition, color = code)) +
  stat_ecdf() + coord_cartesian(xlim = c(0,1500)) + 
  theme_dark() + scale_color_brewer(type = "div",palette = 9)

# Data summary ------------------------------------------------------------
rtdata <- as.data.frame(newdata) %>%
  filter(miss == 1) %>%
  na.omit() %>%
  group_by(code,condition) %>%
  dplyr::summarise(mean_rt = mean(rt), med_rt = median(rt), sd_rt = sd(rt))
  


# Summary Dataframes ------------------------------------------------------
rtdata <- as.data.frame(newdata) %>%
  filter(miss == 1) %>%
  na.omit() 
# dplyr::summarise(mean_rt = mean(rt))

acc <- as.data.frame(newdata) %>%
  mutate(acc = miss*100) 
# ggplot(aes(x = condition, y = mean_rt)) + geom_bar(stat = "identity")

allrt <- summarySEwithin(rtdata, measurevar = "rt", withinvars = c("code","condition"), idvar = "subid")

allac <- summarySEwithin(acc, measurevar = "miss", withinvars = c("code","condition"), idvar = "subid")

# Save Data Frames --------------------------------------------------------


# Plots -------------------------------------------------------------------


rtplot <- ggplot(allrt, aes(x = as.factor(condition), y=rt, shape = code, group=code)) + 
  geom_line() + ylim(350,650) + geom_point(size=3) +
  geom_errorbar(width = .25, aes(ymin=rt-se, ymax=rt+se)) +
  xlab("") + ylab("Reaction Time (ms)") +
  theme_classic() + guides(shape=FALSE) + my.axis.font

accplot <- ggplot(allac, aes(x = as.factor(condition), y=miss, shape = code, group=code)) +
  geom_line() + geom_point(size=3)  + 
  scale_shape_discrete(labels = c("Off Task", "Ones", "Threes", "Baseline"),
                     name = "Load Type") + xlab("") +
  geom_errorbar(width = .25, aes(ymin=miss-se, ymax=miss+se)) +
  theme_classic() + scale_y_continuous(position = "right", name = "Accuracy (0-1)")  +
  coord_cartesian(ylim = c(.85,1)) + my.axis.font + theme(legend.position = "bottom")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(accplot)

lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3))


grid.arrange(rtplot, accplot + guides(shape = FALSE),legend,
             top = "Residual Cognitive Load",bottom = "Driving Difficulty",
             layout_matrix = lay)

