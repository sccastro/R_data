# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_Data.R")
require(tidyverse)
require(gridExtra)
require(lubridate)
require(car)


# Data Cleaning -----------------------------------------------------------


fulldata <- read.csv("Turrill_DRT2_23.csv")
names(fulldata) <- c("sub", "difficulty", "trial","date","file_write","DRT","rt",
                     "code","button_down","button_up", "hit_light","hit_ISO","late","miss","filter","baseline")

flist <- c("sub", "difficulty", "trial","code","hit_light","hit_ISO","late","miss","filter","baseline")
sublist <- unique(fulldata$sub)
trialist <-  unique(fulldata$trial)
codelist <-  unique(fulldata$code)
hitlist  <- unique(fulldata$hit_ISO) 
difficulty <-  unique(fulldata$difficulty)

fulldata <- read_csv("Turrill_DRT2_23.csv", skip = 1,
                     col_names = c("sub", "difficulty", "trial","date","file_write","DRT","rt",
                                   "code","button_down","button_up", "hit_light","hit_ISO","late",
                                   "miss","filter","baseline"),
                     col_types = cols(sub = col_factor(levels = sublist), 
                                      difficulty = col_factor(levels = difficulty),
                                      trial = col_factor(levels = trialist, ordered = TRUE), 
                                      date = col_date(format = "%m/%d/%y"),
                                      file_write = col_character(),
                                      DRT = col_double(),
                                      rt = col_double(),
                                      code = col_factor(levels = codelist),
                                      button_down = col_character(),
                                      button_up = col_character(),
                                      hit_light = col_number(),
                                      hit_ISO = col_number(),
                                      late = col_factor(levels = c(0,1)),
                                      miss = col_factor(levels = c(0,1)),
                                      filter = col_factor(levels = c(0,1)),
                                      baseline = col_double()))


fulldata$difficulty <- factor(fulldata$difficulty, labels = c("easy1","easy2","hard1","hard2"))


# Data Exploration --------------------------------------------------------

(subs <- fulldata %>%
   select(-date, -file_write,-button_down, -button_up, -DRT) %>%
   filter(filter == 1, late == 0) %>%
   group_by(difficulty, miss, code, baseline) %>%
   summarise_each(funs(Count = sum(.), Percent = sum(.)/n()*100), 
                  -sub, -trial, -difficulty, -miss, -late, -code,-filter,-rt))


