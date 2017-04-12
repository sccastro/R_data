# Clear/setup workspace ---------------------------------------------------
rm(list=ls()); par(mfrow = c(1,1))
source("Load_Data.R")
require(gridExtra)
require(lubridate)
require(car)
library(ggrepel)

# 1. Load Data -------------------------------------------------------
rm(list = ls())
source("Load_Data.R")
load("CleanData.Rdata")

# Fix DRT Times ---------------------------------------------------------------
drt.combined.time <- drt.combined %>%
  mutate(subid = as.numeric(as.character(subid), rt = as.numeric(as.character(rt)))) %>%
  mutate(UTC = as.POSIXct(UTC/1000, origin="1970-01-01", tz="UTC")) %>%
  mutate(trial_start =
           as.POSIXct(ifelse(rt != -1, 
                             ifelse(subid != lag(subid, default = 0),
                                    UTC - seconds(5 + rt/1000),
                                    UTC - seconds(rt/1000)),
                             UTC),
                      origin="1970-01-01", tz="UTC")) %>%
  group_by(subid, condition) %>%
  mutate(block = cumsum(
    ifelse(
      as.numeric(difftime(trial_start, 
               lag(
        trial_start, default = 4), 
        units = "secs")) > 50, 1,0))) %>%
  ungroup(subid, condition) %>%
  mutate(stim_onset = UTC - seconds(rt/1000)) %>%
  select(subid,block, trial_start,stim_onset, UTC, rt, condition, response) %>%
  mutate(condition = recode_factor(condition, "1" = "practice", "2" = "green", "3" = "white", "4" = "mixed")) %>%
  filter(condition == "mixed")


drt.combined.time <- as_tibble(drt.combined.time)

drt.combined.time %>%
  group_by(subid) %>%
  dplyr::summarise(blocks = length(unique(block))) %>%
  filter(blocks != 4)



# Fix Visual Times --------------------------------------------------------

vis.time <- vis.combined %>%
  filter(condition == "mixed") %>%
  arrange(subid) %>%
  group_by(subid) %>%
  mutate(tup = cumsum(rt1 + rt2)) %>%
  mutate(tlow = lag(tup, default = 0)) %>%
  mutate(block = cumsum(trial==1)) %>%
  select(subid, block, trial, time, targcol:rt2, tlow, tup) %>%
  ungroup(subid)

vis.time %>%
  group_by(subid) %>%
       dplyr::summarise(blocks = length(unique(block))) %>%
  filter(blocks != 4)
  