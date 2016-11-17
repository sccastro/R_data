
rm(list = ls())
require(doBy)
require(dplyr)
setwd("Documents/R_Data/Nature_DRT")
files <- list.files(path="DRT_all", pattern="*.csv", full.names=T, recursive=FALSE)
subids <- regmatches(files, regexpr("\\d{3}", files))
condition <- as.character(lapply(strsplit(files, "_"), "[",4))

for (i in 1:length(files)) {
  infile <- read.csv(files[i])
  infile$subids <- as.factor(subids[i])
  infile$condition <- as.factor(condition[i])
  
  if(!exists("combined.data")) {
    combined.data <- infile
  }
  else {
    combined.data <- rbind(combined.data,infile)
  }
}


# Explore data ------------------------------------------------------------

summary(combined.data)
str(combined.data)


# Clean up dataset --------------------------------------------------------


clean.data <- combined.data %>%
  select(subid = subids, condition, isPractice,
         r1g0 = LED.Red..1.Red.0.Green., rt = Response.Time)

  
levels(clean.data$condition) <- list(PreTesting = ("PreTesting"), 
                            Testing = ("Testing"),
                            Posttesting = c("Posttesting","Postsetting","Posttetting"))
clean.data$isPractice <- factor(clean.data$isPractice)
levels(clean.data$isPractice) <- c("no","yes")
clean.data$r1g0 <- factor(clean.data$r1g0)
levels(clean.data$r1g0) <- c("g", "r")
clean.data$subid <- factor(clean.data$subid)

clean.data$isPractice[1:100]
clean.data$r1g0[1:100]


# Descriptives ------------------------------------------------------------


summaryBy(rt ~ subid + condition + r1g0, data = clean.data)

clean.data %>%
  group_by(subid, condition, r1g0) %>%
  summarise(total.count = n())


# Plot Reaction Times -----------------------------------------------------


