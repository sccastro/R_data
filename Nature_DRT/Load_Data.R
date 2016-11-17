# Get Data from Dropbox ---------------------------------------------------

require(dplyr)

# Get all files and apply functions to them -------------------------------

#Get filenames
files <- list.files(path="DRT_all", pattern="*.csv", full.names=T, recursive=FALSE)
#Get Subject IDs and conditions
subids <- regmatches(files, regexpr("\\d{3}", files))
condition <- c()
for (i in files) condition[i] <- regmatches(i, regexpr("[PT][roe][a-z]", i))

#vectorize condition and subids
condition <- unname(condition, force = FALSE)
condition <- as.vector(condition)

subids <- as.vector(subids)

newcols <- data.frame(condition,subids)



# Summary of the data types to merge --------------------------------------

#The 3 groups of data are condition, subids, and files. They are all characters with 78 elements. 
#I want to open the file associated with the filename in files, and append subids and condition to 
# every row in the file. Then I want to combine all of the files into one master file.


newcols[7:10,1]
str(newcols)

str(files)








#List of files with all data

#This creates a list of data.frames. datalist is a large matrix of 1092 elements
datalist <- sapply(files, function(x,...){
  read.csv(file=x,header=T)
  
  # I tried looping through the newcols data.frame and appending the element to a file. 
  #It didn't work 
  # for (i in newcols) datalist[[i]]$condition <- rep(newcols[i,1], nrow(datalist[[i]]$UniqueID))
  # for (j in newcols) onefile$subid <- rep(newcols[i,2], nrow(onfile))
  
  })







# merge all of the files together
for (i in 1:78) {
  bind_cols(datalist[[i]],data.frame(condition[[i]]))
}

datalist[[1]]$UniqueID

# combined <- bind_rows(datalist)
# summary(combined)



class(datalist[[1]])

datalist[[1]]
