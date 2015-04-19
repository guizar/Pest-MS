# Install missing packages before running the sequence of scripts

list <- c("gdata","sp","FNN","countrycode","reshape2","plyr","dplyr","classInt","ggplot2","gridExtra","R.matlab","arrayhelpers","raster")

new.packages <- list[!(list %in% installed.packages()[,"Package"])]

print(paste("You are missing the following packages",new.packages,sep = ": "))
print("Please continue evaluating this code to proceed with the installation")

if(length(new.packages)) install.packages(new.packages)
rm(list)
rm(new.packages)
