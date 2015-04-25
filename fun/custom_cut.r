# This function generates breaks using the extreme 5th - 95th percentile range

# by: indicate the NUMBER that will be used to make incremental steps from stats[1]:stats[5]. The size of the step affect the total number of breaks

# digits:  is used to provide an argument for the round function

# works with library(classInt)

custom.cut <- function(x, by =1, digits = 0) {
# generate breaks
wsk = boxplot.stats(x)
brk = seq.int(wsk$stats[1]-by,wsk$stats[5],by=by)
brk = round(brk,digits = digits)

# generate labs
lab = NULL
for (n in seq(length(brk)-1)){
  lb = paste(brk[n]," to ",brk[n+1],sep="")
  lab = c(lab,lb)
}

# lab = c(lab,paste(brk[length(brk)]," < ",sep="")) 

# generate cuts
cut(x, breaks = brk, labels = lab, right = FALSE)
}
