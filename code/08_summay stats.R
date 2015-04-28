library(doBy)

# Load data
load(file.path(wdrdata,"DAT_2c.RData"))
load(file.path(wdrdata,"DAT_4c.RData"))

#median values (global) - 2 degree
summaryBy(value ~ crop + phi + fact, data = DAT_2c, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )
#median values (global) - 4 degree
summaryBy(value ~ crop + phi + fact, data = DAT_4c, FUN = function(x) { c(med = median(x), mean = mean(x), sd = sd(x)) } )

