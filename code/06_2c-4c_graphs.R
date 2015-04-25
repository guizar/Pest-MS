# Load data and initial variables ------

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"

# Load data
setwd(wdrdata)
load("ALL_2c.RData")
load("ALL_4c.RData")

# load libraries used to produce ALL graphs
library(reshape2)
library(classInt)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(R.matlab)
library(arrayhelpers)

# Load functions
source(file.path(wdfun,"cols_gg.r")) 
source(file.path(wdfun,"multiplots_gg.r"))
source(file.path(wdfun,"custom_cut.r"))


# Fig 3) 2 to 4 degree bar for IPM for each crop for each UN region.----
# Variables
# | IPM 2 deg| IPM 4 deg | 
# regions
# crops
# Phi = 0.001
# pending: need to get 4 deg values

