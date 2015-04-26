# Clean space =  leabe only ALL_ and DAT_objects

rm(list = ls()[
  -which(ls() %in% c(ls(pattern ="ALL_"),ls(pattern ="DAT_")))])

# restore paths
wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wddata = "~/R/Pest-MS/data/"
wdrdata = "~/R/Pest-MS/RData/"
wdfun = "~/R/Pest-MS/fun"