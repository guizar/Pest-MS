# Initial variables

rm(list=ls())

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wdtables = "~/R/Pest-MS/tables"
wdfun = "~/R/Pest-MS/fun"
wddata ="~/R/Pest-MS/data/"

# functions
# source(file.path(wdfun,"color bar.R"))
# source(file.path(wdfun,"summarySE function.R"))
# source(file.path(wdfun,"multiplot.R"))


# ---- 1) Input files from Matlab -  Focusing on 2 degrees ----

setwd(paste(wddata,"Impact_2degC",sep="/"))

IPM_M2<-read.table("Impact_popmet_lres_a100b1p2_maize.txt", header= FALSE) 
IPM_M3<-read.table("Impact_popmet_lres_a100b1p3_maize.txt", header= FALSE) 
IPM_M4<-read.table("Impact_popmet_lres_a100b1p4_maize.txt", header= FALSE)

IPM_R2<-read.table("Impact_popmet_lres_a100b1p2_rice.txt", header= FALSE) 
IPM_R3<-read.table("Impact_popmet_lres_a100b1p3_rice.txt", header= FALSE) 
IPM_R4<-read.table("Impact_popmet_lres_a100b1p4_rice.txt", header= FALSE)

IPM_W2<-read.table("Impact_popmet_lres_a100b1p2_wheat.txt", header= FALSE) 
IPM_W3<-read.table("Impact_popmet_lres_a100b1p3_wheat.txt", header= FALSE) 
IPM_W4<-read.table("Impact_popmet_lres_a100b1p4_wheat.txt", header= FALSE)


IY_M2<-read.table("Impact_Yield_lres_a100b1p2_maize.txt", header= FALSE) 
IY_M3<-read.table("Impact_Yield_lres_a100b1p3_maize.txt", header= FALSE) 
IY_M4<-read.table("Impact_Yield_lres_a100b1p4_maize.txt", header= FALSE)

IY_R2<-read.table("Impact_Yield_lres_a100b1p2_rice.txt", header= FALSE) 
IY_R3<-read.table("Impact_Yield_lres_a100b1p3_rice.txt", header= FALSE) 
IY_R4<-read.table("Impact_Yield_lres_a100b1p4_rice.txt", header= FALSE)

IY_W2<-read.table("Impact_Yield_lres_a100b1p2_wheat.txt", header= FALSE) 
IY_W3<-read.table("Impact_Yield_lres_a100b1p3_wheat.txt", header= FALSE) 
IY_W4<-read.table("Impact_Yield_lres_a100b1p4_wheat.txt", header= FALSE)

MET_M2<-read.table("Impact_metabl_lres_a100b1p2_maize.txt", header= FALSE) 
MET_M3<-read.table("Impact_metabl_lres_a100b1p3_maize.txt", header= FALSE) 
MET_M4<-read.table("Impact_metabl_lres_a100b1p4_maize.txt", header= FALSE)

MET_R2<-read.table("Impact_metabl_lres_a100b1p2_rice.txt", header= FALSE) 
MET_R3<-read.table("Impact_metabl_lres_a100b1p3_rice.txt", header= FALSE) 
MET_R4<-read.table("Impact_metabl_lres_a100b1p4_rice.txt", header= FALSE)

MET_W2<-read.table("Impact_metabl_lres_a100b1p2_wheat.txt", header= FALSE) 
MET_W3<-read.table("Impact_metabl_lres_a100b1p3_wheat.txt", header= FALSE) 
MET_W4<-read.table("Impact_metabl_lres_a100b1p4_wheat.txt", header= FALSE)

POP_M2<-read.table("Impact_popul_lres_a100b1p2_maize.txt", header= FALSE) 
POP_M3<-read.table("Impact_popul_lres_a100b1p3_maize.txt", header= FALSE) 
POP_M4<-read.table("Impact_popul_lres_a100b1p4_maize.txt", header= FALSE)

POP_R2<-read.table("Impact_popul_lres_a100b1p2_rice.txt", header= FALSE) 
POP_R3<-read.table("Impact_popul_lres_a100b1p3_rice.txt", header= FALSE) 
POP_R4<-read.table("Impact_popul_lres_a100b1p4_rice.txt", header= FALSE)

POP_W2<-read.table("Impact_popul_lres_a100b1p2_wheat.txt", header= FALSE) 
POP_W3<-read.table("Impact_popul_lres_a100b1p3_wheat.txt", header= FALSE) 
POP_W4<-read.table("Impact_popul_lres_a100b1p4_wheat.txt", header= FALSE)



CLF_M<-read.table("Crop_loss_frac_maize.txt", header= FALSE) #crop fractional loss due to pests currently (from...)
CLF_R<-read.table("Crop_loss_frac_rice.txt", header= FALSE) #crop fractional loss due to pests currently (from...)
CLF_W<-read.table("Crop_loss_frac_wheat.txt", header= FALSE) #crop fractional loss due to pests currently (from...)

CL_M<-read.table("Crop_loss_maize.txt", header= FALSE)#crop loss (in tonnes per hectare of grid area).
CL_R<-read.table("Crop_loss_rice.txt", header= FALSE)#crop loss (in tonnes per hectare of grid area).
CL_W<-read.table("Crop_loss_wheat.txt", header= FALSE)#crop loss (in tonnes per hectare of grid area).

CY_M<-read.table("Crop_yield_maize.txt", header= FALSE) #crop yield (tonnnes per hectare of grid area) 
CY_R<-read.table("Crop_yield_rice.txt", header= FALSE) #crop yield (tonnnes per hectare of grid area) 
CY_W<-read.table("Crop_yield_wheat.txt", header= FALSE) #crop yield (tonnnes per hectare of grid area) 

CA_M<-read.table("Crop_area_maize.txt", header=FALSE) #  LAT LONG and % of grid cell planted  
CA_R<-read.table("Crop_area_rice.txt", header=FALSE) #  LAT LONG and % of grid cell planted  
CA_W<-read.table("Crop_area_wheat.txt", header=FALSE) #  LAT LONG and % of grid cell planted  

CGS_M<-read.table("Crop_grow_maize.txt", header=FALSE) #current growing season for maize.
CGS_R<-read.table("Crop_grow_rice.txt", header=FALSE) #current growing season for maize.
CGS_W<-read.table("Crop_grow_wheat.txt", header= FALSE) #current growing season for maize.

# ---- 2) Remove NAs in imported DF ----
ls = ls()
ls = ls[-which(ls %in% ls(pattern = "^wd") == T,arr.ind=TRUE)]


for (i in ls){
  assign(i, get(i)[complete.cases(i),])
}

rm(ls)
rm(i)

# ---- 3) Rename Variables ----
library(gdata)

# IPM=Insect Population Metabolism (change in pest pressure)
IPM_M2<- rename.vars(IPM_M2, c("V1","V2","V3"), c("LAT","LON","IPM_M2"))
IPM_M3<- rename.vars(IPM_M3, c("V1","V2","V3"), c("LAT","LON","IPM_M3"))
IPM_M4<- rename.vars(IPM_M4, c("V1","V2","V3"), c("LAT","LON","IPM_M4"))

IPM_R2<- rename.vars(IPM_R2, c("V1","V2","V3"), c("LAT","LON","IPM_R2"))
IPM_R3<- rename.vars(IPM_R3, c("V1","V2","V3"), c("LAT","LON","IPM_R3"))
IPM_R4<- rename.vars(IPM_R4, c("V1","V2","V3"), c("LAT","LON","IPM_R4"))

IPM_W2<- rename.vars(IPM_W2, c("V1","V2","V3"), c("LAT","LON","IPM_W2"))
IPM_W3<- rename.vars(IPM_W3, c("V1","V2","V3"), c("LAT","LON","IPM_W3"))
IPM_W4<- rename.vars(IPM_W4, c("V1","V2","V3"), c("LAT","LON","IPM_W4"))

IY_M2<- rename.vars(IY_M2, c("V1","V2","V3"), c("LAT","LON","IY_M2"))
IY_M3<- rename.vars(IY_M3, c("V1","V2","V3"), c("LAT","LON","IY_M3"))
IY_M4<- rename.vars(IY_M4, c("V1","V2","V3"), c("LAT","LON","IY_M4"))

IY_R2<- rename.vars(IY_R2, c("V1","V2","V3"), c("LAT","LON","IY_R2"))
IY_R3<- rename.vars(IY_R3, c("V1","V2","V3"), c("LAT","LON","IY_R3"))
IY_R4<- rename.vars(IY_R4, c("V1","V2","V3"), c("LAT","LON","IY_R4"))

IY_W2<- rename.vars(IY_W2, c("V1","V2","V3"), c("LAT","LON","IY_W2"))
IY_W3<- rename.vars(IY_W3, c("V1","V2","V3"), c("LAT","LON","IY_W3"))
IY_W4<- rename.vars(IY_W4, c("V1","V2","V3"), c("LAT","LON","IY_W4"))

CLF_M<- rename.vars(CLF_M, c("V1","V2","V3"), c("LAT","LON","CLF_M"))
CLF_R<- rename.vars(CLF_R, c("V1","V2","V3"), c("LAT","LON","CLF_R"))
CLF_W<- rename.vars(CLF_W, c("V1","V2","V3"), c("LAT","LON","CLF_W"))

#CLF_M: crop fractional loss due to pests currently for each crop
#Rescaling to percentage
CLF_M$CLF_M<-CLF_M$CLF_M*.01
CLF_R$CLF_R<-CLF_R$CLF_R*.01
CLF_W$CLF_W<-CLF_W$CLF_W*.01

#CL_M:crop loss (in tonnes per hectare) for each crop.
CL_M<- rename.vars(CL_M, c("V1","V2","V3"), c("LAT","LON","CL_M"))
CL_R<- rename.vars(CL_R, c("V1","V2","V3"), c("LAT","LON","CL_R"))
CL_W<- rename.vars(CL_W, c("V1","V2","V3"), c("LAT","LON","CL_W"))

#CY_x:crop yield (tonnnes per hectare) for crop x
CY_M<-rename.vars(CY_M, c("V1","V2","V3"), c("LAT","LON","CY_M"))
CY_R<-rename.vars(CY_R, c("V1","V2","V3"), c("LAT","LON","CY_R"))
CY_W<-rename.vars(CY_W, c("V1","V2","V3"), c("LAT","LON","CY_W"))

#CA_x:percent of grid area planted in crop, for crop x
CA_M<-rename.vars(CA_M, c("V1","V2","V3"), c("LAT","LON","CA_M"))
CA_R<-rename.vars(CA_R, c("V1","V2","V3"), c("LAT","LON","CA_R"))
CA_W<-rename.vars(CA_W, c("V1","V2","V3"), c("LAT","LON","CA_W"))

#CGS_M:current crop growing season (in days) for each crop.
CGS_M<- rename.vars(CGS_M, c("V1","V2","V3","V4","V5"), c("LAT","LON","CGS_M","PLT_M","HRVST_M"))
CGS_R<- rename.vars(CGS_R, c("V1","V2","V3","V4","V5"), c("LAT","LON","CGS_R","PLT_R","HRVST_R"))
CGS_W<- rename.vars(CGS_W, c("V1","V2","V3","V4","V5"), c("LAT","LON","CGS_W","PLT_W","HRVST_W"))

# Rename MET
MET_M2<- rename.vars(MET_M2, c("V1","V2","V3"), c("LAT","LON","MET_M2"))
MET_M3<- rename.vars(MET_M3, c("V1","V2","V3"), c("LAT","LON","MET_M3"))
MET_M4<- rename.vars(MET_M4, c("V1","V2","V3"), c("LAT","LON","MET_M4"))

MET_R2<- rename.vars(MET_R2, c("V1","V2","V3"), c("LAT","LON","MET_R2"))
MET_R3<- rename.vars(MET_R3, c("V1","V2","V3"), c("LAT","LON","MET_R3"))
MET_R4<- rename.vars(MET_R4, c("V1","V2","V3"), c("LAT","LON","MET_R4"))

MET_W2<- rename.vars(MET_W2, c("V1","V2","V3"), c("LAT","LON","MET_W2"))
MET_W3<- rename.vars(MET_W3, c("V1","V2","V3"), c("LAT","LON","MET_W3"))
MET_W4<- rename.vars(MET_W4, c("V1","V2","V3"), c("LAT","LON","MET_W4"))

POP_M2<- rename.vars(POP_M2, c("V1","V2","V3"), c("LAT","LON","POP_M2"))
POP_M3<- rename.vars(POP_M3, c("V1","V2","V3"), c("LAT","LON","POP_M3"))
POP_M4<- rename.vars(POP_M4, c("V1","V2","V3"), c("LAT","LON","POP_M4"))

POP_R2<- rename.vars(POP_R2, c("V1","V2","V3"), c("LAT","LON","POP_R2"))
POP_R3<- rename.vars(POP_R3, c("V1","V2","V3"), c("LAT","LON","POP_R3"))
POP_R4<- rename.vars(POP_R4, c("V1","V2","V3"), c("LAT","LON","POP_R4"))

POP_W2<- rename.vars(POP_W2, c("V1","V2","V3"), c("LAT","LON","POP_W2"))
POP_W3<- rename.vars(POP_W3, c("V1","V2","V3"), c("LAT","LON","POP_W3"))
POP_W4<- rename.vars(POP_W4, c("V1","V2","V3"), c("LAT","LON","POP_W4"))

# ---- 4) Merge DFs ----
#Next map metabolism and population (pick one for main graphics)
METPOP<-merge(MET_M2,MET_M3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,MET_M4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

METPOP<-merge(METPOP,MET_R2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,MET_R3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,MET_R4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

METPOP<-merge(METPOP,MET_W2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,MET_W3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,MET_W4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

METPOP<-merge(METPOP,POP_M2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_M3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_M4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

METPOP<-merge(METPOP,POP_R2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_R3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_R4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

METPOP<-merge(METPOP,POP_W2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_W3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
METPOP<-merge(METPOP,POP_W4, by=c("LAT","LON"),all=TRUE,sort=TRUE)


#### NOW MERGE IPM, IY and COMBINE WITH OTHER VARS
IPM<-merge(IPM_M2,IPM_M3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPM<-merge(IPM,IPM_M4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IPM<-merge(IPM,IPM_R2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPM<-merge(IPM,IPM_R3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPM<-merge(IPM,IPM_R4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IPM<-merge(IPM,IPM_W2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPM<-merge(IPM,IPM_W3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPM<-merge(IPM,IPM_W4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IY<-merge(IY_M2,IY_M3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IY<-merge(IY,IY_M4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IY<-merge(IY,IY_R2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IY<-merge(IY,IY_R3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IY<-merge(IY,IY_R4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IY<-merge(IY,IY_W2, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IY<-merge(IY,IY_W3, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IY<-merge(IY,IY_W4, by=c("LAT","LON"),all=TRUE,sort=TRUE)

IPMIY<-merge(IPM,IY, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CLF_M, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CLF_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CLF_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CL_M, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CL_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CL_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CY_M, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CY_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CY_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CA_M, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CA_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CA_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CGS_M, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CGS_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
IPMIY<-merge(IPMIY,CGS_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)

ALL<-merge(IPMIY,METPOP, by=c("LAT","LON"),all=TRUE,sort=TRUE)

# ---- 5) Delete merged objects ----
rm(list =c("IPM_M2","IPM_M3","IPM_M4","IPM_R2","IPM_R3","IPM_R4","IPM_W2","IPM_W3","IPM_W4","IY_M2","IY_M3","IY_M4","IY_R2","IY_R3","IY_R4","IY_W2","IY_W3","IY_W4","MET_M2","MET_M3","MET_M4","MET_R2","MET_R3","MET_R4","MET_W2","MET_W3","MET_W4","POP_M2","POP_M3","POP_M4","POP_R2","POP_R3","POP_R4","POP_W2","POP_W3","POP_W4","CLF_M","CLF_R","CLF_W","CL_M","CL_R","CL_W","CY_M","CY_R","CY_W","CA_M","CA_R","CA_W","CGS_M","CGS_R","CGS_W","IPM","IY"))

# ---- 6) ADD COUNTRY NAMES TO DATA ----
# to be used to make summary info for regions.
library(sp)
setwd(wddata)
load("WorldPolyCountries.Rdata")

coordinates(ALL) <- c("LON", "LAT")
# proj4string(ALL) = proj4string(WorldPolyCountries)
# identify countries
ALL <- data.frame(ALL, NAME=overlay(WorldPolyCountries, ALL)$NAME)

#Split the data frame up - you could probably do this in one go but this makes things a bit more obvious:
QUERY = ALL[is.na(ALL$NAME),]
NAMES = ALL[!is.na(ALL$NAME),]

library(FNN) # used to find nearest point and assign values
NEIGHS = get.knnx(NAMES[,c("LAT","LON")],QUERY[,c("LAT","LON")],k=1)
#Now insert the replacement colours directly into the data dataframe:
ALL[is.na(ALL$NAME),"NAME"]=NAMES$NAME[NEIGHS$nn.index]
rm(NAMES)
rm(QUERY)
rm(NEIGHS)

# --- 7) Calculate Grid area ----
# Isolate LON LAT from ALL
XY = ALL[,which(c("LON","LAT") %in% names(ALL))]

# are XY coordinates UNIQUE in XY object?
summary(unique(XY[,c(1,2)]) == XY[]) ## yes

# convert XY to SPDF
XY$v = 1
coordinates(XY) = ~LON+LAT

# rasterize XY
library(raster)
ras = raster() # 1 x 1 grid by default
ras = rasterize(XY, ras, field = rep(1, nrow(XY)))

# calculate area from raster object
AREA = area(ras)
AREA[coordinates(!is.na(ras)),] # prints area of non-NA cells

# pass area values to the XY SPDF obj
XY$AREA = extract(AREA,XY)

# convert XY back to DF, test whether LON|LAT cols are IDENTICAL in both XY and ALL objects
XY = as.data.frame(XY)
colnames(XY) = c("LON","LAT","v","AREA")
identical(XY[,c("LON","LAT")],ALL[,c("LON","LAT")]) # TRUE

# pass area values from XY to ALL object. Since area() returns values in KM2, multiply by 100 to convert to HA
ALL$AREA = XY$AREA * 100

# remove unused
rm(XY)
rm(ras)
rm(AREA)

# ---- 8) Merge Country, region and subregion data ----
rm(WorldPolyCountries)

# add counting column
ALL$CELLS<-1

# ADD UN REGIONS
# load country regions and codes
library(countrycode)
data=countrycode_data

# can they be simply merged?
summary(ALL$NAME %in% data$country.name) ## too many missing countries

# looks better if you run countrycode()
summary(is.na(countrycode(ALL$NAME,origin = "country.name",destination = "un")))


# add UN code number to ALL
ALL$code = countrycode(ALL$NAME,origin = "country.name",destination = "un", warn = T)

# which countries are still mising?
ALL$NAME[is.na(ALL$code)]

# add missing country code (to identify region)
ALL$code[is.na(ALL$code)] = unique(ALL$code[ALL$NAME=="China"])

# merge by UN code
t = merge(ALL,data,by.x= "code",by.y = "un",all.x=T,sort=T)

# are there any missing regions?
summary(is.na(t$region)) # looks good to go
unique(t$region)

# correct country.name (Taiwan)
# t$country.name[t$NAME=="Taiwan"] ="Taiwan"

# remove unused cols
t = t[,-which(names(t) %in% c("cowc","cown","fao","fips104","imf","ioc","iso2c","iso3c","iso3n","un","wb","regex","code"))]

# restore and remove
ALL = t
rm(t)
rm(data)
# ///////////////////

# ---- 9) CREATING SUMMARY VARIABLES FOR GRAPHING ----
ALL$MET_AVG2<-rowMeans(ALL[ ,c("MET_M2","MET_R2","MET_W2")],na.rm=TRUE)
ALL$MET_AVG3<-rowMeans(ALL[ ,c("MET_M3","MET_R3","MET_W3")],na.rm=TRUE)
ALL$MET_AVG4<-rowMeans(ALL[ ,c("MET_M4","MET_R4","MET_W4")],na.rm=TRUE)

ALL$POP_AVG2<-rowMeans(ALL[ ,c("POP_M2","POP_R2","POP_W2")],na.rm=TRUE)
ALL$POP_AVG3<-rowMeans(ALL[ ,c("POP_M3","POP_R3","POP_W3")],na.rm=TRUE)
ALL$POP_AVG4<-rowMeans(ALL[ ,c("POP_M4","POP_R4","POP_W4")],na.rm=TRUE)

ALL$IPM_AVG2<-rowMeans(ALL[ ,c("IPM_M2","IPM_R2","IPM_W2")],na.rm=TRUE)
ALL$IPM_AVG3<-rowMeans(ALL[ ,c("IPM_M3","IPM_R3","IPM_W3")],na.rm=TRUE)
ALL$IPM_AVG4<-rowMeans(ALL[ ,c("IPM_M4","IPM_R4","IPM_W4")],na.rm=TRUE)

# ---- 10) CALCULATE other summary variables for graphs and tables ----
########YIELD PER HA PLANTED (YLD_HA_x)
ALL$YLD_HA_M<-ALL$CY_M/ALL$CA_M
ALL$YLD_HA_R<-ALL$CY_R/ALL$CA_R
ALL$YLD_HA_W<-ALL$CY_W/ALL$CA_W

########YIELD (tonnes) TOTAL PER CELL (YLD_TOT_x)
ALL$YLD_TOT_M<-ALL$CY_M*ALL$AREA
ALL$YLD_TOT_R<-ALL$CY_R*ALL$AREA
ALL$YLD_TOT_W<-ALL$CY_W*ALL$AREA

########TOTAL TONNES LOST PER CELL DUE TO CLIMATE INDUCED CHANGES IN INSECT PESTS (CROP LOSS PER CELL- CL2050_xy)
ALL$CL2050_M2<-ALL$IPM_M2*ALL$CLF_M*ALL$CY_M*ALL$AREA
ALL$CL2050_M3<-ALL$IPM_M3*ALL$CLF_M*ALL$CY_M*ALL$AREA
ALL$CL2050_M4<-ALL$IPM_M4*ALL$CLF_M*ALL$CY_M*ALL$AREA

ALL$CL2050_R2<-ALL$IPM_R2*ALL$CLF_R*ALL$CY_R*ALL$AREA
ALL$CL2050_R3<-ALL$IPM_R3*ALL$CLF_R*ALL$CY_R*ALL$AREA
ALL$CL2050_R4<-ALL$IPM_R4*ALL$CLF_R*ALL$CY_R*ALL$AREA

ALL$CL2050_W2<-ALL$IPM_W2*ALL$CLF_W*ALL$CY_W*ALL$AREA
ALL$CL2050_W3<-ALL$IPM_W3*ALL$CLF_W*ALL$CY_W*ALL$AREA
ALL$CL2050_W4<-ALL$IPM_W4*ALL$CLF_W*ALL$CY_W*ALL$AREA

####### TOTAL EXPECTED PROPORTION OF CROP YIELD LOST  PER HA 2  degree worldDUE TO PESTS (CURRENT+ADDEED): CLP2050_xy
ALL$CLP2050_M2<-ALL$CLF_M+ALL$IPM_M2*ALL$CLF_M
ALL$CLP2050_M3<-ALL$CLF_M+ALL$IPM_M3*ALL$CLF_M
ALL$CLP2050_M4<-ALL$CLF_M+ALL$IPM_M4*ALL$CLF_M

ALL$CLP2050_R2<-ALL$CLF_R+ALL$IPM_R2*ALL$CLF_R
ALL$CLP2050_R3<-ALL$CLF_R+ALL$IPM_R3*ALL$CLF_R
ALL$CLP2050_R4<-ALL$CLF_R+ALL$IPM_R4*ALL$CLF_R

ALL$CLP2050_W2<-ALL$CLF_W+ALL$IPM_W2*ALL$CLF_W
ALL$CLP2050_W3<-ALL$CLF_W+ALL$IPM_W3*ALL$CLF_W
ALL$CLP2050_W4<-ALL$CLF_W+ALL$IPM_W4*ALL$CLF_W

#CURRENT CROP YIELD PER HA PLANTED
ALL$CYH_M<-ALL$CY_M/ALL$CA_M
ALL$CYH_R<-ALL$CY_R/ALL$CA_R
ALL$CYH_W<-ALL$CY_W/ALL$CA_W

#PROJECTED FUTURE CROP YIELD PER HA PLANTED
ALL$FCY_M2<-ALL$CYH_M*(1-(ALL$CLF_M*ALL$IPM_M2))
ALL$FCY_M3<-ALL$CYH_M*(1-(ALL$CLF_M*ALL$IPM_M3))
ALL$FCY_M4<-ALL$CYH_M*(1-(ALL$CLF_M*ALL$IPM_M4))

ALL$FCY_R2<-ALL$CYH_R*(1-(ALL$CLF_R*ALL$IPM_R2))
ALL$FCY_R3<-ALL$CYH_R*(1-(ALL$CLF_R*ALL$IPM_R3))
ALL$FCY_R4<-ALL$CYH_R*(1-(ALL$CLF_R*ALL$IPM_R4))

ALL$FCY_W2<-ALL$CYH_W*(1-(ALL$CLF_W*ALL$IPM_W2))
ALL$FCY_W3<-ALL$CYH_W*(1-(ALL$CLF_W*ALL$IPM_W3))
ALL$FCY_W4<-ALL$CYH_W*(1-(ALL$CLF_W*ALL$IPM_W4))

########### YIELD LOSS DO TO CLIMATE CHANGE: IYCC_xy)
# calculated from current crop loss (frational) and deltaIPM.
#example: If current crop loss is 10%, and insect population metabolism is going up 20% in an area,
#then IYCC = 0.02, meaning we expect tthe loss of an additional 2% of crop yield due to insects in 75 years.  
ALL$IYCC_M2<-ALL$IPM_M2*(ALL$CLF_M)
ALL$IYCC_M3<-ALL$IPM_M3*(ALL$CLF_M)
ALL$IYCC_M4<-ALL$IPM_M4*(ALL$CLF_M)

ALL$IYCC_R2<-ALL$IPM_R2*(ALL$CLF_R)
ALL$IYCC_R3<-ALL$IPM_R3*(ALL$CLF_R)
ALL$IYCC_R4<-ALL$IPM_R4*(ALL$CLF_R)

ALL$IYCC_W2<-ALL$IPM_W2*(ALL$CLF_W)
ALL$IYCC_W3<-ALL$IPM_W3*(ALL$CLF_W)
ALL$IYCC_W4<-ALL$IPM_W4*(ALL$CLF_W)

# ---- 11) Bounce master CSV and save RData ----
# CUT ALL ROWS WITH NO LAT AND LONG
summary(is.na(ALL$LAT))
ALL<-subset(ALL,!is.na(ALL$LAT))

setwd(wdtables)

# Write ALL - name it to 2c degrees
# write.csv(ALL, file.path(wdtables,"ALL_2c.csv"),row.names=FALSE)
ALL_2c = ALL
rm(ALL)

# rename and remove
rm(IPMIY)
rm(METPOP)

# Write Rdata
setwd(wd)
rm(list=c("wd","wdpng","wdtables","wdfun","wddata"))
save.image(file ="RData/ALL_2c.RData")
# ///////////////////////

