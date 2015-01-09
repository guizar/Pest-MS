# Initial variables

wd = "~/R/Pest-MS/"
wdpng = "~/R/Pest-MS/png"
wddata = "~/R/Pest-MS/data"
wdfun = "~/R/Pest-MS/functions"

# functions
source(file.path(wdfun,"color bar.R"))
source(file.path(wdfun,"summarySE function.R"))
source(file.path(wdfun,"multiplot.R"))


#Input files from Matlab -  Focusing on 2 degrees

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


# Renaming variables
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

# ////////// Comments //////////
# JTT: I THINK THE CGS DATA ARE NOT CORRECT COMING IN. PLOT THEM...note, we do not expand the growing season with temp change in our model.

# AG re(JJT): How should the CGS data look like? if you try summarizing CGS_x$CGS_x you will see no NA's, however both PLT_x and HRVST_x have a lot of NA's - what are these 2 last variables?
# summary(CGS_M)
# # ///////////////////////////

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

############################ MERGING TO ALL TO FIX HOTSPOTS

ALL<-merge(IPMIY,METPOP, by=c("LAT","LON"),all=TRUE,sort=TRUE)

# /////////////////////////////
# hotspots are areas where IPM and POP predicts insects eating more crop than is in the field.
# this problem is caused primarily by poor data resolution in the crop growing season variable. If we had better data there, we would likely NOT have many cases like these.
# /////////////////////////////

# start by breaking into crops to remove NAs.

ALL_M<-subset(ALL,!is.na(ALL$IPM_M2),select=c(LAT,LON,IPM_M2,IPM_M3,IPM_M4,IY_M2,IY_M3,IY_M4,CLF_M,CL_M,CY_M,CA_M,CGS_M,MET_M2,MET_M3,MET_M4,POP_M2,POP_M3,POP_M4))

# now reset values of IPM in which projected needs is greater than crops present, so that they limit at 

#crops present
ALL_M$IPM_M2[ALL_M$IPM_M2*ALL_M$CLF_M>1]<-1/(ALL_M$CLF_M)[ALL_M$IPM_M2*ALL_M$CLF_M>1]  
ALL_M$IPM_M3[ALL_M$IPM_M3*ALL_M$CLF_M>1]<-1/(ALL_M$CLF_M)[ALL_M$IPM_M3*ALL_M$CLF_M>1]
ALL_M$IPM_M4[ALL_M$IPM_M4*ALL_M$CLF_M>1]<-1/(ALL_M$CLF_M)[ALL_M$IPM_M4*ALL_M$CLF_M>1]

#now reset the values of POP as the difference between the NEW IPM and MET (caps pop at level that would sustain IPM)
ALL_M$POP_M2[ALL_M$IPM_M2*ALL_M$CLF_M1>1]<-(ALL_M$IPM_M2-ALL_M$MET_M2)[ALL_M$IPM_M2*ALL_M$CLF_M>1]
ALL_M$POP_M3[ALL_M$IPM_M3*ALL_M$CLF_M1>1]<-(ALL_M$IPM_M3-ALL_M$MET_M3)[ALL_M$IPM_M3*ALL_M$CLF_M>1]
ALL_M$POP_M4[ALL_M$IPM_M4*ALL_M$CLF_M>1]<-(ALL_M$IPM_M4-ALL_M$MET_M4)[ALL_M$IPM_M4*ALL_M$CLF_M>1]

#Now repeat for rice and wheat
ALL_R<-subset(ALL,!is.na(ALL$IPM_R2),select=c(LAT,LON,IPM_R2,IPM_R3,IPM_R4, IY_R2,IY_R3,IY_R4,CLF_R,CL_R,CY_R,CA_R,CGS_R,MET_R2,MET_R3,MET_R4,POP_R2,POP_R3,POP_R4))

ALL_R$IPM_R2[ALL_R$IPM_R2*ALL_R$CLF_R>1]<-1/(ALL_R$CLF_R)[ALL_R$IPM_R2*ALL_R$CLF_R>1]  
ALL_R$IPM_R3[ALL_R$IPM_R3*ALL_R$CLF_R>1]<-1/(ALL_R$CLF_R)[ALL_R$IPM_R3*ALL_R$CLF_R>1]
ALL_R$IPM_R4[ALL_R$IPM_R4*ALL_R$CLF_R>1]<-1/(ALL_R$CLF_R)[ALL_R$IPM_R4*ALL_R$CLF_R>1]

ALL_R$POP_R2[ALL_R$IPM_R2*ALL_R$CLF_R1>1]<-(ALL_R$IPM_R2-ALL_R$MET_R2)[ALL_R$IPM_R2*ALL_R$CLF_R>1]
ALL_R$POP_R3[ALL_R$IPM_R3*ALL_R$CLF_R1>1]<-(ALL_R$IPM_R3-ALL_R$MET_R3)[ALL_R$IPM_R3*ALL_R$CLF_R>1]
ALL_R$POP_R4[ALL_R$IPM_R4*ALL_R$CLF_R>1]<-(ALL_R$IPM_R4-ALL_R$MET_R4)[ALL_R$IPM_R4*ALL_R$CLF_R>1]

ALL_W<-subset(ALL,!is.na(ALL$IPM_W2),select=c(LAT,LON,IPM_W2,IPM_W3,IPM_W4, IY_W2,IY_W3,IY_W4,CLF_W,CL_W,CY_W,CA_W,CGS_W,MET_W2,MET_W3,MET_W4,POP_W2,POP_W3,POP_W4))

ALL_W$IPM_W2[ALL_W$IPM_W2*ALL_W$CLF_W>1]<-1/(ALL_W$CLF_W)[ALL_W$IPM_W2*ALL_W$CLF_W>1]  
ALL_W$IPM_W3[ALL_W$IPM_W3*ALL_W$CLF_W>1]<-1/(ALL_W$CLF_W)[ALL_W$IPM_W3*ALL_W$CLF_W>1]
ALL_W$IPM_W4[ALL_W$IPM_W4*ALL_W$CLF_W>1]<-1/(ALL_W$CLF_W)[ALL_W$IPM_W4*ALL_W$CLF_W>1]

ALL_W$POP_W2[ALL_W$IPM_W2*ALL_W$CLF_W1>1]<-(ALL_W$IPM_W2-ALL_W$MET_W2)[ALL_W$IPM_W2*ALL_W$CLF_W>1]
ALL_W$POP_W3[ALL_W$IPM_W3*ALL_W$CLF_W1>1]<-(ALL_W$IPM_W3-ALL_W$MET_W3)[ALL_W$IPM_W3*ALL_W$CLF_W>1]
ALL_W$POP_W4[ALL_W$IPM_W4*ALL_W$CLF_W>1]<-(ALL_W$IPM_W4-ALL_W$MET_W4)[ALL_W$IPM_W4*ALL_W$CLF_W>1]

############################ NOW BRINGING THE FILES BACK TOGETHER

ALL<-merge(ALL_M,ALL_R, by=c("LAT","LON"),all=TRUE,sort=TRUE)
ALL<-merge(ALL,ALL_W, by=c("LAT","LON"),all=TRUE,sort=TRUE)

# /////////////////////////////
# CHECK TO MAKE SURE IPM an POP are now correct (turn on to check)
# M2prob <- subset(ALL,IPM_M2*(CLF_M)>1,select=c(LAT,LON,IPM_M2))
# M3prob <- subset(ALL,IPM_M3*(CLF_M)>1,select=c(LAT,LON,IPM_M3))
# M4prob <- subset(ALL,IPM_M4*(CLF_M)>1,select=c(LAT,LON,IPM_M4))

# It turns out that this procedure does not change the nature of the graphics at all, because it simply dampens the largest population
# increases, having no effect on medians, and no effect on colors (as I plot extreme values as all one color, and this just makes them less extreme)

# AG re (JJT): Not sure I'm following, but Mxprob objects result in empty files
# /////////////////////////////

############################ CREATING SUMMARY VARIABLES FOR GRAPHING

ALL$MET_AVG2<-rowMeans(ALL[ ,c("MET_M2","MET_R2","MET_W2")],na.rm=TRUE)
ALL$MET_AVG3<-rowMeans(ALL[ ,c("MET_M3","MET_R3","MET_W3")],na.rm=TRUE)
ALL$MET_AVG4<-rowMeans(ALL[ ,c("MET_M4","MET_R4","MET_W4")],na.rm=TRUE)

ALL$POP_AVG2<-rowMeans(ALL[ ,c("POP_M2","POP_R2","POP_W2")],na.rm=TRUE)
ALL$POP_AVG3<-rowMeans(ALL[ ,c("POP_M3","POP_R3","POP_W3")],na.rm=TRUE)
ALL$POP_AVG4<-rowMeans(ALL[ ,c("POP_M4","POP_R4","POP_W4")],na.rm=TRUE)

ALL$IPM_AVG2<-rowMeans(ALL[ ,c("IPM_M2","IPM_R2","IPM_W2")],na.rm=TRUE)
ALL$IPM_AVG3<-rowMeans(ALL[ ,c("IPM_M3","IPM_R3","IPM_W3")],na.rm=TRUE)
ALL$IPM_AVG4<-rowMeans(ALL[ ,c("IPM_M4","IPM_R4","IPM_W4")],na.rm=TRUE)



#################################### ADD COUNTRY NAMES TO DATA
# to be used to make summary info for regions.
setwd(wddata)
load("WorldPolyCountries.Rdata")

library(sp)
coordinates(ALL) <- c("LON", "LAT")

# identify countries
ALL <- data.frame(ALL, NAME=overlay(WorldPolyCountries, ALL)$NAME)

#Split the data frame up - you could probably do this in one go but this makes things a bit more obvious:
QUERY = ALL[is.na(ALL$NAME),]
NAMES = ALL[!is.na(ALL$NAME),]

library(FNN) # used to find nearest point and assign values
NEIGHS = get.knnx(NAMES[,c("LAT","LON")],QUERY[,c("LAT","LON")],k=1)
#Now insert the replacement colours directly into the data dataframe:
ALL[is.na(ALL$NAME),"NAME"]=NAMES$NAME[NEIGHS$nn.index]



############################### CALCULATING CELL AREA
ALL$AREA= ((6371^2)*abs((ALL$LON*pi/180)-((ALL$LON-0.1667)*pi/180))*abs(sin(ALL$LAT*pi/180)-sin((ALL$LAT-0.1667)*pi/180)))*100
#This is in HA (note the multiplication by 100.)

############################### CALCULATING other summary variables for graphs and tables


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

# Merge Country to region data -- MISSING CSV IN THIS SECTION
# setwd(wddata)
# C_TO_R<-data.frame(read.csv("COUNTRY_TO_REGION.csv", header=TRUE))
# ALL<-merge(C_TO_R,ALL,"NAME",all=TRUE,sort=TRUE)
# #add counting column
# ALL$CELLS<-1
# ///////////////////


# /////////////////// ADD UN REGIONS
# not used anymore: unsub = read.csv("unsubregions2.csv",header = T, stringsAsFactor=F) -  data from Ahttp://unstats.un.org/unsd/methods/m49/m49regin.htm

# load country regions and codes
library(countrycode)
data=countrycode_data

# can they be simply merged?
summary(ALL$NAME %in% data$country.name) ## too many missing countries

# looks better if you run countrycode()
summary(is.na(countrycode(ALL$NAME,origin = "country.name",destination = "un")))


# add UN code number to ALL
ALL$code = countrycode(ALL$NAME,origin = "country.name",destination = "un", warn = T)

# add missing country code (to identify region)
ALL$code[is.na(ALL$code)] = unique(ALL$code[ALL$NAME=="China"])
# merge by UN code
t = merge(ALL,data,by.x= "code",by.y = "un",all.x=T,sort=T)

# are there any missing regions?
summary(is.na(t$region)) # looks good to go
unique(t$region)

# correct country.name (Taiwan)
t$country.name[t$NAME=="Taiwan"] ="Taiwan"

# remove unused cols
t = t[,-which(names(t) %in% c("cowc","cown","fao","fips104","imf","ioc","iso2c","iso3c","iso3n","un","wb","regex"))]

# restore and remove
ALL = t
rm(t)
# ///////////////////



# /////////////////////// This is the master csv that will be used to produce the plots and tables

# CUT ALL ROWS WITH NO LAT AND LONG
summary(is.na(ALL$LAT))
ALL<-subset(ALL,!is.na(ALL$LAT))

setwd(wddata)

# Write ALL
write.table(ALL, file = "ALL.csv", sep = ",", col.names = NA)

# Load ALL
ALL<-read.csv("ALL.csv", header=TRUE,na.strings = c("#VALUE!", "#N/A", "N/A", "NA", ""))
# ///////////////////////



# Total impact on crop production : these produce outputs FOR PAPER, not in a table

# /////////////////////////////
# CLEAN:  THIS SECTION NEEDS TO BE CLEANED - NOT CORRECT - do not / by sum, div. by count, or take mean...
# AG re (JJT): Where are you allocating these vals? did you changed sum by count or mean?
# /////////////////////////////

sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M2,na.rm=T)/sum(ALL$YLD_TOT_M,na.rm=T)
sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R2,na.rm=T)/sum(ALL$YLD_TOT_R,na.rm=T)
sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W2,na.rm=T)/sum(ALL$YLD_TOT_W,na.rm=T)

sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M3,na.rm=T)/sum(ALL$YLD_TOT_M,na.rm=T)
sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R3,na.rm=T)/sum(ALL$YLD_TOT_R,na.rm=T)
sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W3,na.rm=T)/sum(ALL$YLD_TOT_W,na.rm=T)

sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M4,na.rm=T)/sum(ALL$YLD_TOT_M,na.rm=T)
sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R4,na.rm=T)/sum(ALL$YLD_TOT_R,na.rm=T)
sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W4,na.rm=T)/sum(ALL$YLD_TOT_W,na.rm=T)

sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M3,na.rm=T)
sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R3,na.rm=T)
sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W3,na.rm=T)


# Maize cost
#sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M2,na.rm=T)*321
#sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M3,na.rm=T)*321
#sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M4,na.rm=T)*321

# Rice cost
#sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R2,na.rm=T)*583
#sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R3,na.rm=T)*583
#sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R4,na.rm=T)*583

# Wheat cost
#sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W2,na.rm=T)*305
#sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W3,na.rm=T)*305
#sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W4,na.rm=T)*305



sum(ALL$YLD_TOT_M*ALL$CLF_M*ALL$IPM_M4,na.rm=T)/sum(ALL$YLD_TOT_M,na.rm=T)
sum(ALL$YLD_TOT_R*ALL$CLF_R*ALL$IPM_R4,na.rm=T)/sum(ALL$YLD_TOT_R,na.rm=T)
sum(ALL$YLD_TOT_W*ALL$CLF_W*ALL$IPM_W4,na.rm=T)/sum(ALL$YLD_TOT_W,na.rm=T)


# /////////////////////// AG: Graphics start here, I will comment everything from this point onwards and work on the plots in the graphs.r file

# #### SETUP MAPPING
# # world box
# xy.sp= cbind(ALL$LON, ALL$LAT)
# xy.cc= coordinates(xy.sp)
# bbox(xy.sp)
# #setup for plotting (making grid)
# grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = ALL[,c(3:110)], tolerance=0.75)
# #ONLY NEED THIS ONCE IF WE USE THE ALL FILE
# wrld.m= map("world", plot=FALSE, fill=FALSE)
# wrld.p= map2SpatialLines(wrld.m, proj4string=CRS("+proj=longlat +proj=lcc"))  #turn map into spatial lines
# wrld<- list("sp.lines", wrld.p, col="black")
# 
# ##
# #First a map of world to place UN subregion boundaries on
# world<-map("world")
# #now export that
# map('worldHires')
# ############ FIG 1 GRAPHICS
# 
# ############# FIG 1A: METABOLISM GRAPHIC
# quants<-c(-Inf,seq(0,.3,by=.02),Inf) #to specify quant cutoffs
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c('midnightblue',jet.colors,"#7F0000")
# txt0<-list("sp.text", c(20,-40), "MET_AVG")
# p_met <- spplot(grd1, 'MET_AVG3', at=quants, sp.layout=list(wrld,txt0), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9))
# 
# png(file="~/Dropbox/climate change/pest MS shared/METABOLISM.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p_met)
# dev.off()
# 
# ############# FIG 1B: POPULATION GRAPHIC
# txt0<-list("sp.text", c(20,-40), "POP_AVG")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# pop <- spplot(grd1, 'POP_AVG2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,
# colorkey=list(space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# png(file="~/Dropbox/climate change/pest MS shared/POPULATION AVG P3.png", width=900, height=1400) # CHECK RESOLUTION 
# print(pop)
# dev.off()
# 
# ############################# FIG 1C IPM AVERAGE ACROSS ALL THREE CROPS####
# txt0<-list("sp.text", c(20,-40), "IPM (change in Pest pressure)_AVG")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# pop <- spplot(grd1, 'IPM_AVG3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,
#               colorkey=list(space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# png(file="~/Dropbox/climate change/pest MS shared/IPM_AVG P3.png", width=900, height=1400) # CHECK RESOLUTION 
# print(pop)
# dev.off()
# 
# #############################LATITUDE MEDEAN LINES GRAPHIC
# #############################Fig 1A, left hand side: METABOLISM AND POP AT PHI=3
# #-180  to -20 = new world, -20 to 60 = europe, and 60 up = asia.
# #could create lines for different areas using ...subset(IPMIY,LON<-20)., in ddply, but starting for world.
# #round latitude to get medians for 1 degree belts
# ALL$LAT2<-round(ALL$LAT)
# #counting cells for ddply.
# ALL$CELLS_M<-ifelse(is.na(ALL$POP_M3),0,1)
# ALL$CELLS_R<-ifelse(is.na(ALL$POP_R3),0,1)
# ALL$CELLS_W<-ifelse(is.na(ALL$POP_W3),0,1)
# LATMEDALL<-ddply(ALL,.(LAT2),summarise,
#                  METM3=median(MET_M3,na.rm=TRUE),
#                  METR3=median(MET_R3,na.rm=TRUE),
#                  METW3=median(MET_W3,na.rm=TRUE),
#                  POPM3=median(POP_M3,na.rm=TRUE),
#                  POPR3=median(POP_R3,na.rm=TRUE),
#                  POPW3=median(POP_W3,na.rm=TRUE),
#                  IPMM3=median(IPM_M3,na.rm=TRUE),
#                  IPMR3=median(IPM_R3,na.rm=TRUE),
#                  IPMW3=median(IPM_W3,na.rm=TRUE),
#                  N_M=sum(CELLS_M),
#                  N_R=sum(CELLS_R),
#                  N_W=sum(CELLS_W))
# 
# LINE1<- melt(LATMEDALL[ ,c(1:4)], id="LAT2")  # convert to long format
# 
# ###ALEX - the next three graphs are kind of a pain - we need to cli
# 
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p5<-ggplot(data=LINE1,aes(x=LAT2, y=value,group=variable)) +
#   geom_line(aes(colour=variable), size = 1.5)+
#   coord_flip() +
#   scale_x_continuous(limits=c(-60,60),name="Latitude") +
#   theme_bw() +
#   scale_colour_manual(values=cbbPalette, 
#                       name="                 Crop",
#                       breaks=c("METM3", "METR3","METW3"),               
#                       labels=c("Maize", "Rice","Wheat")) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = -60)
# png(file="~/Dropbox/climate change/pest MS shared/LAT_MEAN_MET.png", width=400, height=1000) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# LINE1<- melt(LATMEDALL[ ,c(1,5:7)], id="LAT2")  # convert to long format
# 
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p5<-ggplot(data=LINE1,aes(x=LAT2, y=value,group=variable)) +
#   geom_line(aes(colour=variable), size = 1.5)+
#   coord_flip() +
#   scale_x_continuous(limits=c(-60,60),name="Latitude") +
#   theme_bw() +
#   scale_colour_manual(values=cbbPalette, 
#                       name="                 Crop",
#                       breaks=c("METM3", "METR3","METW3"),               
#                       labels=c("Maize", "Rice","Wheat")) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = -60)
# png(file="~/Dropbox/climate change/pest MS shared/LAT_MEAN_POP.png", width=400, height=1000) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# LINE1<- melt(LATMEDALL[ ,c(1,8:10)], id="LAT2")  # convert to long format
# 
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# p5<-ggplot(data=LINE1,aes(x=LAT2, y=value,group=variable)) +
#   geom_line(aes(colour=variable), size = 1.5)+
#   coord_flip() +
#   scale_x_continuous(limits=c(-60,60),name="Latitude") +
#   theme_bw() +
#   scale_colour_manual(values=cbbPalette, 
#                       name="                 Crop",
#                       breaks=c("METM3", "METR3","METW3"),               
#                       labels=c("Maize", "Rice","Wheat")) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = -60)
# png(file="~/Dropbox/climate change/pest MS shared/LAT_MEAN_IPM.png", width=400, height=1000) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# # that was kind of a hack job, but I can clean it up more if needed.
# 
# 
# 
# ################################IPM and YIELD GRAPHS FOR MAIN AND SUP.##################
# 
# #Fig 2 maps
# #######deltaIPM for each crop at phi = 0.001
# 
# 
# #Define quantiles
# #quants<-c(-1,seq(0,1,by=.025),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# #jet.colors <-
# #  colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# #cols<-c("grey",jet.colors,"#7F0000")
# 
# #wrld.m= map("world", plot=FALSE, fill=FALSE)
# #wrld.p= map2SpatialLines(wrld.m, proj4string=CRS("+proj=longlat +proj=lcc"))  #turn map into spatial lines
# #wrld<- list("sp.lines", wrld.p, col="darkgray")
# 
# ### 
# txt0<-list("sp.text", c(20,-40), "pest pressure change")  # WHY HIDE?
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001") # WHY HIDE?
# p2 <- spplot(grd1, 'IPM_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p3 <- spplot(grd1, 'IPM_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p4 <- spplot(grd1, 'IPM_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Pressure change_IPM_allthreecrops_phi001.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# # Supp Figures - all three crops run individually at all runs of PHI
# #Maize IPM all runs of PHI
# txt0<-list("sp.text", c(20,-40), "pest pressure change")
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p2 <- spplot(grd1, 'IPM_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p5 <- spplot(grd1, 'IPM_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p8 <- spplot(grd1, 'IPM_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p8,p5,p2, layout=c(2,4))
# png(file="~/Dropbox/climate change/pest MS shared/Pressure change_IPM_MAIZEfullrange.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# #RICE IPM all runs of PHI
# txt0<-list("sp.text", c(20,-40), "pest pressure change")
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p2 <- spplot(grd1, 'IPM_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p5 <- spplot(grd1, 'IPM_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p8 <- spplot(grd1, 'IPM_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p8,p5,p2, layout=c(2,4))
# png(file="~/Dropbox/climate change/pest MS shared/Pressure change_IPM_RICEfullrange.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# #Wheat IPM all runs of PHI
# txt0<-list("sp.text", c(20,-40), "pest pressure change")
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p2 <- spplot(grd1, 'IPM_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p5 <- spplot(grd1, 'IPM_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p8 <- spplot(grd1, 'IPM_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p8,p5,p2, layout=c(2,4))
# png(file="~/Dropbox/climate change/pest MS shared/Pressure change_IPM_WHEATfullrange.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ### Supplimental figures 9 panel 3 crops and 3 phi IPM
# 
# # Supp Figures - all three crops run individually at all runs of PHI
# #Maize IPM all runs of PHI
# txt0<-list("sp.text", c(20,-40), "pest pressure change")
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p1 <- spplot(grd1, 'IPM_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p2 <- spplot(grd1, 'IPM_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p3 <- spplot(grd1, 'IPM_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# 
# #RICE IPM all runs of PHI
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p4 <- spplot(grd1, 'IPM_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p5 <- spplot(grd1, 'IPM_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p6 <- spplot(grd1, 'IPM_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# 
# #Wheat IPM all runs of PHI
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p7 <- spplot(grd1, 'IPM_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p8 <- spplot(grd1, 'IPM_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p9 <- spplot(grd1, 'IPM_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p10 <- c(p7,p8,p9,p4,p5,p6,p1,p2,p3, layout=c(3,3))
# png(file="~/Dropbox/climate change/pest MS shared/Pressure change_IPM_3by3.png", width=2000, height=1500) # CHECK RESOLUTION 
# print(p10)
# dev.off()
# 
# ####
# 
# #############POPULATION GRAPHIC
# 
# ##########Define quantiles and Colors for log scale color ramps going positive and negative
# a<-lseq(.001,.5,length=22)
# quants<-c(rev(a*-1),0,a,Inf)
# 
# #ticks<-c(-1,-.5,-.1,0,.1,.5,1)
# #quants<-c(seq(-.5,.5,by=.025),Inf)
# 
# #specify colors and split in the center
# gray.colors<-
#   gray.colors(length(quants/2), start = 0.1, end = 0.9, gamma = 2.2)
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c(gray.colors,jet.colors,"#7F0000")
# txt0<-list("sp.text", c(20,-40), "POP_AVG")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_AVG2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3log <- spplot(grd1, 'POP_AVG3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_AVG4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3log,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POPAVG_log.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ################POP_MAIZE Log
# txt0<-list("sp.text", c(20,-40), "POP_MAIZE_log")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_MAIZE_log.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# ################POP_RICE Log
# txt0<-list("sp.text", c(20,-40), "POP_RICE_log")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# ### when ready to export:
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_RICE_log.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ################POP_WHEAT LOG
# txt0<-list("sp.text", c(20,-40), "POP_WHEAT_log")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# ### when ready to export:
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_WHEAT_log.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# ##########ALTERNATIVE COLOR SCALE LINEAR
# #Define quantiles
# b<-seq(.025,.5,by=.025)
# quants<-c(rev(b*-1),0,b,Inf)
# 
# gray.colors<-
#   gray.colors(length(quants/2), start = 0.1, end = 0.9, gamma = 2.2)
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c(gray.colors, jet.colors,"#7F0000")
# #quants<-c(-1,seq(0,.5,by=.025),Inf) #to specify quant cutoffs
# #specify colors
# #n<-length(quants)
# #cols<-rainbow(n,start=.4,end=0, alpha=1)
# #cols<-colorRampPalette(c("yellow", "red"),
# #               space = "Lab")
# #jet.colors <-
# #  colorRampPalette(c("blue", "#007FFF",
# #                     "#7FFF7F", "yellow", "#FF7F00", "red"))(length(quants))
# #cols<-c("grey",jet.colors,"#7F0000")
# 
# 
# txt0<-list("sp.text", c(20,-40), "POP_AVG")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_AVG2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3lin <- spplot(grd1, 'POP_AVG3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_AVG4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# ### when ready to export:
# p5 <- c(p4,p3lin,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POPAVG_linear.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ################POP_MAIZE LINEAR
# txt0<-list("sp.text", c(20,-40), "POP_MAIZE")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# ### when ready to export:
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_MAIZE_linear.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# ################POP_RICE LINEAR
# txt0<-list("sp.text", c(20,-40), "POP_RICE")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# ### when ready to export:
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_RICE_linear.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ################POP_WHEAT LINEAR
# txt0<-list("sp.text", c(20,-40), "POP_WHEAT")
# txt1<-list("sp.text", c(20,-50), "Phi = 0.01")
# p2 <- spplot(grd1, 'POP_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.001")
# p3 <- spplot(grd1, 'POP_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Phi = 0.0001")
# p4 <- spplot(grd1, 'POP_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# 
# ### when ready to export:
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/POP_WHEAT_linear.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# ############################# SCRIPT FOR MAKING A CUSTOM COLOR BAR (STILL NEEDS WORK, FIX OR SKIP)
# # this is still not correct.  The ticks need to be LOG SCALE.
# #ticks=quants
# #et.colors2 <-
# #  colorRampPalette(c("blue", "#007FFF",
# #                     "#7FFF7F", "yellow", "#FF7F00", "red"))(20)
# #gray.colors2<-gray.colors(20, start = 0.1, end = 0.9, gamma = 2.2)
# #col2<-c(gray.colors2,jet.colors2,"#7F0000","#7F0000")
# #color.bar(col2,-1.02,1.075, ticks=ticks)
# #note color bar is a custom function (saved in functions folder)
# 
# 
# 
# 
# quants<-c(-1,seq(0,500,by=10),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c("grey",jet.colors,"#7F0000")
# txt0<-list("sp.text", c(20,-40), "crop loss (tonnes / cell)")
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p2 <- spplot(grd1, 'CL2050_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p3 <- spplot(grd1, 'CL2050_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p4 <- spplot(grd1, 'CL2050_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Crop_loss_per_cell_MAIZE.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# txt0<-list("sp.text", c(20,-40), "crop loss (tonnes / cell)")
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p2 <- spplot(grd1, 'CL2050_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p3 <- spplot(grd1, 'CL2050_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p4 <- spplot(grd1, 'CL2050_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Crop_loss_per_cell_Rice.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# txt0<-list("sp.text", c(20,-40), "crop loss (tonnes / cell)")
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p2 <- spplot(grd1, 'CL2050_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p3 <- spplot(grd1, 'CL2050_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p4 <- spplot(grd1, 'CL2050_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Crop_loss_per_cell_Wheat.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# #### Crop loss due to climate change (proportion), IYCC
# quants<-c(-1,seq(0,.1,by=0.0001),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c("grey",jet.colors,"#7F0000")
# 
# txt0<-list("sp.text", c(20,-40), "% current crop lost due to CC (IYCC)")
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p2 <- spplot(grd1, 'IYCC_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p3 <- spplot(grd1, 'IYCC_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p4 <- spplot(grd1, 'IYCC_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss_prop_dueto_CC_MAIZE.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# txt0<-list("sp.text", c(20,-40), "% current crop lost due to CC (IYCC)")
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p2 <- spplot(grd1, 'IYCC_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p3 <- spplot(grd1, 'IYCC_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p4 <- spplot(grd1, 'IYCC_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss_prop_dueto_CC_Rice.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# txt0<-list("sp.text", c(20,-40), "% current crop lost due to CC (IYCC)")
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p2 <- spplot(grd1, 'IYCC_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p3 <- spplot(grd1, 'IYCC_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p4 <- spplot(grd1, 'IYCC_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss_prop_dueto_CC_Wheat.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# #### Current crop loss due to pests#
# quants<-c(seq(0,.2,by=0.01),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c(jet.colors,"#7F0000")
# 
# #txt0<-list("sp.text", c(20,-40), "% current crop loss")
# #txt1<-list("sp.text", c(20,-50), "Maize")
# p2 <- spplot(grd1, 'CLF_M', at=quants, sp.layout=list(wrld), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# #txt1<-list("sp.text", c(20,-50), "Rice")
# p3 <- spplot(grd1, 'CLF_R', at=quants, sp.layout=list(wrld), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# #txt1<-list("sp.text", c(20,-50), "Wheat")
# p4 <- spplot(grd1, 'CLF_W', at=quants, sp.layout=list(wrld), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Current crop loss.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# #### Future crop loss due to pests#
# txt0<-list("sp.text", c(20,-40), "crop loss 2dC increase")
# #Change to 4DC when rerun.
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p2 <- spplot(grd1, 'CLP2050_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p3 <- spplot(grd1, 'CLP2050_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p4 <- spplot(grd1, 'CLP2050_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss2dC_Maize.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# txt0<-list("sp.text", c(20,-40), "crop loss 2dC increase")
# #Change to 4DC when rerun w4dC model
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p2 <- spplot(grd1, 'CLP2050_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p3 <- spplot(grd1, 'CLP2050_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p4 <- spplot(grd1, 'CLP2050_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss2dC_Rice.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# txt0<-list("sp.text", c(20,-40), "crop loss 2dC increase")
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p2 <- spplot(grd1, 'CLP2050_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p3 <- spplot(grd1, 'CLP2050_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p4 <- spplot(grd1, 'CLP2050_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/croploss2dC_Wheat.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ###Crop growing season
# 
# quants<-c(-1,seq(0,220,by=10),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c("grey",jet.colors,"#7F0000")
# 
# txt0<-list("sp.text", c(20,-40), "crop growing season length")
# txt1<-list("sp.text", c(20,-50), "Maize")
# p2 <- spplot(grd1, 'CGS_M', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice")
# p3 <- spplot(grd1, 'CGS_R', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat")
# p4 <- spplot(grd1, 'CGS_W', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/growing season.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ###Current Yield
# 
# quants<-c(-1,seq(0,10,by=.1),Inf) #THIS IS A LINEAR QUANT SCALE
# 
# jet.colors <-
#   colorRampPalette(c("blue", "#007FFF",
#                      "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# cols<-c("grey",jet.colors,"#7F0000")
# 
# txt0<-list("sp.text", c(20,-40), "yield (tonnes / ha")
# txt1<-list("sp.text", c(20,-50), "Maize")
# p2 <- spplot(grd1, 'CYH_M', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice")
# p3 <- spplot(grd1, 'CYH_R', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat")
# p4 <- spplot(grd1, 'CYH_W', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Current_yield_per_ha.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# ###FUTURE YIELD
# 
# txt0<-list("sp.text", c(20,-40), "predicted yield (tonnes / ha")
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.01")
# p2 <- spplot(grd1, 'FCY_M2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.001")
# p3 <- spplot(grd1, 'FCY_M3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Maize_phi=0.0001")
# p4 <- spplot(grd1, 'FCY_M4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Future_yield_per_HA_MAIZE.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# txt0<-list("sp.text", c(20,-40), "predicted yield (tonnes / ha")
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.01")
# p2 <- spplot(grd1, 'FCY_R2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.001")
# p3 <- spplot(grd1, 'FCY_R3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Rice_phi=0.0001")
# p4 <- spplot(grd1, 'FCY_R4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Future_yield_per_HA_RICE.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# txt0<-list("sp.text", c(20,-40), "predicted yield (tonnes / ha")
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.01")
# p2 <- spplot(grd1, 'FCY_W2', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.001")
# p3 <- spplot(grd1, 'FCY_W3', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=FALSE)
# txt1<-list("sp.text", c(20,-50), "Wheat_phi=0.0001")
# p4 <- spplot(grd1, 'FCY_W4', at=quants, sp.layout=list(wrld,txt0,txt1), col.regions=cols, cuts=length(cols)-1,colorkey=list(
#   space='bottom', width=1, height=0.9, tick.number=9, labels = list(cex = 1.5)))
# 
# p5 <- c(p4,p3,p2, layout=c(1,3))
# png(file="~/Dropbox/climate change/pest MS shared/Future_yield_per_HA_WHEAT.png", width=900, height=1400) # CHECK RESOLUTION 
# print(p5)
# dev.off()
# 
# 
# # summary tables ----------------------------------------------------------
# 
# 
# 
# setwd("~/Dropbox/climate change/food security/climate and crop pressure MS/data/ascii_crops_hires")
# 
# C_TO_R<-data.frame(read.csv("COUNTRY_TO_REGION.csv", header=TRUE))
# 
# # now just need to use plyr and ddply to aggregate by region and subregion!
# setwd("~/Dropbox/climate change/pest MS shared")
# 
# MED_c<-ddply(ALL,.(NAME),summarise,
#              IPM_M2=median(IPM_M2,na.rm=TRUE),
#              IPM_M3=median(IPM_M3,na.rm=TRUE),
#              IPM_M4=median(IPM_M4,na.rm=TRUE),  
#              
#              IPM_R2=median(IPM_R2,na.rm=TRUE),
#              IPM_R3=median(IPM_R3,na.rm=TRUE),
#              IPM_R4=median(IPM_R4,na.rm=TRUE),
#              
#              IPM_W2=median(IPM_W2,na.rm=TRUE),
#              IPM_W3=median(IPM_W3,na.rm=TRUE),
#              IPM_W4=median(IPM_W4,na.rm=TRUE),
#              
#              IY_M2=median(IY_M2,na.rm=TRUE),
#              IY_M3=median(IY_M3,na.rm=TRUE),
#              IY_M4=median(IY_M4,na.rm=TRUE),  
#              
#              IY_R2=median(IY_R2,na.rm=TRUE),
#              IY_R3=median(IY_R3,na.rm=TRUE),
#              IY_R4=median(IY_R4,na.rm=TRUE),
#              
#              IY_W2=median(IY_W2,na.rm=TRUE),
#              IY_W3=median(IY_W3,na.rm=TRUE),
#              IY_W4=median(IY_W4,na.rm=TRUE),
#              
#              CLF_M=median(CLF_M,na.rm=TRUE),
#              CL_M=median(CL_M,na.rm=TRUE),
#              CY_M=median(CY_M,na.rm=TRUE),
#              CA_M=median(CA_M,na.rm=TRUE),
#              CGS_M=median(CGS_M,na.rm=TRUE),
#              
#              CLF_R=median(CLF_R,na.rm=TRUE),
#              CL_R=median(CL_R,na.rm=TRUE),
#              CY_R=median(CY_R,na.rm=TRUE),
#              CA_R=median(CA_R,na.rm=TRUE),
#              CGS_R=median(CGS_R,na.rm=TRUE),
#              
#              CLF_W=median(CLF_W,na.rm=TRUE),
#              CL_W=median(CL_W,na.rm=TRUE),
#              CY_W=median(CY_W,na.rm=TRUE),
#              CA_W=median(CA_W,na.rm=TRUE),
#              CGS_W=median(CGS_W,na.rm=TRUE),
#              
#              MET_M2=median(MET_M2,na.rm=TRUE),
#              MET_M3=median(MET_M3,na.rm=TRUE),
#              MET_M4=median(MET_M4,na.rm=TRUE),  
#              
#              MET_R2=median(MET_R2,na.rm=TRUE),
#              MET_R3=median(MET_R3,na.rm=TRUE),
#              MET_R4=median(MET_R4,na.rm=TRUE),
#              
#              MET_W2=median(MET_W2,na.rm=TRUE),
#              MET_W3=median(MET_W3,na.rm=TRUE),
#              MET_W4=median(MET_W4,na.rm=TRUE),
#              
#              POP_M2=median(POP_M2,na.rm=TRUE),
#              POP_M3=median(POP_M3,na.rm=TRUE),
#              POP_M4=median(POP_M4,na.rm=TRUE),  
#              
#              POP_R2=median(POP_R2,na.rm=TRUE),
#              POP_R3=median(POP_R3,na.rm=TRUE),
#              POP_R4=median(POP_R4,na.rm=TRUE),
#              
#              POP_W2=median(POP_W2,na.rm=TRUE),
#              POP_W3=median(POP_W3,na.rm=TRUE),
#              POP_W4=median(POP_W4,na.rm=TRUE),
#              
#              MET_AVG2=median(MET_AVG2,na.rm=TRUE),
#              MET_AVG3=median(MET_AVG3,na.rm=TRUE),
#              MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
#              
#              POP_AVG2=median(POP_AVG2,na.rm=TRUE),
#              POP_AVG3=median(POP_AVG3,na.rm=TRUE),
#              POP_AVG4=median(POP_AVG4,na.rm=TRUE),
#              
#              IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
#              IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
#              IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
#              
#              YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#              YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#              YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#              YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
#              YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
#              YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
#              
#              
#              CL2050_M2=median(CL2050_M2,na.rm=TRUE),
#              CL2050_M3=median(CL2050_M3,na.rm=TRUE),
#              CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
#              
#              CL2050_R2=median(CL2050_R2,na.rm=TRUE),
#              CL2050_R3=median(CL2050_R3,na.rm=TRUE),
#              CL2050_R4=median(CL2050_R4,na.rm=TRUE),
#              
#              CL2050_W2=median(CL2050_W2,na.rm=TRUE),
#              CL2050_W3=median(CL2050_W3,na.rm=TRUE),
#              CL2050_W4=median(CL2050_W4,na.rm=TRUE),
#              
#              CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
#              CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
#              CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
#              
#              CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
#              CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
#              CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
#              
#              CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
#              CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
#              CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
#              
#              IYCC_M2=median(IYCC_M2,na.rm=TRUE),
#              IYCC_M3=median(IYCC_M3,na.rm=TRUE),
#              IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
#              
#              IYCC_R2=median(IYCC_R2,na.rm=TRUE),
#              IYCC_R3=median(IYCC_R3,na.rm=TRUE),
#              IYCC_R4=median(IYCC_R4,na.rm=TRUE),
#              
#              IYCC_W2=median(IYCC_W2,na.rm=TRUE),
#              IYCC_W3=median(IYCC_W3,na.rm=TRUE),
#              IYCC_W4=median(IYCC_W4,na.rm=TRUE),
#              
#              YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#              YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#              YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#              
#              YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#              YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#              YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#              
#              YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#              YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#              YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))
# 
# 
# #MEDIEANS BY SUBREGION
# 
# MED_sr<-ddply(ALL,.(Subregion),summarise,
#               IPM_M2=median(IPM_M2,na.rm=TRUE),
#               IPM_M3=median(IPM_M3,na.rm=TRUE),
#               IPM_M4=median(IPM_M4,na.rm=TRUE),  
#               
#               IPM_R2=median(IPM_R2,na.rm=TRUE),
#               IPM_R3=median(IPM_R3,na.rm=TRUE),
#               IPM_R4=median(IPM_R4,na.rm=TRUE),
#               
#               IPM_W2=median(IPM_W2,na.rm=TRUE),
#               IPM_W3=median(IPM_W3,na.rm=TRUE),
#               IPM_W4=median(IPM_W4,na.rm=TRUE),
#               
#               IY_M2=median(IY_M2,na.rm=TRUE),
#               IY_M3=median(IY_M3,na.rm=TRUE),
#               IY_M4=median(IY_M4,na.rm=TRUE),  
#               
#               IY_R2=median(IY_R2,na.rm=TRUE),
#               IY_R3=median(IY_R3,na.rm=TRUE),
#               IY_R4=median(IY_R4,na.rm=TRUE),
#               
#               IY_W2=median(IY_W2,na.rm=TRUE),
#               IY_W3=median(IY_W3,na.rm=TRUE),
#               IY_W4=median(IY_W4,na.rm=TRUE),
#               
#               CLF_M=median(CLF_M,na.rm=TRUE),
#               CL_M=median(CL_M,na.rm=TRUE),
#               CY_M=median(CY_M,na.rm=TRUE),
#               CA_M=median(CA_M,na.rm=TRUE),
#               CGS_M=median(CGS_M,na.rm=TRUE),
#               
#               CLF_R=median(CLF_R,na.rm=TRUE),
#               CL_R=median(CL_R,na.rm=TRUE),
#               CY_R=median(CY_R,na.rm=TRUE),
#               CA_R=median(CA_R,na.rm=TRUE),
#               CGS_R=median(CGS_R,na.rm=TRUE),
#               
#               CLF_W=median(CLF_W,na.rm=TRUE),
#               CL_W=median(CL_W,na.rm=TRUE),
#               CY_W=median(CY_W,na.rm=TRUE),
#               CA_W=median(CA_W,na.rm=TRUE),
#               CGS_W=median(CGS_W,na.rm=TRUE),
#               
#               MET_M2=median(MET_M2,na.rm=TRUE),
#               MET_M3=median(MET_M3,na.rm=TRUE),
#               MET_M4=median(MET_M4,na.rm=TRUE),  
#               
#               MET_R2=median(MET_R2,na.rm=TRUE),
#               MET_R3=median(MET_R3,na.rm=TRUE),
#               MET_R4=median(MET_R4,na.rm=TRUE),
#               
#               MET_W2=median(MET_W2,na.rm=TRUE),
#               MET_W3=median(MET_W3,na.rm=TRUE),
#               MET_W4=median(MET_W4,na.rm=TRUE),
#               
#               POP_M2=median(POP_M2,na.rm=TRUE),
#               POP_M3=median(POP_M3,na.rm=TRUE),
#               POP_M4=median(POP_M4,na.rm=TRUE),  
#               
#               POP_R2=median(POP_R2,na.rm=TRUE),
#               POP_R3=median(POP_R3,na.rm=TRUE),
#               POP_R4=median(POP_R4,na.rm=TRUE),
#               
#               POP_W2=median(POP_W2,na.rm=TRUE),
#               POP_W3=median(POP_W3,na.rm=TRUE),
#               POP_W4=median(POP_W4,na.rm=TRUE),
#               
#               MET_AVG2=median(MET_AVG2,na.rm=TRUE),
#               MET_AVG3=median(MET_AVG3,na.rm=TRUE),
#               MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
#               
#               POP_AVG2=median(POP_AVG2,na.rm=TRUE),
#               POP_AVG3=median(POP_AVG3,na.rm=TRUE),
#               POP_AVG4=median(POP_AVG4,na.rm=TRUE),
#               
#               IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
#               IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
#               IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
#               
#               YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#               YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#               YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#               YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
#               YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
#               YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
#               
#               
#               CL2050_M2=median(CL2050_M2,na.rm=TRUE),
#               CL2050_M3=median(CL2050_M3,na.rm=TRUE),
#               CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
#               
#               CL2050_R2=median(CL2050_R2,na.rm=TRUE),
#               CL2050_R3=median(CL2050_R3,na.rm=TRUE),
#               CL2050_R4=median(CL2050_R4,na.rm=TRUE),
#               
#               CL2050_W2=median(CL2050_W2,na.rm=TRUE),
#               CL2050_W3=median(CL2050_W3,na.rm=TRUE),
#               CL2050_W4=median(CL2050_W4,na.rm=TRUE),
#               
#               CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
#               CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
#               CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
#               
#               CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
#               CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
#               CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
#               
#               CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
#               CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
#               CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
#               
#               IYCC_M2=median(IYCC_M2,na.rm=TRUE),
#               IYCC_M3=median(IYCC_M3,na.rm=TRUE),
#               IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
#               
#               IYCC_R2=median(IYCC_R2,na.rm=TRUE),
#               IYCC_R3=median(IYCC_R3,na.rm=TRUE),
#               IYCC_R4=median(IYCC_R4,na.rm=TRUE),
#               
#               IYCC_W2=median(IYCC_W2,na.rm=TRUE),
#               IYCC_W3=median(IYCC_W3,na.rm=TRUE),
#               IYCC_W4=median(IYCC_W4,na.rm=TRUE),
# 
#               YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#               YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#               YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#               
#               YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#               YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               
# YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
# YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
# YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))
# 
# 
# #MEDIANS by region
# 
# MED_r<-ddply(ALL,.(Region),summarise,
#              IPM_M2=median(IPM_M2,na.rm=TRUE),
#              IPM_M3=median(IPM_M3,na.rm=TRUE),
#              IPM_M4=median(IPM_M4,na.rm=TRUE),  
#              
#              IPM_R2=median(IPM_R2,na.rm=TRUE),
#              IPM_R3=median(IPM_R3,na.rm=TRUE),
#              IPM_R4=median(IPM_R4,na.rm=TRUE),
#              
#              IPM_W2=median(IPM_W2,na.rm=TRUE),
#              IPM_W3=median(IPM_W3,na.rm=TRUE),
#              IPM_W4=median(IPM_W4,na.rm=TRUE),
#              
#              IY_M2=median(IY_M2,na.rm=TRUE),
#              IY_M3=median(IY_M3,na.rm=TRUE),
#              IY_M4=median(IY_M4,na.rm=TRUE),  
#              
#              IY_R2=median(IY_R2,na.rm=TRUE),
#              IY_R3=median(IY_R3,na.rm=TRUE),
#              IY_R4=median(IY_R4,na.rm=TRUE),
#              
#              IY_W2=median(IY_W2,na.rm=TRUE),
#              IY_W3=median(IY_W3,na.rm=TRUE),
#              IY_W4=median(IY_W4,na.rm=TRUE),
#              
#              CLF_M=median(CLF_M,na.rm=TRUE),
#              CL_M=median(CL_M,na.rm=TRUE),
#              CY_M=median(CY_M,na.rm=TRUE),
#              CA_M=median(CA_M,na.rm=TRUE),
#              CGS_M=median(CGS_M,na.rm=TRUE),
#              
#              CLF_R=median(CLF_R,na.rm=TRUE),
#              CL_R=median(CL_R,na.rm=TRUE),
#              CY_R=median(CY_R,na.rm=TRUE),
#              CA_R=median(CA_R,na.rm=TRUE),
#              CGS_R=median(CGS_R,na.rm=TRUE),
#              
#              CLF_W=median(CLF_W,na.rm=TRUE),
#              CL_W=median(CL_W,na.rm=TRUE),
#              CY_W=median(CY_W,na.rm=TRUE),
#              CA_W=median(CA_W,na.rm=TRUE),
#              CGS_W=median(CGS_W,na.rm=TRUE),
#              
#              MET_M2=median(MET_M2,na.rm=TRUE),
#              MET_M3=median(MET_M3,na.rm=TRUE),
#              MET_M4=median(MET_M4,na.rm=TRUE),  
#              
#              MET_R2=median(MET_R2,na.rm=TRUE),
#              MET_R3=median(MET_R3,na.rm=TRUE),
#              MET_R4=median(MET_R4,na.rm=TRUE),
#              
#              MET_W2=median(MET_W2,na.rm=TRUE),
#              MET_W3=median(MET_W3,na.rm=TRUE),
#              MET_W4=median(MET_W4,na.rm=TRUE),
#              
#              POP_M2=median(POP_M2,na.rm=TRUE),
#              POP_M3=median(POP_M3,na.rm=TRUE),
#              POP_M4=median(POP_M4,na.rm=TRUE),  
#              
#              POP_R2=median(POP_R2,na.rm=TRUE),
#              POP_R3=median(POP_R3,na.rm=TRUE),
#              POP_R4=median(POP_R4,na.rm=TRUE),
#              
#              POP_W2=median(POP_W2,na.rm=TRUE),
#              POP_W3=median(POP_W3,na.rm=TRUE),
#              POP_W4=median(POP_W4,na.rm=TRUE),
#              
#              MET_AVG2=median(MET_AVG2,na.rm=TRUE),
#              MET_AVG3=median(MET_AVG3,na.rm=TRUE),
#              MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
#              
#              POP_AVG2=median(POP_AVG2,na.rm=TRUE),
#              POP_AVG3=median(POP_AVG3,na.rm=TRUE),
#              POP_AVG4=median(POP_AVG4,na.rm=TRUE),
#              
#              IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
#              IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
#              IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
#              
#              YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#              YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#              YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#              YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
#              YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
#              YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
#              
#              
#              CL2050_M2=median(CL2050_M2,na.rm=TRUE),
#              CL2050_M3=median(CL2050_M3,na.rm=TRUE),
#              CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
#              
#              CL2050_R2=median(CL2050_R2,na.rm=TRUE),
#              CL2050_R3=median(CL2050_R3,na.rm=TRUE),
#              CL2050_R4=median(CL2050_R4,na.rm=TRUE),
#              
#              CL2050_W2=median(CL2050_W2,na.rm=TRUE),
#              CL2050_W3=median(CL2050_W3,na.rm=TRUE),
#              CL2050_W4=median(CL2050_W4,na.rm=TRUE),
#              
#              CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
#              CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
#              CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
#              
#              CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
#              CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
#              CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
#              
#              CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
#              CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
#              CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
#              
#              IYCC_M2=median(IYCC_M2,na.rm=TRUE),
#              IYCC_M3=median(IYCC_M3,na.rm=TRUE),
#              IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
#              
#              IYCC_R2=median(IYCC_R2,na.rm=TRUE),
#              IYCC_R3=median(IYCC_R3,na.rm=TRUE),
#              IYCC_R4=median(IYCC_R4,na.rm=TRUE),
#              
#              IYCC_W2=median(IYCC_W2,na.rm=TRUE),
#              IYCC_W3=median(IYCC_W3,na.rm=TRUE),
#              IYCC_W4=median(IYCC_W4,na.rm=TRUE),
# 
# YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
# YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
# YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
# 
# YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
# YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
# YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
# 
# YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
# YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
# YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))
# 
# CL_c<-ddply(ALL,.(NAME),summarise,                       
#             MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#             MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#             MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#             TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
#             TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
#             TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
#             TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
#             TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
#             TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
#             TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
#             TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
#             TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))
# 
# 
# #MEDIEANS BY SUBREGION
# 
# CL_sr<-ddply(ALL,.(Subregion),summarise,
#              MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#              MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#              MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#              TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
#              TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
#              TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
#              TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
#              TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
#              TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
#              TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
#              TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
#              TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
#              TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
#              TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
#              TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))
# 
# #MEDIANS by region
# 
# CL_r<-ddply(ALL,.(Region),summarise,
#             MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#             MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#             MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#             TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
#             TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
#             TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
#             TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
#             TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
#             TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
#             TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
#             TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
#             TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
#             TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))
# 
# MED_r<-merge(MED_r,CL_r,by="Region",all=TRUE,sort=TRUE)
# MED_sr<-merge(MED_sr,CL_sr,by="Subregion",all=TRUE,sort=TRUE)
# MED_c<-merge(MED_c,CL_c,by="NAME",all=TRUE,sort=TRUE)
# 
# 
# # country tables ----------------------------------------------------------
# 
# 
# #making a datatables for IPM and yield impact by country
# 
# 
# 
# MAIZE<-ALL[complete.cases(ALL[ ,7:39]),] # subsets ALL to give only maize data
# #Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
# TOT_YLD_MAIZE<-ddply(MAIZE,.(NAME),summarise,                       
#                TOT_YLD_TOT_M=sum(YLD_TOT_M))
# 
# MAIZE<-merge(MAIZE,TOT_YLD_MAIZE, by=c("NAME"),all=TRUE,sort=TRUE)
# MAIZE_c<-ddply(MAIZE,.(NAME),summarise,                       
#                MEAN_LAT_M=mean(LAT),
#                WMEAN_LAT_M=sum((LAT*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MED_YLD_HA_M=median(YLD_HA_M),
#                MEAN_YLD_HA_M=mean(YLD_HA_M),
#                WMEAN_YLD_HA_M=sum((YLD_HA_M*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_YLD_HA_M=min(YLD_HA_M),
#                MAX_YLD_HA_M=max(YLD_HA_M),
#                
#                TOT_YLD_TOT_M=sum(YLD_TOT_M),
#                
#                TOT_CL2050_M2=sum(CL2050_M2), 
#                TOT_CL2050_M3=sum(CL2050_M3),
#                TOT_CL2050_M4=sum(CL2050_M4),  
#                
#                MED_IPM_M2=median(IPM_M2),
#                MEAN_IPM_M2=mean(IPM_M2),
#                WMEAN_IPM_M2=sum((IPM_M2*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_IPM_M2=min(IPM_M2),
#                MAX_IPM_M2=max(IPM_M2),
#                
#                MED_IPM_M3=median(IPM_M3),
#                MEAN_IPM_M3=mean(IPM_M3),
#                WMEAN_IPM_M3=sum((IPM_M3*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_IPM_M3=min(IPM_M3),
#                MAX_IPM_M3=max(IPM_M3),
#                
#                MED_IPM_M4=median(IPM_M4),
#                MEAN_IPM_M4=mean(IPM_M4),
#                WMEAN_IPM_M4=sum((IPM_M4*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_IPM_M4=min(IPM_M4),
#                MAX_IPM_M4=max(IPM_M4),
#                
#                MED_IYCC_M2=median(IYCC_M2),
#                MEAN_IYCC_M2=mean(IYCC_M2),
#                MIN_IYCC_M2=min(IYCC_M2),
#                MAX_IYCC_M2=max(IYCC_M2),
#                
#                MED_IYCC_M3=median(IYCC_M3),
#                MEAN_IYCC_M3=mean(IYCC_M3),
#                MIN_IYCC_M3=min(IYCC_M3),
#                MAX_IYCC_M3=max(IYCC_M3),
#                
#                MED_IYCC_M4=median(IYCC_M4),
#                MEAN_IYCC_M4=mean(IYCC_M4),
#                MIN_IYCC_M4=min(IYCC_M4),
#                MAX_IYCC_M4=max(IYCC_M4),   
#                
#                MED_YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M)), 
#                MEAN_YLPH_M2=mean(IPM_M2*CLF_M*(CY_M/CA_M)),
#                WMEAN_YLPH_M2=sum(((IPM_M2*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_YLPH_M2=min(IPM_M2*CLF_M*(CY_M/CA_M)),
#                MAX_YLPH_M2=max(IPM_M2*CLF_M*(CY_M/CA_M)),
#                
#                MED_YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M)),
#                MEAN_YLPH_M3=mean(IPM_M3*CLF_M*(CY_M/CA_M)),
#                WMEAN_YLPH_M3=sum(((IPM_M3*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_YLPH_M3=min(IPM_M3*CLF_M*(CY_M/CA_M)),
#                MAX_YLPH_M3=max(IPM_M3*CLF_M*(CY_M/CA_M)),
#                
#                MED_YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M)),
#                MEAN_YLPH_M4=mean(IPM_M4*CLF_M*(CY_M/CA_M)),
#                WMEAN_YLPH_M4=sum(((IPM_M4*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
#                MIN_YLPH_M4=min(IPM_M4*CLF_M*(CY_M/CA_M)),
#                MAX_YLPH_M4=max(IPM_M4*CLF_M*(CY_M/CA_M)),
#       
#                CELLS_M=sum(CELLS))   
# 
# 
# RICE<-ALL[complete.cases(ALL[ ,40:72]),] # subsets ALL to give only maize data
# #Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
# TOT_YLD_RICE<-ddply(RICE,.(NAME),summarise,                       
#                      TOT_YLD_TOT_R=sum(YLD_TOT_R))
# 
# RICE<-merge(RICE,TOT_YLD_RICE, by=c("NAME"),all=TRUE,sort=TRUE)
# 
# 
# RICE_c<-ddply(RICE,.(NAME),summarise,                       
#                MEAN_LAT_R=mean(LAT),
#                WMEAN_LAT_R=sum((LAT*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MED_YLD_HA_R=median(YLD_HA_R),
#                MEAN_YLD_HA_R=mean(YLD_HA_R),
#                WMEAN_YLD_HA_R=sum((YLD_HA_R*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_YLD_HA_R=min(YLD_HA_R),
#                MAX_YLD_HA_R=max(YLD_HA_R),
#                
#                TOT_YLD_TOT_R=sum(YLD_TOT_R),
#                
#                TOT_CL2050_R2=sum(CL2050_R2), 
#                TOT_CL2050_R3=sum(CL2050_R3),
#                TOT_CL2050_R4=sum(CL2050_R4),  
#                
#                MED_IPM_R2=median(IPM_R2),
#                MEAN_IPM_R2=mean(IPM_R2),
#                WMEAN_IPM_R2=sum((IPM_R2*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_IPM_R2=min(IPM_R2),
#                MAX_IPM_R2=max(IPM_R2),
#                
#                MED_IPM_R3=median(IPM_R3),
#                MEAN_IPM_R3=mean(IPM_R3),
#                WMEAN_IPM_R3=sum((IPM_R3*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_IPM_R3=min(IPM_R3),
#                MAX_IPM_R3=max(IPM_R3),
#                
#                MED_IPM_R4=median(IPM_R4),
#                MEAN_IPM_R4=mean(IPM_R4),
#                WMEAN_IPM_R4=sum((IPM_R4*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_IPM_R4=min(IPM_R4),
#                MAX_IPM_R4=max(IPM_R4),
#                
#                MED_IYCC_R2=median(IYCC_R2),
#                MEAN_IYCC_R2=mean(IYCC_R2),
#                MIN_IYCC_R2=min(IYCC_R2),
#                MAX_IYCC_R2=max(IYCC_R2),
#                
#                MED_IYCC_R3=median(IYCC_R3),
#                MEAN_IYCC_R3=mean(IYCC_R3),
#                MIN_IYCC_R3=min(IYCC_R3),
#                MAX_IYCC_R3=max(IYCC_R3),
#                
#                MED_IYCC_R4=median(IYCC_R4),
#                MEAN_IYCC_R4=mean(IYCC_R4),
#                MIN_IYCC_R4=min(IYCC_R4),
#                MAX_IYCC_R4=max(IYCC_R4),   
#                
#                MED_YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R)), 
#                MEAN_YLPH_R2=mean(IPM_R2*CLF_R*(CY_R/CA_R)),
#                WMEAN_YLPH_R2=sum(((IPM_R2*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_YLPH_R2=min(IPM_R2*CLF_R*(CY_R/CA_R)),
#                MAX_YLPH_R2=max(IPM_R2*CLF_R*(CY_R/CA_R)),
#                
#                MED_YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R)),
#                MEAN_YLPH_R3=mean(IPM_R3*CLF_R*(CY_R/CA_R)),
#                WMEAN_YLPH_R3=sum(((IPM_R3*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_YLPH_R3=min(IPM_R3*CLF_R*(CY_R/CA_R)),
#                MAX_YLPH_R3=max(IPM_R3*CLF_R*(CY_R/CA_R)),
#                
#                MED_YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R)),
#                MEAN_YLPH_R4=mean(IPM_R4*CLF_R*(CY_R/CA_R)),
#                WMEAN_YLPH_R4=sum(((IPM_R4*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
#                MIN_YLPH_R4=min(IPM_R4*CLF_R*(CY_R/CA_R)),
#                MAX_YLPH_R4=max(IPM_R4*CLF_R*(CY_R/CA_R)),
#                
#                CELLS_R=sum(CELLS))
# 
# WHEAT<-ALL[complete.cases(ALL[ ,73:105]),] # subsets ALL to give only maize data
# 
# #Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
# #Not sure what the last line means.... hmm
# 
# TOT_YLD_WHEAT<-ddply(WHEAT,.(NAME),summarise,                       
#                     TOT_YLD_TOT_W=sum(YLD_TOT_W))
# 
# WHEAT<-merge(WHEAT,TOT_YLD_WHEAT, by=c("NAME"),all=TRUE,sort=TRUE)
# 
# 
# 
# WHEAT_c<-ddply(WHEAT,.(NAME),summarise,                       
#               MEAN_LAT_W=mean(LAT),
#               WMEAN_LAT_W=sum((LAT*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MED_YLD_HA_W=median(YLD_HA_W),
#               MEAN_YLD_HA_W=mean(YLD_HA_W),
#               WMEAN_YLD_HA_W=sum((YLD_HA_W*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_YLD_HA_W=min(YLD_HA_W),
#               MAX_YLD_HA_W=max(YLD_HA_W),
#               
#               TOT_YLD_TOT_W=sum(YLD_TOT_W),
#               
#               TOT_CL2050_W2=sum(CL2050_W2), 
#               TOT_CL2050_W3=sum(CL2050_W3),
#               TOT_CL2050_W4=sum(CL2050_W4),  
#               
#               MED_IPM_W2=median(IPM_W2),
#               MEAN_IPM_W2=mean(IPM_W2),
#               WMEAN_IPM_W2=sum((IPM_W2*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_IPM_W2=min(IPM_W2),
#               MAX_IPM_W2=max(IPM_W2),
#               
#               MED_IPM_W3=median(IPM_W3),
#               MEAN_IPM_W3=mean(IPM_W3),
#               WMEAN_IPM_W3=sum((IPM_W3*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_IPM_W3=min(IPM_W3),
#               MAX_IPM_W3=max(IPM_W3),
#               
#               MED_IPM_W4=median(IPM_W4),
#               MEAN_IPM_W4=mean(IPM_W4),
#               WMEAN_IPM_W4=sum((IPM_W4*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_IPM_W4=min(IPM_W4),
#               MAX_IPM_W4=max(IPM_W4),
#               
#               MED_IYCC_W2=median(IYCC_W2),
#               MEAN_IYCC_W2=mean(IYCC_W2),
#               MIN_IYCC_W2=min(IYCC_W2),
#               MAX_IYCC_W2=max(IYCC_W2),
#               
#               MED_IYCC_W3=median(IYCC_W3),
#               MEAN_IYCC_W3=mean(IYCC_W3),
#               MIN_IYCC_W3=min(IYCC_W3),
#               MAX_IYCC_W3=max(IYCC_W3),
#               
#               MED_IYCC_W4=median(IYCC_W4),
#               MEAN_IYCC_W4=mean(IYCC_W4),
#               MIN_IYCC_W4=min(IYCC_W4),
#               MAX_IYCC_W4=max(IYCC_W4),   
#               
#               MED_YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W)), 
#               MEAN_YLPH_W2=mean(IPM_W2*CLF_W*(CY_W/CA_W)),
#               WMEAN_YLPH_W2=sum(((IPM_W2*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_YLPH_W2=min(IPM_W2*CLF_W*(CY_W/CA_W)),
#               MAX_YLPH_W2=max(IPM_W2*CLF_W*(CY_W/CA_W)),
#               
#               MED_YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W)),
#               MEAN_YLPH_W3=mean(IPM_W3*CLF_W*(CY_W/CA_W)),
#               WMEAN_YLPH_W3=sum(((IPM_W3*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_YLPH_W3=min(IPM_W3*CLF_W*(CY_W/CA_W)),
#               MAX_YLPH_W3=max(IPM_W3*CLF_W*(CY_W/CA_W)),
#               
#               MED_YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W)),
#               MEAN_YLPH_W4=mean(IPM_W4*CLF_W*(CY_W/CA_W)),
#               WMEAN_YLPH_W4=sum(((IPM_W4*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
#               MIN_YLPH_W4=min(IPM_W4*CLF_W*(CY_W/CA_W)),
#               MAX_YLPH_W4=max(IPM_W4*CLF_W*(CY_W/CA_W)),
#               
#               CELLS_W=sum(CELLS))
# 
# 
# ####BY SUBREGION
# 
# 
# 
# 
# MAIZE_r<-ddply(ALL,.(Region),summarise,                       
#                MED_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
#                MEAN_YLD_HA_M=mean(YLD_HA_M,na.rm=TRUE),
#                MIN_YLD_HA_M=min(YLD_HA_M,na.rm=TRUE),
#                MAX_YLD_HA_M=max(YLD_HA_M,na.rm=TRUE),
#                
#                TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
#                
#                TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
#                TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
#                TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
#                
#                MED_IPM_M2=median(IPM_M2,na.rm=TRUE),
#                MEAN_IPM_M2=mean(IPM_M2,na.rm=TRUE),
#                MIN_IPM_M2=min(IPM_M2,na.rm=TRUE),
#                MAX_IPM_M2=max(IPM_M2,na.rm=TRUE),
#                
#                MED_IPM_M3=median(IPM_M3,na.rm=TRUE),
#                MEAN_IPM_M3=mean(IPM_M3,na.rm=TRUE),
#                MIN_IPM_M3=min(IPM_M3,na.rm=TRUE),
#                MAX_IPM_M3=max(IPM_M3,na.rm=TRUE),
#                
#                MED_IPM_M4=median(IPM_M4,na.rm=TRUE),
#                MEAN_IPM_M4=mean(IPM_M4,na.rm=TRUE),
#                MIN_IPM_M4=min(IPM_M4,na.rm=TRUE),
#                MAX_IPM_M4=max(IPM_M4,na.rm=TRUE),
#                
#                MED_IYCC_M2=median(IYCC_M2,na.rm=TRUE),
#                MEAN_IYCC_M2=mean(IYCC_M2,na.rm=TRUE),
#                MIN_IYCC_M2=min(IYCC_M2,na.rm=TRUE),
#                MAX_IYCC_M2=max(IYCC_M2,na.rm=TRUE),
#                
#                MED_IYCC_M3=median(IYCC_M3,na.rm=TRUE),
#                MEAN_IYCC_M3=mean(IYCC_M3,na.rm=TRUE),
#                MIN_IYCC_M3=min(IYCC_M3,na.rm=TRUE),
#                MAX_IYCC_M3=max(IYCC_M3,na.rm=TRUE),
#                
#                MED_IYCC_M4=median(IYCC_M4,na.rm=TRUE),
#                MEAN_IYCC_M4=mean(IYCC_M4,na.rm=TRUE),
#                MIN_IYCC_M4=min(IYCC_M4,na.rm=TRUE),
#                MAX_IYCC_M4=max(IYCC_M4,na.rm=TRUE),   
#                
#                MED_YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#                MEAN_YLPH_M2=mean(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MIN_YLPH_M2=min(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MAX_YLPH_M2=max(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                
#                MED_YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MEAN_YLPH_M3=mean(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MIN_YLPH_M3=min(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MAX_YLPH_M3=max(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                
#                MED_YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MEAN_YLPH_M4=mean(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MIN_YLPH_M4=min(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                MAX_YLPH_M4=max(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
#                
#                CELLS_M=sum(CELLS,na.rm=TRUE))   
# 
# 
# RICE_r<-ddply(ALL,.(Region),summarise,                       
#               MED_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
#               MEAN_YLD_HA_R=mean(YLD_HA_R,na.rm=TRUE),
#               MIN_YLD_HA_R=min(YLD_HA_R,na.rm=TRUE),
#               MAX_YLD_HA_R=max(YLD_HA_R,na.rm=TRUE),
#               
#               TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
#               
#               TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
#               TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
#               TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),  
#               
#               MED_IPM_R2=median(IPM_R2,na.rm=TRUE),
#               MEAN_IPM_R2=mean(IPM_R2,na.rm=TRUE),
#               MIN_IPM_R2=min(IPM_R2,na.rm=TRUE),
#               MAX_IPM_R2=max(IPM_R2,na.rm=TRUE),
#               
#               MED_IPM_R3=median(IPM_R3,na.rm=TRUE),
#               MEAN_IPM_R3=mean(IPM_R3,na.rm=TRUE),
#               MIN_IPM_R3=min(IPM_R3,na.rm=TRUE),
#               MAX_IPM_R3=max(IPM_R3,na.rm=TRUE),
#               
#               MED_IPM_R4=median(IPM_R4,na.rm=TRUE),
#               MEAN_IPM_R4=mean(IPM_R4,na.rm=TRUE),
#               MIN_IPM_R4=min(IPM_R4,na.rm=TRUE),
#               MAX_IPM_R4=max(IPM_R4,na.rm=TRUE),
#               
#               MED_IYCC_R2=median(IYCC_R2,na.rm=TRUE),
#               MEAN_IYCC_R2=mean(IYCC_R2,na.rm=TRUE),
#               MIN_IYCC_R2=min(IYCC_R2,na.rm=TRUE),
#               MAX_IYCC_R2=max(IYCC_R2,na.rm=TRUE),
#               
#               MED_IYCC_R3=median(IYCC_R3,na.rm=TRUE),
#               MEAN_IYCC_R3=mean(IYCC_R3,na.rm=TRUE),
#               MIN_IYCC_R3=min(IYCC_R3,na.rm=TRUE),
#               MAX_IYCC_R3=max(IYCC_R3,na.rm=TRUE),
#               
#               MED_IYCC_R4=median(IYCC_R4,na.rm=TRUE),
#               MEAN_IYCC_R4=mean(IYCC_R4,na.rm=TRUE),
#               MIN_IYCC_R4=min(IYCC_R4,na.rm=TRUE),
#               MAX_IYCC_R4=max(IYCC_R4,na.rm=TRUE),   
#               
#               MED_YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#               MEAN_YLPH_R2=mean(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MIN_YLPH_R2=min(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MAX_YLPH_R2=max(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               
#               MED_YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MEAN_YLPH_R3=mean(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MIN_YLPH_R3=min(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MAX_YLPH_R3=max(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               
#               MED_YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MEAN_YLPH_R4=mean(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MIN_YLPH_R4=min(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               MAX_YLPH_R4=max(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
#               
#               CELLS_R=sum(CELLS,na.rm=TRUE))
# 
# WHEAT_r<-ddply(ALL,.(Region),summarise,                       
#                MED_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
#                MEAN_YLD_HA_W=mean(YLD_HA_W,na.rm=TRUE),
#                MIN_YLD_HA_W=min(YLD_HA_W,na.rm=TRUE),
#                MAX_YLD_HA_W=max(YLD_HA_W,na.rm=TRUE),
#                
#                TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
#                
#                TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
#                TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
#                TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE),  
#                
#                MED_IPM_W2=median(IPM_W2,na.rm=TRUE),
#                MEAN_IPM_W2=mean(IPM_W2,na.rm=TRUE),
#                MIN_IPM_W2=min(IPM_W2,na.rm=TRUE),
#                MAX_IPM_W2=max(IPM_W2,na.rm=TRUE),
#                
#                MED_IPM_W3=median(IPM_W3,na.rm=TRUE),
#                MEAN_IPM_W3=mean(IPM_W3,na.rm=TRUE),
#                MIN_IPM_W3=min(IPM_W3,na.rm=TRUE),
#                MAX_IPM_W3=max(IPM_W3,na.rm=TRUE),
#                
#                MED_IPM_W4=median(IPM_W4,na.rm=TRUE),
#                MEAN_IPM_W4=mean(IPM_W4,na.rm=TRUE),
#                MIN_IPM_W4=min(IPM_W4,na.rm=TRUE),
#                MAX_IPM_W4=max(IPM_W4,na.rm=TRUE),
#                
#                MED_IYCC_W2=median(IYCC_W2,na.rm=TRUE),
#                MEAN_IYCC_W2=mean(IYCC_W2,na.rm=TRUE),
#                MIN_IYCC_W2=min(IYCC_W2,na.rm=TRUE),
#                MAX_IYCC_W2=max(IYCC_W2,na.rm=TRUE),
#                
#                MED_IYCC_W3=median(IYCC_W3,na.rm=TRUE),
#                MEAN_IYCC_W3=mean(IYCC_W3,na.rm=TRUE),
#                MIN_IYCC_W3=min(IYCC_W3,na.rm=TRUE),
#                MAX_IYCC_W3=max(IYCC_W3,na.rm=TRUE),
#                
#                MED_IYCC_W4=median(IYCC_W4,na.rm=TRUE),
#                MEAN_IYCC_W4=mean(IYCC_W4,na.rm=TRUE),
#                MIN_IYCC_W4=min(IYCC_W4,na.rm=TRUE),
#                MAX_IYCC_W4=max(IYCC_W4,na.rm=TRUE),   
#                
#                MED_YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
#                MEAN_YLPH_W2=mean(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MIN_YLPH_W2=min(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MAX_YLPH_W2=max(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                
#                MED_YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MEAN_YLPH_W3=mean(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MIN_YLPH_W3=min(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MAX_YLPH_W3=max(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                
#                MED_YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MEAN_YLPH_W4=mean(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MIN_YLPH_W4=min(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                MAX_YLPH_W4=max(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
#                
#                CELLS_W=sum(CELLS,na.rm=TRUE))
# 
# 
# setwd("~/Dropbox/climate change/food security/climate and crop pressure MS/data/ascii_crops_hires")
# C_TO_R<-data.frame(read.csv("COUNTRY_TO_REGION.csv", header=TRUE))
# MAIZE_c<-merge(C_TO_R,MAIZE_c, by=c("NAME"),all=TRUE,sort=TRUE)
# RICE_c<-merge(C_TO_R,RICE_c, by=c("NAME"),all=TRUE,sort=TRUE)
# WHEAT_c<-merge(C_TO_R,WHEAT_c, by=c("NAME"),all=TRUE,sort=TRUE)
# C_TO_R$NUM<-1
# SR_TO_R<-summarySE(C_TO_R,measurevar = "NUM",groupvars=c("Region","Subregion"))
# SR_TO_R<-SR_TO_R[1:2]
# 
# MED_sr<-merge(SR_TO_R,MED_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
# MAIZE_sr<-merge(SR_TO_R,MAIZE_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
# RICE_sr<-merge(SR_TO_R,RICE_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
# WHEAT_sr<-merge(SR_TO_R,WHEAT_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
# 
# setwd("~/Dropbox/climate change/pest MS shared")
# 
# write.table(MED_r, file = "Medians by region.csv", sep = ",", col.names = NA)
# write.table(MED_sr, file = "Medians by subregion.csv", sep = ",", col.names = NA)
# write.table(MED_c, file = "Medians by country.csv", sep = ",", col.names = NA)
# 
# write.table(MAIZE_c, file = "MAIZE by country.csv", sep = ",", col.names = NA)
# write.table(MAIZE_sr, file = "MAIZE by subregion.csv", sep = ",", col.names = NA)
# write.table(MAIZE_r, file = "MAIZE by region.csv", sep = ",", col.names = NA)
# write.table(RICE_c, file = "RICE by country.csv", sep = ",", col.names = NA)
# write.table(RICE_sr, file = "RICE by subregion.csv", sep = ",", col.names = NA)
# write.table(RICE_r, file = "RICE by region.csv", sep = ",", col.names = NA)
# write.table(WHEAT_c, file = "WHEAT by country.csv", sep = ",", col.names = NA)
# write.table(WHEAT_sr, file = "WHEAT by subregion.csv", sep = ",", col.names = NA)
# write.table(WHEAT_r, file = "WHEAT by region.csv", sep = ",", col.names = NA)
# 
# #################  GRAPHICS FROM SUMMARY DATA
# 
# 
# IPMscale<-c(0,.6)
# IYCCscale<-c(0,.07)
# CLP2050scale<-c(0,.3)
# TOT_CL2050scale<-c(0,200000)
# YLPHscale<-c(0,.2)
# IPMaxis<-"median change in Insect Population Metabolism"
# IYCCaxis<-"median proportion of crop lost due to climate change"
# CLP2050axis<-"predicted median proportion of crop lost due to pests [2050-2070]"
# TOT_CL2050axis<-"predicted annual total regional loss of crops (in tonnes) due to climate change alone [2050-2070]"
# YLPHaxis<-"predicted yeild lost due to climate change alone, in tonnes per ha planted [2050-2070]"
# #cycle through variable picks for different graphics (eg. change IPM to IYCC ), use scaling above.
# 
# 
# #MAIZE
# myvars<-c("Region","Subregion","TOT_CL2050_M2","TOT_CL2050_M3","TOT_CL2050_M4")
# SUM1<-MED_sr[myvars]
# SUM1<- rename.vars(SUM1, c("TOT_CL2050_M2","TOT_CL2050_M3","TOT_CL2050_M4"), c("ymin","ymed","ymax"))
# SUM1$Subregion<-with(SUM1,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sem <- ggplot(SUM1, aes(colour=Region, y=ymed, x=Subregion)) 
# sem <- sem+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = TOT_CL2050scale)
# sem <- sem + xlab("Subregion") + ylab(TOT_CL2050axis) 
# sem
# + # Set axis labels
#   #theme_bw() +    
#   #opts(legend.position=c(.9, .6),axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
#   #     axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),
#   #     legend.text=theme_text(size=14),legend.title=theme_text(size=14))
# 
# 
# #RICE
# myvars<-c("Region","Subregion","TOT_CL2050_R2","TOT_CL2050_R3","TOT_CL2050_R4")
# SUM2<-MED_sr[myvars]
# SUM2<- rename.vars(SUM2, c("TOT_CL2050_R2","TOT_CL2050_R3","TOT_CL2050_R4"), c("ymin","ymed","ymax"))
# SUM2$Subregion<-with(SUM2,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# ser <- ggplot(SUM2, aes(colour=Region, y=ymed, x=Subregion)) 
# ser <- ser+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = TOT_CL2050scale)
# ser <- ser + xlab("Subregion") + ylab(TOT_CL2050axis) 
# ser
# + # Set axis labels
#   theme_bw() +    
#   opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
#        axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),legend.position="none")
# 
# #WHEAT
# myvars<-c("Region","Subregion","TOT_CL2050_W2","TOT_CL2050_W3","TOT_CL2050_W4")
# SUM3<-MED_sr[myvars]
# SUM3<- rename.vars(SUM3, c("TOT_CL2050_W2","TOT_CL2050_W3","TOT_CL2050_W4"), c("ymin","ymed","ymax"))
# SUM3$Subregion<-with(SUM3,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sew <- ggplot(SUM3, aes(colour=Region, y=ymed, x=Subregion)) 
# sew <- sew+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = TOT_CL2050scale)
# sew <- sew + xlab("Subregion") + ylab(TOT_CL2050axis)
# sew
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #      axis.title.y  = theme_blank(),axis.title.x  = theme_text(size=14),legend.position="none")
# 
# 
# dev.off()
# png(file="~/Dropbox/climate change/pest MS shared/sr_med_crop_loss_2dC.png", width=900, height=1400) # CHECK RESOLUTION 
# multiplot(sem,ser,sew, cols = 1)
# dev.off()
# #to set range
# 
# #MAIZE
# myvars<-c("Region","Subregion","IPM_M2","IPM_M3","IPM_M4")
# SUM1<-MED_sr[myvars]
# SUM1<- rename.vars(SUM1, c("IPM_M2","IPM_M3","IPM_M4"), c("ymin","ymed","ymax"))
# SUM1$Subregion<-with(SUM1,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sem <- ggplot(SUM1, aes(colour=Region, y=ymed, x=Subregion)) 
# sem <- sem+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IPMscale)
# sem <- sem + xlab("Subregion") + ylab(IPMaxis)
# sem + # Set axis labels
# #  theme_bw() +    
# #  opts(legend.position=c(.9, .6),axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),
# #       legend.text=theme_text(size=14),legend.title=theme_text(size=14))
# 
# #RICE
# myvars<-c("Region","Subregion","IPM_R2","IPM_R3","IPM_R4")
# SUM2<-MED_sr[myvars]
# SUM2<- rename.vars(SUM2, c("IPM_R2","IPM_R3","IPM_R4"), c("ymin","ymed","ymax"))
# SUM2$Subregion<-with(SUM2,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# ser <- ggplot(SUM2, aes(colour=Region, y=ymed, x=Subregion)) 
# ser <- ser+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IPMscale)
# ser <- ser + xlab("Subregion") + ylab(IPMaxis)  # Set axis labels
# ser
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),legend.position="none")
# 
# #WHEAT
# myvars<-c("Region","Subregion","IPM_W2","IPM_W3","IPM_W4")
# SUM3<-MED_sr[myvars]
# SUM3<- rename.vars(SUM3, c("IPM_W2","IPM_W3","IPM_W4"), c("ymin","ymed","ymax"))
# SUM3$Subregion<-with(SUM3,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sew <- ggplot(SUM3, aes(colour=Region, y=ymed, x=Subregion)) 
# sew <- sew+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IPMscale)
# sew <- sew + xlab("Subregion") + ylab(IPMaxis) # Set axis labels
# sew
# #theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_text(size=14),legend.position="none")
# 
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/sr_med_deltaIPM_pest_pressure.png", width=900, height=1400) # CHECK RESOLUTION 
# multiplot(sem,ser,sew, cols = 1)
# dev.off()
# 
# #MAIZE
# myvars<-c("Region","Subregion","IYCC_M2","IYCC_M3","IYCC_M4")
# SUM1<-MED_sr[myvars]
# SUM1<- rename.vars(SUM1, c("IYCC_M2","IYCC_M3","IYCC_M4"), c("ymin","ymed","ymax"))
# SUM1$Subregion<-with(SUM1,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sem <- ggplot(SUM1, aes(colour=Region, y=ymed, x=Subregion)) 
# sem <- sem+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IYCCscale)
# sem <- sem + xlab("Subregion") + ylab(IYCCaxis) # Set axis labels
# sem
# theme_bw() +    
#   opts(legend.position=c(.9, .6),axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
#        axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),
#        legend.text=theme_text(size=14),legend.title=theme_text(size=14))
# 
# #RICE
# myvars<-c("Region","Subregion","IYCC_R2","IYCC_R3","IYCC_R4")
# SUM2<-MED_sr[myvars]
# SUM2<- rename.vars(SUM2, c("IYCC_R2","IYCC_R3","IYCC_R4"), c("ymin","ymed","ymax"))
# SUM2$Subregion<-with(SUM2,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# ser <- ggplot(SUM2, aes(colour=Region, y=ymed, x=Subregion)) 
# ser <- ser+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IYCCscale)
# ser <- ser + xlab("Subregion") + ylab(IYCCaxis)
# ser
# + # Set axis labels
#   theme_bw() +    
#   opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
#        axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),legend.position="none")
# 
# #WHEAT
# myvars<-c("Region","Subregion","IYCC_W2","IYCC_W3","IYCC_W4")
# SUM3<-MED_sr[myvars]
# SUM3<- rename.vars(SUM3, c("IYCC_W2","IYCC_W3","IYCC_W4"), c("ymin","ymed","ymax"))
# SUM3$Subregion<-with(SUM3,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sew <- ggplot(SUM3, aes(colour=Region, y=ymed, x=Subregion)) 
# sew <- sew+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = IYCCscale)
# sew <- sew + xlab("Subregion") + ylab(IYCCaxis) 
# sew
# + # Set axis labels
#   theme_bw() +    
#   opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
#        axis.title.y  = theme_blank(),axis.title.x  = theme_text(size=14),legend.position="none")
# 
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/sr_med_IYCC_PROP_LOST_CC.png", width=900, height=1400) # CHECK RESOLUTION 
# multiplot(sem,ser,sew, cols = 1)
# dev.off()
# 
# 
# #MAIZE
# myvars<-c("Region","Subregion","CLP2050_M2","CLP2050_M3","CLP2050_M4")
# SUM1<-MED_sr[myvars]
# SUM1<- rename.vars(SUM1, c("CLP2050_M2","CLP2050_M3","CLP2050_M4"), c("ymin","ymed","ymax"))
# SUM1$Subregion<-with(SUM1,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sem <- ggplot(SUM1, aes(colour=Region, y=ymed, x=Subregion)) 
# sem <- sem+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = CLP2050scale)
# sem <- sem + xlab("Subregion") + ylab(CLP2050axis)
# sem
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(legend.position=c(.9, .6),axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),
# #       legend.text=theme_text(size=14),legend.title=theme_text(size=14))
# 
# #RICE
# myvars<-c("Region","Subregion","CLP2050_R2","CLP2050_R3","CLP2050_R4")
# SUM2<-MED_sr[myvars]
# SUM2<- rename.vars(SUM2, c("CLP2050_R2","CLP2050_R3","CLP2050_R4"), c("ymin","ymed","ymax"))
# SUM2$Subregion<-with(SUM2,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# ser <- ggplot(SUM2, aes(colour=Region, y=ymed, x=Subregion)) 
# ser <- ser+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = CLP2050scale)
# ser <- ser + xlab("Subregion") + ylab(CLP2050axis)
# ser
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),legend.position="none")
# 
# #WHEAT
# myvars<-c("Region","Subregion","CLP2050_W2","CLP2050_W3","CLP2050_W4")
# SUM3<-MED_sr[myvars]
# SUM3<- rename.vars(SUM3, c("CLP2050_W2","CLP2050_W3","CLP2050_W4"), c("ymin","ymed","ymax"))
# SUM3$Subregion<-with(SUM3,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sew <- ggplot(SUM3, aes(colour=Region, y=ymed, x=Subregion)) 
# sew <- sew+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = CLP2050scale)
# sew <- sew + xlab("Subregion") + ylab(CLP2050axis)
# sew
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_text(size=14),legend.position="none")
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/sr_med_prop_yield_loss_2dC.png", width=900, height=1400) # CHECK RESOLUTION 
# multiplot(sem,ser,sew, cols = 1)
# dev.off()
# 
# #GOOD TO RUN YLPH as well, but not set up for that currently.
# #MAIZE
# myvars<-c("Region","Subregion","YLPH_M2","YLPH_M3","YLPH_M4")
# SUM1<-MED_sr[myvars]
# SUM1<- rename.vars(SUM1, c("YLPH_M2","YLPH_M3","YLPH_M4"), c("ymin","ymed","ymax"))
# SUM1$Subregion<-with(SUM1,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sem <- ggplot(SUM1, aes(colour=Region, y=ymed, x=Subregion)) 
# sem <- sem+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = YLPHscale)
# sem <- sem + xlab("Subregion") + ylab(YLPHaxis)
# sem
# + # Set axis labels
# #  theme_bw() +    
# #  opts(legend.position=c(.9, .6),axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),
# #       legend.text=theme_text(size=14),legend.title=theme_text(size=14))
# 
# #RICE
# myvars<-c("Region","Subregion","YLPH_R2","YLPH_R3","YLPH_R4")
# SUM2<-MED_sr[myvars]
# SUM2<- rename.vars(SUM2, c("YLPH_R2","YLPH_R3","YLPH_R4"), c("ymin","ymed","ymax"))
# SUM2$Subregion<-with(SUM2,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# ser <- ggplot(SUM2, aes(colour=Region, y=ymed, x=Subregion)) 
# ser <- ser+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = YLPHscale)
# ser <- ser + xlab("Subregion") + ylab(YLPHaxis)
# ser
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_blank(),legend.position="none")
# 
# #WHEAT
# myvars<-c("Region","Subregion","YLPH_W2","YLPH_W3","YLPH_W4")
# SUM3<-MED_sr[myvars]
# SUM3<- rename.vars(SUM3, c("YLPH_W2","YLPH_W3","YLPH_W4"), c("ymin","ymed","ymax"))
# SUM3$Subregion<-with(SUM3,factor(Subregion,levels = c('Melanesia','Australia and New Zealand','Western Europe','Southern Europe','Northern Europe','Eastern Europe','Western Asia',
#                                                       'Southern Asia','South-Eastern Asia','Eastern Asia','Central Asia','South America','Northern America','Central America','Caribbean',
#                                                       'Western Africa','Southern Africa','Northern Africa','Middle Africa','East Africa')))
# limits <- aes(ymax = ymax, ymin = ymin)
# sew <- ggplot(SUM3, aes(colour=Region, y=ymed, x=Subregion)) 
# sew <- sew+geom_crossbar(limits,width=0.6)+coord_flip() + scale_y_continuous(limits = YLPHscale)
# sew <- sew + xlab("Subregion") + ylab(YLPHaxis)
# sew
# #+ # Set axis labels
# #  theme_bw() +    
# #  opts(axis.text.x  = theme_text(size=14),axis.text.y  = theme_text(size=14,hjust=1),
# #       axis.title.y  = theme_blank(),axis.title.x  = theme_text(size=14),legend.position="none")
# 
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/subregion_medians_Yield_loss_tonnes_per_hectare.png", width=900, height=1400) # CHECK RESOLUTION 
# multiplot(sem,ser,sew, cols = 1)
# dev.off()
# 
# #### Bubble charts
# 
# setwd("~/Dropbox/climate change/pest MS shared")
# 
# MAIZE_c<-read.csv("MAIZE by country.csv", header=TRUE)
# MAIZE_c<-MAIZE_c[complete.cases(MAIZE_c),]
# RICE_c<-read.csv("RICE by country.csv", header=TRUE)
# RICE_c<-RICE_c[complete.cases(RICE_c),]
# WHEAT_c<-read.csv("WHEAT by country.csv", header=TRUE)
# WHEAT_c<-WHEAT_c[complete.cases(WHEAT_c),]
# 
# 
# MAIZE2<-subset(MAIZE_c,TOT_YLD_TOT_M >= max(TOT_YLD_TOT_M*.001)) # only including countries with total yield = .001 of top yield
# MAIZE2$radius <- sqrt(MAIZE2$TOT_YLD_TOT_M/ pi )
# 
# RICE2<-subset(RICE_c,TOT_YLD_TOT_R >= max(TOT_YLD_TOT_R*.001)) # only including countries with total yield = .001 of top yield
# RICE2$radius <- sqrt(RICE2$TOT_YLD_TOT_R/ pi )
# 
# 
# WHEAT2<-subset(WHEAT_c,TOT_YLD_TOT_W >= max(TOT_YLD_TOT_W*.001)) # only including countries with total yield = .001 of top yield
# WHEAT2$radius <- sqrt(WHEAT2$TOT_YLD_TOT_W/ pi )
# #bubble plot code here: http://flowingdata.com/2010/11/23/how-to-make-bubble-charts/
# 
# #setting all bubles proportional to eachother.
# max(MAIZE2$TOT_YLD_TOT_M)
# max(RICE2$TOT_YLD_TOT_R)
# max(WHEAT2$TOT_YLD_TOT_W)
# 
# #this part simply provides size info for making circles for the legand for the plot
# #legand is made in InDesign or Illustrator.
# scale legond1<-max(sqrt(300000000/pi))/max(sqrt(MAIZE2$TOT_YLD_TOT_M/pi))
# scale_legond2<-max(sqrt(100000000/pi))/max(sqrt(300000000/pi))
# scale_legond3<-max(sqrt(30000000/pi))/max(sqrt(300000000/pi))
# scale_legond4<-max(sqrt(3000000/pi))/max(sqrt(300000000/pi))
# 
# R_scale<-max(sqrt(RICE2$TOT_YLD_TOT_R/pi))/max(sqrt(MAIZE2$TOT_YLD_TOT_M/pi))
# W_scale<-max(sqrt(WHEAT2$TOT_YLD_TOT_W/pi))/max(sqrt(MAIZE2$TOT_YLD_TOT_M/pi))
# 
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/MAIZE_Bubble_IPM.png", width=600, height=900)
# symbols(MAIZE2$WMEAN_YLD_HA_M, MAIZE2$WMEAN_IPM_M3, circles=MAIZE2$radius, inches=0.6, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops",
#         bg=ifelse(abs(MAIZE2$MEAN_LAT_M) >= 40,"red", ifelse(abs(MAIZE2$MEAN_LAT_M) >= 30 & abs(MAIZE2$MEAN_LAT_M) < 40,"orange",
#                                                              ifelse(abs(MAIZE2$MEAN_LAT_M) >= 20 & abs(MAIZE2$MEAN_LAT_M) < 30,"yellow",ifelse(abs(MAIZE2$MEAN_LAT_M) >= 10 & abs (MAIZE2$MEAN_LAT_M) < 20,"green","blue")))))  
# #text(MAIZE2$MEAN_YLD_HA_M, MAIZE2$MEAN_IPM_M3, MAIZE2$NAME, cex=0.5)
# dev.off()
# 
# png(file="~/Dropbox/climate change/pest MS shared/MAIZE_Bubble_TOT_LOSS.png", width=600, height=900)
# symbols(MAIZE2$WMEAN_YLD_HA_M, MAIZE2$TOT_CL2050_M3, circles=MAIZE2$radius, inches=0.6, fg="white", xlab="Yield per HA", ylab="crop loss (%) due to climate change",
#         bg=ifelse(abs(MAIZE2$MEAN_LAT_M) >= 40,"red", ifelse(abs(MAIZE2$MEAN_LAT_M) >= 30 & abs(MAIZE2$MEAN_LAT_M) < 40,"orange",
#                                                              ifelse(abs(MAIZE2$MEAN_LAT_M) >= 20 & abs(MAIZE2$MEAN_LAT_M) < 30,"yellow",ifelse(abs(MAIZE2$MEAN_LAT_M) >= 10 & abs (MAIZE2$MEAN_LAT_M) < 20,"green","blue")))))  
# #text(MAIZE2$MEAN_YLD_HA_M, MAIZE2$MEAN_IPM_M3, MAIZE2$NAME, cex=0.5)
# dev.off()
# 
# png(file="~/Dropbox/climate change/pest MS shared/MAIZE_Bubble_CL2015.png", width=600, height=900)
# symbols(MAIZE2$WMEAN_YLD_HA_M, MAIZE2$WMEAN_CL2050_M3, circles=MAIZE2$radius, inches=0.6, fg="white", xlab="Yield per HA", ylab="total crop loss (%) in a 2 deg world",
#         bg=ifelse(abs(MAIZE2$MEAN_LAT_M) >= 40,"red", ifelse(abs(MAIZE2$MEAN_LAT_M) >= 30 & abs(MAIZE2$MEAN_LAT_M) < 40,"orange",
#                                                              ifelse(abs(MAIZE2$MEAN_LAT_M) >= 20 & abs(MAIZE2$MEAN_LAT_M) < 30,"yellow",ifelse(abs(MAIZE2$MEAN_LAT_M) >= 10 & abs (MAIZE2$MEAN_LAT_M) < 20,"green","blue")))))  
# #text(MAIZE2$MEAN_YLD_HA_M, MAIZE2$MEAN_IPM_M3, MAIZE2$NAME, cex=0.5)
# dev.off()
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/MAIZE_Bubble2.png", width=600, height=900)
# symbols(MAIZE2$WMEAN_YLD_HA_M, MAIZE2$WMEAN_IPM_M3, circles=MAIZE2$radius, inches=0.6, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops",
#         bg=ifelse(abs(MAIZE2$MEAN_LAT_M) >= 40,"red", ifelse(abs(MAIZE2$MEAN_LAT_M) >= 30 & abs(MAIZE2$MEAN_LAT_M) < 40,"orange",
#                                                              ifelse(abs(MAIZE2$MEAN_LAT_M) >= 20 & abs(MAIZE2$MEAN_LAT_M) < 30,"yellow",ifelse(abs(MAIZE2$MEAN_LAT_M) >= 10 & abs (MAIZE2$MEAN_LAT_M) < 20,"green","blue")))))  
# dev.off()
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/RICE_Bubble.png", width=600, height=900)
# symbols(RICE2$WMEAN_YLD_HA_R, RICE2$WMEAN_IPM_R3, circles=RICE2$radius, inches=0.6*R_scale, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops",ylim = c(0,.5), xlim = c(0,10),
#         bg=ifelse(abs(RICE2$MEAN_LAT_R) >= 40,"red", ifelse(abs(RICE2$MEAN_LAT_R) >= 30 & abs(RICE2$MEAN_LAT_R) < 40,"orange",
#                                                             ifelse(abs(RICE2$MEAN_LAT_R) >= 20 & abs(RICE2$MEAN_LAT_R) < 30,"yellow",ifelse(abs(RICE2$MEAN_LAT_R) >= 10 & abs (RICE2$MEAN_LAT_R) < 20,"green","blue")))))  
# dev.off()
# png(file="~/Dropbox/climate change/pest MS shared/RICE_Bubble2.png", width=600, height=900)
# symbols(RICE2$WMEAN_YLD_HA_R, RICE2$WMEAN_IPM_R3, circles=RICE2$radius, inches=0.6*R_scale, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops",
#         bg=ifelse(abs(RICE2$MEAN_LAT_R) >= 40,"red", ifelse(abs(RICE2$MEAN_LAT_R) >= 30 & abs(RICE2$MEAN_LAT_R) < 40,"orange",
#                                                             ifelse(abs(RICE2$MEAN_LAT_R) >= 20 & abs(RICE2$MEAN_LAT_R) < 30,"yellow",ifelse(abs(RICE2$MEAN_LAT_R) >= 10 & abs (RICE2$MEAN_LAT_R) < 20,"green","blue")))))  
# dev.off()
# 
# 
# png(file="~/Dropbox/climate change/pest MS shared/WHEAT_Bubble.png", width=600, height=900)
# symbols(WHEAT2$WMEAN_YLD_HA_W, WHEAT2$WMEAN_IPM_W3, circles=WHEAT2$radius, inches=0.6*W_scale, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops", ylim = c(0,.5),xlim = c(0,10),
#         bg=ifelse(abs(WHEAT2$MEAN_LAT_W) >= 40,"red", ifelse(abs(WHEAT2$MEAN_LAT_W) >= 30 & abs(WHEAT2$MEAN_LAT_W) < 40,"orange",
#                                                              ifelse(abs(WHEAT2$MEAN_LAT_W) >= 20 & abs(WHEAT2$MEAN_LAT_W) < 30,"yellow",ifelse(abs(WHEAT2$MEAN_LAT_W) >= 10 & abs (WHEAT2$MEAN_LAT_W) < 20,"green","blue")))))
# 
# dev.off()
# png(file="~/Dropbox/climate change/pest MS shared/WHEAT_Bubble2.png", width=600, height=900)
# symbols(WHEAT2$WMEAN_YLD_HA_W, WHEAT2$WMEAN_IPM_W3, circles=WHEAT2$radius, inches=0.6*W_scale, fg="white", xlab="Yield per HA", ylab="Increase in insect damage to crops", ylim = c(0,.6),
#         bg=ifelse(abs(WHEAT2$MEAN_LAT_W) >= 40,"red", ifelse(abs(WHEAT2$MEAN_LAT_W) >= 30 & abs(WHEAT2$MEAN_LAT_W) < 40,"orange",
#                                                              ifelse(abs(WHEAT2$MEAN_LAT_W) >= 20 & abs(WHEAT2$MEAN_LAT_W) < 30,"yellow",ifelse(abs(WHEAT2$MEAN_LAT_W) >= 10 & abs (WHEAT2$MEAN_LAT_W) < 20,"green","blue")))))  
# text(WHEAT2$MED_YLD_HA_W, WHEAT2$MED_IPM_W3, WHEAT2$NAME, cex=0.5)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# quants
# quants<-c(seq(-50,50,by=10),Inf)                                                                           
# jet.colors <-colorRampPalette(c("blue", "#007FFF", "#7FFF7F", "yellow", "#FF7F00", "red"),space = "Lab")(length(quants))
# 
# text(MAIZE_c$MED_YLD_HA_M, MAIZE_c$MED_IPM_M3, MAIZE_c$NAME, cex=0.5)
# 
# 
# 
# 
# install.packages("googleVis")  ## If you need to install the package
# library(googleVis)
# M<-gvisMotionChart(Fruits, "Fruit", "Year")
# plot(M)
# cat(M$html$chart, file = "tmp.html")
