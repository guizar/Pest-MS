# summary tables ----------------------------------------------------------



setwd("~/Dropbox/climate change/food security/climate and crop pressure MS/data/ascii_crops_hires")

C_TO_R<-data.frame(read.csv("COUNTRY_TO_REGION.csv", header=TRUE))

# now just need to use plyr and ddply to aggregate by region and subregion!
setwd("~/Dropbox/climate change/pest MS shared")

MED_c<-ddply(ALL,.(NAME),summarise,
             IPM_M2=median(IPM_M2,na.rm=TRUE),
             IPM_M3=median(IPM_M3,na.rm=TRUE),
             IPM_M4=median(IPM_M4,na.rm=TRUE),  
             
             IPM_R2=median(IPM_R2,na.rm=TRUE),
             IPM_R3=median(IPM_R3,na.rm=TRUE),
             IPM_R4=median(IPM_R4,na.rm=TRUE),
             
             IPM_W2=median(IPM_W2,na.rm=TRUE),
             IPM_W3=median(IPM_W3,na.rm=TRUE),
             IPM_W4=median(IPM_W4,na.rm=TRUE),
             
             IY_M2=median(IY_M2,na.rm=TRUE),
             IY_M3=median(IY_M3,na.rm=TRUE),
             IY_M4=median(IY_M4,na.rm=TRUE),  
             
             IY_R2=median(IY_R2,na.rm=TRUE),
             IY_R3=median(IY_R3,na.rm=TRUE),
             IY_R4=median(IY_R4,na.rm=TRUE),
             
             IY_W2=median(IY_W2,na.rm=TRUE),
             IY_W3=median(IY_W3,na.rm=TRUE),
             IY_W4=median(IY_W4,na.rm=TRUE),
             
             CLF_M=median(CLF_M,na.rm=TRUE),
             CL_M=median(CL_M,na.rm=TRUE),
             CY_M=median(CY_M,na.rm=TRUE),
             CA_M=median(CA_M,na.rm=TRUE),
             CGS_M=median(CGS_M,na.rm=TRUE),
             
             CLF_R=median(CLF_R,na.rm=TRUE),
             CL_R=median(CL_R,na.rm=TRUE),
             CY_R=median(CY_R,na.rm=TRUE),
             CA_R=median(CA_R,na.rm=TRUE),
             CGS_R=median(CGS_R,na.rm=TRUE),
             
             CLF_W=median(CLF_W,na.rm=TRUE),
             CL_W=median(CL_W,na.rm=TRUE),
             CY_W=median(CY_W,na.rm=TRUE),
             CA_W=median(CA_W,na.rm=TRUE),
             CGS_W=median(CGS_W,na.rm=TRUE),
             
             MET_M2=median(MET_M2,na.rm=TRUE),
             MET_M3=median(MET_M3,na.rm=TRUE),
             MET_M4=median(MET_M4,na.rm=TRUE),  
             
             MET_R2=median(MET_R2,na.rm=TRUE),
             MET_R3=median(MET_R3,na.rm=TRUE),
             MET_R4=median(MET_R4,na.rm=TRUE),
             
             MET_W2=median(MET_W2,na.rm=TRUE),
             MET_W3=median(MET_W3,na.rm=TRUE),
             MET_W4=median(MET_W4,na.rm=TRUE),
             
             POP_M2=median(POP_M2,na.rm=TRUE),
             POP_M3=median(POP_M3,na.rm=TRUE),
             POP_M4=median(POP_M4,na.rm=TRUE),  
             
             POP_R2=median(POP_R2,na.rm=TRUE),
             POP_R3=median(POP_R3,na.rm=TRUE),
             POP_R4=median(POP_R4,na.rm=TRUE),
             
             POP_W2=median(POP_W2,na.rm=TRUE),
             POP_W3=median(POP_W3,na.rm=TRUE),
             POP_W4=median(POP_W4,na.rm=TRUE),
             
             MET_AVG2=median(MET_AVG2,na.rm=TRUE),
             MET_AVG3=median(MET_AVG3,na.rm=TRUE),
             MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
             
             POP_AVG2=median(POP_AVG2,na.rm=TRUE),
             POP_AVG3=median(POP_AVG3,na.rm=TRUE),
             POP_AVG4=median(POP_AVG4,na.rm=TRUE),
             
             IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
             IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
             IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
             
             YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
             YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
             YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
             YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
             YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
             YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
             
             
             CL2050_M2=median(CL2050_M2,na.rm=TRUE),
             CL2050_M3=median(CL2050_M3,na.rm=TRUE),
             CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
             
             CL2050_R2=median(CL2050_R2,na.rm=TRUE),
             CL2050_R3=median(CL2050_R3,na.rm=TRUE),
             CL2050_R4=median(CL2050_R4,na.rm=TRUE),
             
             CL2050_W2=median(CL2050_W2,na.rm=TRUE),
             CL2050_W3=median(CL2050_W3,na.rm=TRUE),
             CL2050_W4=median(CL2050_W4,na.rm=TRUE),
             
             CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
             CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
             CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
             
             CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
             CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
             CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
             
             CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
             CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
             CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
             
             IYCC_M2=median(IYCC_M2,na.rm=TRUE),
             IYCC_M3=median(IYCC_M3,na.rm=TRUE),
             IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
             
             IYCC_R2=median(IYCC_R2,na.rm=TRUE),
             IYCC_R3=median(IYCC_R3,na.rm=TRUE),
             IYCC_R4=median(IYCC_R4,na.rm=TRUE),
             
             IYCC_W2=median(IYCC_W2,na.rm=TRUE),
             IYCC_W3=median(IYCC_W3,na.rm=TRUE),
             IYCC_W4=median(IYCC_W4,na.rm=TRUE),
             
             YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
             YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
             YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
             
             YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
             YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
             YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
             
             YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
             YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
             YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))


#MEDIEANS BY SUBREGION

MED_sr<-ddply(ALL,.(Subregion),summarise,
              IPM_M2=median(IPM_M2,na.rm=TRUE),
              IPM_M3=median(IPM_M3,na.rm=TRUE),
              IPM_M4=median(IPM_M4,na.rm=TRUE),  
              
              IPM_R2=median(IPM_R2,na.rm=TRUE),
              IPM_R3=median(IPM_R3,na.rm=TRUE),
              IPM_R4=median(IPM_R4,na.rm=TRUE),
              
              IPM_W2=median(IPM_W2,na.rm=TRUE),
              IPM_W3=median(IPM_W3,na.rm=TRUE),
              IPM_W4=median(IPM_W4,na.rm=TRUE),
              
              IY_M2=median(IY_M2,na.rm=TRUE),
              IY_M3=median(IY_M3,na.rm=TRUE),
              IY_M4=median(IY_M4,na.rm=TRUE),  
              
              IY_R2=median(IY_R2,na.rm=TRUE),
              IY_R3=median(IY_R3,na.rm=TRUE),
              IY_R4=median(IY_R4,na.rm=TRUE),
              
              IY_W2=median(IY_W2,na.rm=TRUE),
              IY_W3=median(IY_W3,na.rm=TRUE),
              IY_W4=median(IY_W4,na.rm=TRUE),
              
              CLF_M=median(CLF_M,na.rm=TRUE),
              CL_M=median(CL_M,na.rm=TRUE),
              CY_M=median(CY_M,na.rm=TRUE),
              CA_M=median(CA_M,na.rm=TRUE),
              CGS_M=median(CGS_M,na.rm=TRUE),
              
              CLF_R=median(CLF_R,na.rm=TRUE),
              CL_R=median(CL_R,na.rm=TRUE),
              CY_R=median(CY_R,na.rm=TRUE),
              CA_R=median(CA_R,na.rm=TRUE),
              CGS_R=median(CGS_R,na.rm=TRUE),
              
              CLF_W=median(CLF_W,na.rm=TRUE),
              CL_W=median(CL_W,na.rm=TRUE),
              CY_W=median(CY_W,na.rm=TRUE),
              CA_W=median(CA_W,na.rm=TRUE),
              CGS_W=median(CGS_W,na.rm=TRUE),
              
              MET_M2=median(MET_M2,na.rm=TRUE),
              MET_M3=median(MET_M3,na.rm=TRUE),
              MET_M4=median(MET_M4,na.rm=TRUE),  
              
              MET_R2=median(MET_R2,na.rm=TRUE),
              MET_R3=median(MET_R3,na.rm=TRUE),
              MET_R4=median(MET_R4,na.rm=TRUE),
              
              MET_W2=median(MET_W2,na.rm=TRUE),
              MET_W3=median(MET_W3,na.rm=TRUE),
              MET_W4=median(MET_W4,na.rm=TRUE),
              
              POP_M2=median(POP_M2,na.rm=TRUE),
              POP_M3=median(POP_M3,na.rm=TRUE),
              POP_M4=median(POP_M4,na.rm=TRUE),  
              
              POP_R2=median(POP_R2,na.rm=TRUE),
              POP_R3=median(POP_R3,na.rm=TRUE),
              POP_R4=median(POP_R4,na.rm=TRUE),
              
              POP_W2=median(POP_W2,na.rm=TRUE),
              POP_W3=median(POP_W3,na.rm=TRUE),
              POP_W4=median(POP_W4,na.rm=TRUE),
              
              MET_AVG2=median(MET_AVG2,na.rm=TRUE),
              MET_AVG3=median(MET_AVG3,na.rm=TRUE),
              MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
              
              POP_AVG2=median(POP_AVG2,na.rm=TRUE),
              POP_AVG3=median(POP_AVG3,na.rm=TRUE),
              POP_AVG4=median(POP_AVG4,na.rm=TRUE),
              
              IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
              IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
              IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
              
              YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
              YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
              YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
              YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
              YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
              YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
              
              
              CL2050_M2=median(CL2050_M2,na.rm=TRUE),
              CL2050_M3=median(CL2050_M3,na.rm=TRUE),
              CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
              
              CL2050_R2=median(CL2050_R2,na.rm=TRUE),
              CL2050_R3=median(CL2050_R3,na.rm=TRUE),
              CL2050_R4=median(CL2050_R4,na.rm=TRUE),
              
              CL2050_W2=median(CL2050_W2,na.rm=TRUE),
              CL2050_W3=median(CL2050_W3,na.rm=TRUE),
              CL2050_W4=median(CL2050_W4,na.rm=TRUE),
              
              CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
              CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
              CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
              
              CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
              CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
              CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
              
              CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
              CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
              CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
              
              IYCC_M2=median(IYCC_M2,na.rm=TRUE),
              IYCC_M3=median(IYCC_M3,na.rm=TRUE),
              IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
              
              IYCC_R2=median(IYCC_R2,na.rm=TRUE),
              IYCC_R3=median(IYCC_R3,na.rm=TRUE),
              IYCC_R4=median(IYCC_R4,na.rm=TRUE),
              
              IYCC_W2=median(IYCC_W2,na.rm=TRUE),
              IYCC_W3=median(IYCC_W3,na.rm=TRUE),
              IYCC_W4=median(IYCC_W4,na.rm=TRUE),

              YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
              YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
              YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
              
              YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
              YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              
YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))


#MEDIANS by region

MED_r<-ddply(ALL,.(Region),summarise,
             IPM_M2=median(IPM_M2,na.rm=TRUE),
             IPM_M3=median(IPM_M3,na.rm=TRUE),
             IPM_M4=median(IPM_M4,na.rm=TRUE),  
             
             IPM_R2=median(IPM_R2,na.rm=TRUE),
             IPM_R3=median(IPM_R3,na.rm=TRUE),
             IPM_R4=median(IPM_R4,na.rm=TRUE),
             
             IPM_W2=median(IPM_W2,na.rm=TRUE),
             IPM_W3=median(IPM_W3,na.rm=TRUE),
             IPM_W4=median(IPM_W4,na.rm=TRUE),
             
             IY_M2=median(IY_M2,na.rm=TRUE),
             IY_M3=median(IY_M3,na.rm=TRUE),
             IY_M4=median(IY_M4,na.rm=TRUE),  
             
             IY_R2=median(IY_R2,na.rm=TRUE),
             IY_R3=median(IY_R3,na.rm=TRUE),
             IY_R4=median(IY_R4,na.rm=TRUE),
             
             IY_W2=median(IY_W2,na.rm=TRUE),
             IY_W3=median(IY_W3,na.rm=TRUE),
             IY_W4=median(IY_W4,na.rm=TRUE),
             
             CLF_M=median(CLF_M,na.rm=TRUE),
             CL_M=median(CL_M,na.rm=TRUE),
             CY_M=median(CY_M,na.rm=TRUE),
             CA_M=median(CA_M,na.rm=TRUE),
             CGS_M=median(CGS_M,na.rm=TRUE),
             
             CLF_R=median(CLF_R,na.rm=TRUE),
             CL_R=median(CL_R,na.rm=TRUE),
             CY_R=median(CY_R,na.rm=TRUE),
             CA_R=median(CA_R,na.rm=TRUE),
             CGS_R=median(CGS_R,na.rm=TRUE),
             
             CLF_W=median(CLF_W,na.rm=TRUE),
             CL_W=median(CL_W,na.rm=TRUE),
             CY_W=median(CY_W,na.rm=TRUE),
             CA_W=median(CA_W,na.rm=TRUE),
             CGS_W=median(CGS_W,na.rm=TRUE),
             
             MET_M2=median(MET_M2,na.rm=TRUE),
             MET_M3=median(MET_M3,na.rm=TRUE),
             MET_M4=median(MET_M4,na.rm=TRUE),  
             
             MET_R2=median(MET_R2,na.rm=TRUE),
             MET_R3=median(MET_R3,na.rm=TRUE),
             MET_R4=median(MET_R4,na.rm=TRUE),
             
             MET_W2=median(MET_W2,na.rm=TRUE),
             MET_W3=median(MET_W3,na.rm=TRUE),
             MET_W4=median(MET_W4,na.rm=TRUE),
             
             POP_M2=median(POP_M2,na.rm=TRUE),
             POP_M3=median(POP_M3,na.rm=TRUE),
             POP_M4=median(POP_M4,na.rm=TRUE),  
             
             POP_R2=median(POP_R2,na.rm=TRUE),
             POP_R3=median(POP_R3,na.rm=TRUE),
             POP_R4=median(POP_R4,na.rm=TRUE),
             
             POP_W2=median(POP_W2,na.rm=TRUE),
             POP_W3=median(POP_W3,na.rm=TRUE),
             POP_W4=median(POP_W4,na.rm=TRUE),
             
             MET_AVG2=median(MET_AVG2,na.rm=TRUE),
             MET_AVG3=median(MET_AVG3,na.rm=TRUE),
             MET_AVG4=median(MET_AVG4,na.rm=TRUE),  
             
             POP_AVG2=median(POP_AVG2,na.rm=TRUE),
             POP_AVG3=median(POP_AVG3,na.rm=TRUE),
             POP_AVG4=median(POP_AVG4,na.rm=TRUE),
             
             IPM_AVG2=median(IPM_AVG2,na.rm=TRUE),
             IPM_AVG3=median(IPM_AVG3,na.rm=TRUE),
             IPM_AVG4=median(IPM_AVG4,na.rm=TRUE),
             
             YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
             YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
             YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
             YLD_TOT_M=median(YLD_TOT_M,na.rm=TRUE),
             YLD_TOT_R=median(YLD_TOT_R,na.rm=TRUE),
             YLD_TOT_W=median(YLD_TOT_W,na.rm=TRUE),
             
             
             CL2050_M2=median(CL2050_M2,na.rm=TRUE),
             CL2050_M3=median(CL2050_M3,na.rm=TRUE),
             CL2050_M4=median(CL2050_M4,na.rm=TRUE),  
             
             CL2050_R2=median(CL2050_R2,na.rm=TRUE),
             CL2050_R3=median(CL2050_R3,na.rm=TRUE),
             CL2050_R4=median(CL2050_R4,na.rm=TRUE),
             
             CL2050_W2=median(CL2050_W2,na.rm=TRUE),
             CL2050_W3=median(CL2050_W3,na.rm=TRUE),
             CL2050_W4=median(CL2050_W4,na.rm=TRUE),
             
             CLP2050_M2=median(CLP2050_M2,na.rm=TRUE),
             CLP2050_M3=median(CLP2050_M3,na.rm=TRUE),
             CLP2050_M4=median(CLP2050_M4,na.rm=TRUE),  
             
             CLP2050_R2=median(CLP2050_R2,na.rm=TRUE),
             CLP2050_R3=median(CLP2050_R3,na.rm=TRUE),
             CLP2050_R4=median(CLP2050_R4,na.rm=TRUE),
             
             CLP2050_W2=median(CLP2050_W2,na.rm=TRUE),
             CLP2050_W3=median(CLP2050_W3,na.rm=TRUE),
             CLP2050_W4=median(CLP2050_W4,na.rm=TRUE),
             
             IYCC_M2=median(IYCC_M2,na.rm=TRUE),
             IYCC_M3=median(IYCC_M3,na.rm=TRUE),
             IYCC_M4=median(IYCC_M4,na.rm=TRUE),  
             
             IYCC_R2=median(IYCC_R2,na.rm=TRUE),
             IYCC_R3=median(IYCC_R3,na.rm=TRUE),
             IYCC_R4=median(IYCC_R4,na.rm=TRUE),
             
             IYCC_W2=median(IYCC_W2,na.rm=TRUE),
             IYCC_W3=median(IYCC_W3,na.rm=TRUE),
             IYCC_W4=median(IYCC_W4,na.rm=TRUE),

YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),

YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),

YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE))

CL_c<-ddply(ALL,.(NAME),summarise,                       
            MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
            MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
            MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
            TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
            TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
            TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
            TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
            TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
            TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
            TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
            TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
            TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))


#MEDIEANS BY SUBREGION

CL_sr<-ddply(ALL,.(Subregion),summarise,
             MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
             MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
             MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
             TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
             TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
             TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
             TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
             TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
             TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
             TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
             TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
             TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
             TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
             TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
             TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))

#MEDIANS by region

CL_r<-ddply(ALL,.(Region),summarise,
            MEAN_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
            MEAN_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
            MEAN_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
            TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
            TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
            TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
            TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
            TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
            TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
            TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),         
            TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
            TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
            TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE))

MED_r<-merge(MED_r,CL_r,by="Region",all=TRUE,sort=TRUE)
MED_sr<-merge(MED_sr,CL_sr,by="Subregion",all=TRUE,sort=TRUE)
MED_c<-merge(MED_c,CL_c,by="NAME",all=TRUE,sort=TRUE)


# country tables ----------------------------------------------------------


#making a datatables for IPM and yield impact by country




MAIZE<-ALL[complete.cases(ALL[ ,7:39]),] # subsets ALL to give only maize data

#Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
TOT_YLD_MAIZE<-ddply(MAIZE,.(NAME),summarise,                       
               TOT_YLD_TOT_M=sum(YLD_TOT_M))

MAIZE<-merge(MAIZE,TOT_YLD_MAIZE, by=c("NAME"),all=TRUE,sort=TRUE)
MAIZE_c<-ddply(MAIZE,.(NAME),summarise,                       
               MEAN_LAT_M=mean(LAT),
               WMEAN_LAT_M=sum((LAT*YLD_TOT_M)/TOT_YLD_TOT_M),
               MED_YLD_HA_M=median(YLD_HA_M),
               MEAN_YLD_HA_M=mean(YLD_HA_M),
               WMEAN_YLD_HA_M=sum((YLD_HA_M*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_YLD_HA_M=min(YLD_HA_M),
               MAX_YLD_HA_M=max(YLD_HA_M),
               
               TOT_YLD_TOT_M=sum(YLD_TOT_M),
               
               TOT_CL2050_M2=sum(CL2050_M2), 
               TOT_CL2050_M3=sum(CL2050_M3),
               TOT_CL2050_M4=sum(CL2050_M4),  
               
               MED_IPM_M2=median(IPM_M2),
               MEAN_IPM_M2=mean(IPM_M2),
               WMEAN_IPM_M2=sum((IPM_M2*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_IPM_M2=min(IPM_M2),
               MAX_IPM_M2=max(IPM_M2),
               
               MED_IPM_M3=median(IPM_M3),
               MEAN_IPM_M3=mean(IPM_M3),
               WMEAN_IPM_M3=sum((IPM_M3*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_IPM_M3=min(IPM_M3),
               MAX_IPM_M3=max(IPM_M3),
               
               MED_IPM_M4=median(IPM_M4),
               MEAN_IPM_M4=mean(IPM_M4),
               WMEAN_IPM_M4=sum((IPM_M4*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_IPM_M4=min(IPM_M4),
               MAX_IPM_M4=max(IPM_M4),
               
               MED_IYCC_M2=median(IYCC_M2),
               MEAN_IYCC_M2=mean(IYCC_M2),
               MIN_IYCC_M2=min(IYCC_M2),
               MAX_IYCC_M2=max(IYCC_M2),
               
               MED_IYCC_M3=median(IYCC_M3),
               MEAN_IYCC_M3=mean(IYCC_M3),
               MIN_IYCC_M3=min(IYCC_M3),
               MAX_IYCC_M3=max(IYCC_M3),
               
               MED_IYCC_M4=median(IYCC_M4),
               MEAN_IYCC_M4=mean(IYCC_M4),
               MIN_IYCC_M4=min(IYCC_M4),
               MAX_IYCC_M4=max(IYCC_M4),   
               
               MED_YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M)), 
               MEAN_YLPH_M2=mean(IPM_M2*CLF_M*(CY_M/CA_M)),
               WMEAN_YLPH_M2=sum(((IPM_M2*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_YLPH_M2=min(IPM_M2*CLF_M*(CY_M/CA_M)),
               MAX_YLPH_M2=max(IPM_M2*CLF_M*(CY_M/CA_M)),
               
               MED_YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M)),
               MEAN_YLPH_M3=mean(IPM_M3*CLF_M*(CY_M/CA_M)),
               WMEAN_YLPH_M3=sum(((IPM_M3*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_YLPH_M3=min(IPM_M3*CLF_M*(CY_M/CA_M)),
               MAX_YLPH_M3=max(IPM_M3*CLF_M*(CY_M/CA_M)),
               
               MED_YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M)),
               MEAN_YLPH_M4=mean(IPM_M4*CLF_M*(CY_M/CA_M)),
               WMEAN_YLPH_M4=sum(((IPM_M4*CLF_M*(CY_M/CA_M))*YLD_TOT_M)/TOT_YLD_TOT_M),
               MIN_YLPH_M4=min(IPM_M4*CLF_M*(CY_M/CA_M)),
               MAX_YLPH_M4=max(IPM_M4*CLF_M*(CY_M/CA_M)),
      
               CELLS_M=sum(CELLS))   


RICE<-ALL[complete.cases(ALL[ ,40:72]),] # subsets ALL to give only maize data
#Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
TOT_YLD_RICE<-ddply(RICE,.(NAME),summarise,                       
                     TOT_YLD_TOT_R=sum(YLD_TOT_R))

RICE<-merge(RICE,TOT_YLD_RICE, by=c("NAME"),all=TRUE,sort=TRUE)


RICE_c<-ddply(RICE,.(NAME),summarise,                       
               MEAN_LAT_R=mean(LAT),
               WMEAN_LAT_R=sum((LAT*YLD_TOT_R)/TOT_YLD_TOT_R),
               MED_YLD_HA_R=median(YLD_HA_R),
               MEAN_YLD_HA_R=mean(YLD_HA_R),
               WMEAN_YLD_HA_R=sum((YLD_HA_R*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_YLD_HA_R=min(YLD_HA_R),
               MAX_YLD_HA_R=max(YLD_HA_R),
               
               TOT_YLD_TOT_R=sum(YLD_TOT_R),
               
               TOT_CL2050_R2=sum(CL2050_R2), 
               TOT_CL2050_R3=sum(CL2050_R3),
               TOT_CL2050_R4=sum(CL2050_R4),  
               
               MED_IPM_R2=median(IPM_R2),
               MEAN_IPM_R2=mean(IPM_R2),
               WMEAN_IPM_R2=sum((IPM_R2*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_IPM_R2=min(IPM_R2),
               MAX_IPM_R2=max(IPM_R2),
               
               MED_IPM_R3=median(IPM_R3),
               MEAN_IPM_R3=mean(IPM_R3),
               WMEAN_IPM_R3=sum((IPM_R3*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_IPM_R3=min(IPM_R3),
               MAX_IPM_R3=max(IPM_R3),
               
               MED_IPM_R4=median(IPM_R4),
               MEAN_IPM_R4=mean(IPM_R4),
               WMEAN_IPM_R4=sum((IPM_R4*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_IPM_R4=min(IPM_R4),
               MAX_IPM_R4=max(IPM_R4),
               
               MED_IYCC_R2=median(IYCC_R2),
               MEAN_IYCC_R2=mean(IYCC_R2),
               MIN_IYCC_R2=min(IYCC_R2),
               MAX_IYCC_R2=max(IYCC_R2),
               
               MED_IYCC_R3=median(IYCC_R3),
               MEAN_IYCC_R3=mean(IYCC_R3),
               MIN_IYCC_R3=min(IYCC_R3),
               MAX_IYCC_R3=max(IYCC_R3),
               
               MED_IYCC_R4=median(IYCC_R4),
               MEAN_IYCC_R4=mean(IYCC_R4),
               MIN_IYCC_R4=min(IYCC_R4),
               MAX_IYCC_R4=max(IYCC_R4),   
               
               MED_YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R)), 
               MEAN_YLPH_R2=mean(IPM_R2*CLF_R*(CY_R/CA_R)),
               WMEAN_YLPH_R2=sum(((IPM_R2*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_YLPH_R2=min(IPM_R2*CLF_R*(CY_R/CA_R)),
               MAX_YLPH_R2=max(IPM_R2*CLF_R*(CY_R/CA_R)),
               
               MED_YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R)),
               MEAN_YLPH_R3=mean(IPM_R3*CLF_R*(CY_R/CA_R)),
               WMEAN_YLPH_R3=sum(((IPM_R3*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_YLPH_R3=min(IPM_R3*CLF_R*(CY_R/CA_R)),
               MAX_YLPH_R3=max(IPM_R3*CLF_R*(CY_R/CA_R)),
               
               MED_YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R)),
               MEAN_YLPH_R4=mean(IPM_R4*CLF_R*(CY_R/CA_R)),
               WMEAN_YLPH_R4=sum(((IPM_R4*CLF_R*(CY_R/CA_R))*YLD_TOT_R)/TOT_YLD_TOT_R),
               MIN_YLPH_R4=min(IPM_R4*CLF_R*(CY_R/CA_R)),
               MAX_YLPH_R4=max(IPM_R4*CLF_R*(CY_R/CA_R)),
               
               CELLS_R=sum(CELLS))

WHEAT<-ALL[complete.cases(ALL[ ,73:105]),] # subsets ALL to give only maize data

#Got to first calculate sums, then add them back into maize with merge, and then do these calcs.
#Not sure what the last line means.... hmm

TOT_YLD_WHEAT<-ddply(WHEAT,.(NAME),summarise,                       
                    TOT_YLD_TOT_W=sum(YLD_TOT_W))

WHEAT<-merge(WHEAT,TOT_YLD_WHEAT, by=c("NAME"),all=TRUE,sort=TRUE)



WHEAT_c<-ddply(WHEAT,.(NAME),summarise,                       
              MEAN_LAT_W=mean(LAT),
              WMEAN_LAT_W=sum((LAT*YLD_TOT_W)/TOT_YLD_TOT_W),
              MED_YLD_HA_W=median(YLD_HA_W),
              MEAN_YLD_HA_W=mean(YLD_HA_W),
              WMEAN_YLD_HA_W=sum((YLD_HA_W*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_YLD_HA_W=min(YLD_HA_W),
              MAX_YLD_HA_W=max(YLD_HA_W),
              
              TOT_YLD_TOT_W=sum(YLD_TOT_W),
              
              TOT_CL2050_W2=sum(CL2050_W2), 
              TOT_CL2050_W3=sum(CL2050_W3),
              TOT_CL2050_W4=sum(CL2050_W4),  
              
              MED_IPM_W2=median(IPM_W2),
              MEAN_IPM_W2=mean(IPM_W2),
              WMEAN_IPM_W2=sum((IPM_W2*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_IPM_W2=min(IPM_W2),
              MAX_IPM_W2=max(IPM_W2),
              
              MED_IPM_W3=median(IPM_W3),
              MEAN_IPM_W3=mean(IPM_W3),
              WMEAN_IPM_W3=sum((IPM_W3*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_IPM_W3=min(IPM_W3),
              MAX_IPM_W3=max(IPM_W3),
              
              MED_IPM_W4=median(IPM_W4),
              MEAN_IPM_W4=mean(IPM_W4),
              WMEAN_IPM_W4=sum((IPM_W4*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_IPM_W4=min(IPM_W4),
              MAX_IPM_W4=max(IPM_W4),
              
              MED_IYCC_W2=median(IYCC_W2),
              MEAN_IYCC_W2=mean(IYCC_W2),
              MIN_IYCC_W2=min(IYCC_W2),
              MAX_IYCC_W2=max(IYCC_W2),
              
              MED_IYCC_W3=median(IYCC_W3),
              MEAN_IYCC_W3=mean(IYCC_W3),
              MIN_IYCC_W3=min(IYCC_W3),
              MAX_IYCC_W3=max(IYCC_W3),
              
              MED_IYCC_W4=median(IYCC_W4),
              MEAN_IYCC_W4=mean(IYCC_W4),
              MIN_IYCC_W4=min(IYCC_W4),
              MAX_IYCC_W4=max(IYCC_W4),   
              
              MED_YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W)), 
              MEAN_YLPH_W2=mean(IPM_W2*CLF_W*(CY_W/CA_W)),
              WMEAN_YLPH_W2=sum(((IPM_W2*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_YLPH_W2=min(IPM_W2*CLF_W*(CY_W/CA_W)),
              MAX_YLPH_W2=max(IPM_W2*CLF_W*(CY_W/CA_W)),
              
              MED_YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W)),
              MEAN_YLPH_W3=mean(IPM_W3*CLF_W*(CY_W/CA_W)),
              WMEAN_YLPH_W3=sum(((IPM_W3*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_YLPH_W3=min(IPM_W3*CLF_W*(CY_W/CA_W)),
              MAX_YLPH_W3=max(IPM_W3*CLF_W*(CY_W/CA_W)),
              
              MED_YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W)),
              MEAN_YLPH_W4=mean(IPM_W4*CLF_W*(CY_W/CA_W)),
              WMEAN_YLPH_W4=sum(((IPM_W4*CLF_W*(CY_W/CA_W))*YLD_TOT_W)/TOT_YLD_TOT_W),
              MIN_YLPH_W4=min(IPM_W4*CLF_W*(CY_W/CA_W)),
              MAX_YLPH_W4=max(IPM_W4*CLF_W*(CY_W/CA_W)),
              
              CELLS_W=sum(CELLS))


####BY SUBREGION




MAIZE_r<-ddply(ALL,.(Region),summarise,                       
               MED_YLD_HA_M=median(YLD_HA_M,na.rm=TRUE),
               MEAN_YLD_HA_M=mean(YLD_HA_M,na.rm=TRUE),
               MIN_YLD_HA_M=min(YLD_HA_M,na.rm=TRUE),
               MAX_YLD_HA_M=max(YLD_HA_M,na.rm=TRUE),
               
               TOT_YLD_TOT_M=sum(YLD_TOT_M,na.rm=TRUE),
               
               TOT_CL2050_M2=sum(CL2050_M2,na.rm=TRUE), #total crop loss due to climate change in country
               TOT_CL2050_M3=sum(CL2050_M3,na.rm=TRUE),
               TOT_CL2050_M4=sum(CL2050_M4,na.rm=TRUE),  
               
               MED_IPM_M2=median(IPM_M2,na.rm=TRUE),
               MEAN_IPM_M2=mean(IPM_M2,na.rm=TRUE),
               MIN_IPM_M2=min(IPM_M2,na.rm=TRUE),
               MAX_IPM_M2=max(IPM_M2,na.rm=TRUE),
               
               MED_IPM_M3=median(IPM_M3,na.rm=TRUE),
               MEAN_IPM_M3=mean(IPM_M3,na.rm=TRUE),
               MIN_IPM_M3=min(IPM_M3,na.rm=TRUE),
               MAX_IPM_M3=max(IPM_M3,na.rm=TRUE),
               
               MED_IPM_M4=median(IPM_M4,na.rm=TRUE),
               MEAN_IPM_M4=mean(IPM_M4,na.rm=TRUE),
               MIN_IPM_M4=min(IPM_M4,na.rm=TRUE),
               MAX_IPM_M4=max(IPM_M4,na.rm=TRUE),
               
               MED_IYCC_M2=median(IYCC_M2,na.rm=TRUE),
               MEAN_IYCC_M2=mean(IYCC_M2,na.rm=TRUE),
               MIN_IYCC_M2=min(IYCC_M2,na.rm=TRUE),
               MAX_IYCC_M2=max(IYCC_M2,na.rm=TRUE),
               
               MED_IYCC_M3=median(IYCC_M3,na.rm=TRUE),
               MEAN_IYCC_M3=mean(IYCC_M3,na.rm=TRUE),
               MIN_IYCC_M3=min(IYCC_M3,na.rm=TRUE),
               MAX_IYCC_M3=max(IYCC_M3,na.rm=TRUE),
               
               MED_IYCC_M4=median(IYCC_M4,na.rm=TRUE),
               MEAN_IYCC_M4=mean(IYCC_M4,na.rm=TRUE),
               MIN_IYCC_M4=min(IYCC_M4,na.rm=TRUE),
               MAX_IYCC_M4=max(IYCC_M4,na.rm=TRUE),   
               
               MED_YLPH_M2=median(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
               MEAN_YLPH_M2=mean(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MIN_YLPH_M2=min(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MAX_YLPH_M2=max(IPM_M2*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               
               MED_YLPH_M3=median(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MEAN_YLPH_M3=mean(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MIN_YLPH_M3=min(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MAX_YLPH_M3=max(IPM_M3*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               
               MED_YLPH_M4=median(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MEAN_YLPH_M4=mean(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MIN_YLPH_M4=min(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               MAX_YLPH_M4=max(IPM_M4*CLF_M*(CY_M/CA_M),na.rm=TRUE),
               
               CELLS_M=sum(CELLS,na.rm=TRUE))   


RICE_r<-ddply(ALL,.(Region),summarise,                       
              MED_YLD_HA_R=median(YLD_HA_R,na.rm=TRUE),
              MEAN_YLD_HA_R=mean(YLD_HA_R,na.rm=TRUE),
              MIN_YLD_HA_R=min(YLD_HA_R,na.rm=TRUE),
              MAX_YLD_HA_R=max(YLD_HA_R,na.rm=TRUE),
              
              TOT_YLD_TOT_R=sum(YLD_TOT_R,na.rm=TRUE),
              
              TOT_CL2050_R2=sum(CL2050_R2,na.rm=TRUE), #total crop loss due to climate change in country
              TOT_CL2050_R3=sum(CL2050_R3,na.rm=TRUE),
              TOT_CL2050_R4=sum(CL2050_R4,na.rm=TRUE),  
              
              MED_IPM_R2=median(IPM_R2,na.rm=TRUE),
              MEAN_IPM_R2=mean(IPM_R2,na.rm=TRUE),
              MIN_IPM_R2=min(IPM_R2,na.rm=TRUE),
              MAX_IPM_R2=max(IPM_R2,na.rm=TRUE),
              
              MED_IPM_R3=median(IPM_R3,na.rm=TRUE),
              MEAN_IPM_R3=mean(IPM_R3,na.rm=TRUE),
              MIN_IPM_R3=min(IPM_R3,na.rm=TRUE),
              MAX_IPM_R3=max(IPM_R3,na.rm=TRUE),
              
              MED_IPM_R4=median(IPM_R4,na.rm=TRUE),
              MEAN_IPM_R4=mean(IPM_R4,na.rm=TRUE),
              MIN_IPM_R4=min(IPM_R4,na.rm=TRUE),
              MAX_IPM_R4=max(IPM_R4,na.rm=TRUE),
              
              MED_IYCC_R2=median(IYCC_R2,na.rm=TRUE),
              MEAN_IYCC_R2=mean(IYCC_R2,na.rm=TRUE),
              MIN_IYCC_R2=min(IYCC_R2,na.rm=TRUE),
              MAX_IYCC_R2=max(IYCC_R2,na.rm=TRUE),
              
              MED_IYCC_R3=median(IYCC_R3,na.rm=TRUE),
              MEAN_IYCC_R3=mean(IYCC_R3,na.rm=TRUE),
              MIN_IYCC_R3=min(IYCC_R3,na.rm=TRUE),
              MAX_IYCC_R3=max(IYCC_R3,na.rm=TRUE),
              
              MED_IYCC_R4=median(IYCC_R4,na.rm=TRUE),
              MEAN_IYCC_R4=mean(IYCC_R4,na.rm=TRUE),
              MIN_IYCC_R4=min(IYCC_R4,na.rm=TRUE),
              MAX_IYCC_R4=max(IYCC_R4,na.rm=TRUE),   
              
              MED_YLPH_R2=median(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
              MEAN_YLPH_R2=mean(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MIN_YLPH_R2=min(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MAX_YLPH_R2=max(IPM_R2*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              
              MED_YLPH_R3=median(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MEAN_YLPH_R3=mean(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MIN_YLPH_R3=min(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MAX_YLPH_R3=max(IPM_R3*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              
              MED_YLPH_R4=median(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MEAN_YLPH_R4=mean(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MIN_YLPH_R4=min(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              MAX_YLPH_R4=max(IPM_R4*CLF_R*(CY_R/CA_R),na.rm=TRUE),
              
              CELLS_R=sum(CELLS,na.rm=TRUE))

WHEAT_r<-ddply(ALL,.(Region),summarise,                       
               MED_YLD_HA_W=median(YLD_HA_W,na.rm=TRUE),
               MEAN_YLD_HA_W=mean(YLD_HA_W,na.rm=TRUE),
               MIN_YLD_HA_W=min(YLD_HA_W,na.rm=TRUE),
               MAX_YLD_HA_W=max(YLD_HA_W,na.rm=TRUE),
               
               TOT_YLD_TOT_W=sum(YLD_TOT_W,na.rm=TRUE),
               
               TOT_CL2050_W2=sum(CL2050_W2,na.rm=TRUE), #total crop loss due to climate change in country
               TOT_CL2050_W3=sum(CL2050_W3,na.rm=TRUE),
               TOT_CL2050_W4=sum(CL2050_W4,na.rm=TRUE),  
               
               MED_IPM_W2=median(IPM_W2,na.rm=TRUE),
               MEAN_IPM_W2=mean(IPM_W2,na.rm=TRUE),
               MIN_IPM_W2=min(IPM_W2,na.rm=TRUE),
               MAX_IPM_W2=max(IPM_W2,na.rm=TRUE),
               
               MED_IPM_W3=median(IPM_W3,na.rm=TRUE),
               MEAN_IPM_W3=mean(IPM_W3,na.rm=TRUE),
               MIN_IPM_W3=min(IPM_W3,na.rm=TRUE),
               MAX_IPM_W3=max(IPM_W3,na.rm=TRUE),
               
               MED_IPM_W4=median(IPM_W4,na.rm=TRUE),
               MEAN_IPM_W4=mean(IPM_W4,na.rm=TRUE),
               MIN_IPM_W4=min(IPM_W4,na.rm=TRUE),
               MAX_IPM_W4=max(IPM_W4,na.rm=TRUE),
               
               MED_IYCC_W2=median(IYCC_W2,na.rm=TRUE),
               MEAN_IYCC_W2=mean(IYCC_W2,na.rm=TRUE),
               MIN_IYCC_W2=min(IYCC_W2,na.rm=TRUE),
               MAX_IYCC_W2=max(IYCC_W2,na.rm=TRUE),
               
               MED_IYCC_W3=median(IYCC_W3,na.rm=TRUE),
               MEAN_IYCC_W3=mean(IYCC_W3,na.rm=TRUE),
               MIN_IYCC_W3=min(IYCC_W3,na.rm=TRUE),
               MAX_IYCC_W3=max(IYCC_W3,na.rm=TRUE),
               
               MED_IYCC_W4=median(IYCC_W4,na.rm=TRUE),
               MEAN_IYCC_W4=mean(IYCC_W4,na.rm=TRUE),
               MIN_IYCC_W4=min(IYCC_W4,na.rm=TRUE),
               MAX_IYCC_W4=max(IYCC_W4,na.rm=TRUE),   
               
               MED_YLPH_W2=median(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE), # yield (tonnes) lost per ha of planted area due to climate
               MEAN_YLPH_W2=mean(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MIN_YLPH_W2=min(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MAX_YLPH_W2=max(IPM_W2*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               
               MED_YLPH_W3=median(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MEAN_YLPH_W3=mean(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MIN_YLPH_W3=min(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MAX_YLPH_W3=max(IPM_W3*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               
               MED_YLPH_W4=median(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MEAN_YLPH_W4=mean(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MIN_YLPH_W4=min(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               MAX_YLPH_W4=max(IPM_W4*CLF_W*(CY_W/CA_W),na.rm=TRUE),
               
               CELLS_W=sum(CELLS,na.rm=TRUE))


setwd("~/Dropbox/climate change/food security/climate and crop pressure MS/data/ascii_crops_hires")
C_TO_R<-data.frame(read.csv("COUNTRY_TO_REGION.csv", header=TRUE))
MAIZE_c<-merge(C_TO_R,MAIZE_c, by=c("NAME"),all=TRUE,sort=TRUE)
RICE_c<-merge(C_TO_R,RICE_c, by=c("NAME"),all=TRUE,sort=TRUE)
WHEAT_c<-merge(C_TO_R,WHEAT_c, by=c("NAME"),all=TRUE,sort=TRUE)
C_TO_R$NUM<-1
SR_TO_R<-summarySE(C_TO_R,measurevar = "NUM",groupvars=c("Region","Subregion"))
SR_TO_R<-SR_TO_R[1:2]

MED_sr<-merge(SR_TO_R,MED_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
MAIZE_sr<-merge(SR_TO_R,MAIZE_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
RICE_sr<-merge(SR_TO_R,RICE_sr, by=c("Subregion"),all=TRUE,sort=TRUE)
WHEAT_sr<-merge(SR_TO_R,WHEAT_sr, by=c("Subregion"),all=TRUE,sort=TRUE)

setwd("~/Dropbox/climate change/pest MS shared")

write.table(MED_r, file = "Medians by region.csv", sep = ",", col.names = NA)
write.table(MED_sr, file = "Medians by subregion.csv", sep = ",", col.names = NA)
write.table(MED_c, file = "Medians by country.csv", sep = ",", col.names = NA)

write.table(MAIZE_c, file = "MAIZE by country.csv", sep = ",", col.names = NA)
write.table(MAIZE_sr, file = "MAIZE by subregion.csv", sep = ",", col.names = NA)
write.table(MAIZE_r, file = "MAIZE by region.csv", sep = ",", col.names = NA)
write.table(RICE_c, file = "RICE by country.csv", sep = ",", col.names = NA)
write.table(RICE_sr, file = "RICE by subregion.csv", sep = ",", col.names = NA)
write.table(RICE_r, file = "RICE by region.csv", sep = ",", col.names = NA)
write.table(WHEAT_c, file = "WHEAT by country.csv", sep = ",", col.names = NA)
write.table(WHEAT_sr, file = "WHEAT by subregion.csv", sep = ",", col.names = NA)
write.table(WHEAT_r, file = "WHEAT by region.csv", sep = ",", col.names = NA)
