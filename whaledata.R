library(raster)
library(gdalUtils)
library(gdal)

setwd("F:/whales")


###avhrr

#download the avhrr files
urlsfeb <- XML::getHTMLLinks(
  "https://www-ncei-noaa-gov.eres.library.manoa.hawaii.edu/data/oceans/pathfinder/Version5.3/L3C/1985/data/", 
  xpQuery = "//a/@href['198502'=substring(., string-length(.) - 3)]"
)


#open avhrr dataset. sst climatology for 1985-2001
#convert HDF files to TIF
sds <- get_subdatasets("F:/whales/avhrr/month02_combined.hdf")
gdal_translate(sds[4], dst_dataset = "F:/whales/avhrr/month02_combined_hdfout.tif")

r <- raster()
extent(r) <- extent(-180,180,-90,90)
r <- disaggregate(r,fact=24)

avhrrfeb <- raster("F:/whales/avhrr/month02_combined_hdfout.tif")
avhrrFeb <- projectRaster(avhrrfeb,r,method='bilinear',filename="F:/whales/avhrr/month02_combined.tif",overwrite=T)

sds <- get_subdatasets("F:/whales/avhrr/month08_combined.hdf")
gdal_translate(sds[4], dst_dataset = "F:/whales/avhrr/month08_combined_hdfout.tif")

avhrraug <- raster("F:/whales/avhrr/month08_combined_hdfout.tif")
avhrrAug <- projectRaster(avhrraug,r,method='bilinear',filename="F:/whales/avhrr/month08_combined.tif",overwrite=T)

##########









###hist 

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/hist"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
febs <- c(2+12*(0:16)) #17 years between 1985-2001
augs <- c(8+12*(0:16))

#skip #5 CAS-ESM2-0 because not equally spaced grid
#skip #13,14,16,17 GFDL-CM4, GFDL-ESM4, MIROC6,MRI-ESM2-0 for same reason
for(i in 1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #some models are broken up into multiple files. make model into single brick with all years
  r <- raster()
  for(j in 1:length(modelfiles)){
    modelbrick <- brick(modelfiles[j],varname="tos")
    r <- stack(r,modelbrick)
  }
  #keep just the years 1985-2001
  r8501 <- dropLayer(r,c(1:1620,(nlayers(r)-155):nlayers(r))) #CHECK THIS
  #extract just the febs
  febstack <- stack(r8501[[febs[1:17]]])
  #calculate mean feb temp for 1985-2001 (base period) for each model. save
  calc(febstack,fun=mean,filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb.tif"))
  
  #extract just the augs
  augstack <- stack(r8501[[augs[1:17]]])
  #calculate mean aug temp for 1985-2001 (base period) for each model. save
  calc(augstack,fun=mean,filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug.tif"))

}



###ssp245

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/ssp245"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0","IPSL-CM6A-LR","NESM3","NorESM2-LM","NorESM2-MM")
#skip #5 CAS-ESM2-0 because not equally spaced grid
#skip #13,14,16,17 GFDL-CM4, GFDL-ESM4, MIROC6,MRI-ESM2-0 for same reason

febs <- c(2+12*(0:85)) #86 years between 2015-2100 
augs <- c(8+12*(0:85))

for(i in 1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #make single brick with all model years 2015-2100
  r <- raster()
  for(j in 1:length(modelfiles)){
    modelbrick <- brick(modelfiles[j],varname="tos")
    r <- stack(r,modelbrick)
  }
  #extract just the febs. this will be feb temps 2015-2100 (rasterbrick with 86 layers)
  febstack <- stack(r[[febs[1:86]]])
  #calculate mean decadal temps (e.g. 2015-2024, 2025-2034, ...). will result in a 9 layer stack
  index <- c(1:nlayers(febstack))
  febdecadal <- raster()
  for(k in 0:8){
    ten <- c((k*10+1):(k*10+10))
    window <- intersect(index,ten)
    decade <- stack(febstack[[window[1:length(window)]]])
    decademean <- calc(decade,fun=mean,na.rm=T)
    febdecadal <- stack(febdecadal,decademean)
  }
  #calculate delta  - difference between model hist base period and decadal mean. save
  febbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb.tif"))
  overlay(febdecadal,febbase, fun=function(x,y,na.rm=T){return(x-y)},
                      filename=paste0("F:/whales/cmip6/ssp245/",modelnames[i],"_decadedelta_feb.tif"))
  
  #extract just the augs. this will be aug temps 2015-2100 (rasterbrick with 86 layers)
  augstack <- stack(r[[augs[1:86]]])
  #calculate mean decadal temps (e.g. 2016-2025, 2026-2035, ...). will result in a 9 layer stack
  index <- c(1:nlayers(augstack))
  augdecadal <- raster()
  for(l in 0:8){
    ten <- c((l*10+1):(l*10+10))
    window <- intersect(index,ten)
    decade <- stack(augstack[[window[1:length(window)]]])
    decademean <- calc(decade,fun=mean,na.rm=T)
    augdecadal <- stack(augdecadal,decademean)
  }
  #calculate delta  - difference between model hist base period and decadal mean. save
  augbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug.tif"))
  overlay(augdecadal,augbase, fun=function(x,y,na.rm=T){return(x-y)},
                      filename=paste0("F:/whales/cmip6/ssp245/",modelnames[i],"_decadedelta_aug.tif"))
}

extent(rtemplate) <- extent(0,360,-90,90)
#calculate model median delta
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/ssp245",pattern="feb",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("feb245.", i)
  assign(object,r)
}
deltafeb <- overlay(feb245.1,feb245.2,feb245.3,feb245.4,feb245.5,feb245.6,feb245.7,feb245.8,feb245.9,feb245.10,
        feb245.11,feb245.12,feb245.13,feb245.14,feb245.15,
        fun=median,na.rm=T,filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_feb.tif")
#for aug
augmodelfiles <- list.files("F:/whales/cmip6/ssp245",pattern="aug",full.names=T)
for(i in 1:length(augmodelfiles)){
  object <- paste0("aug245.", i)
  r <- brick(augmodelfiles[i])
  assign(object,r)
}
deltaaug <- overlay(aug245.1,aug245.2,aug245.3,aug245.4,aug245.5,aug245.6,aug245.7,aug245.8,aug245.9,aug245.10,
        aug245.11,aug245.12,aug245.13,aug245.14,aug245.15,aug245.16,aug245.17,aug245.18,aug245.19,aug245.20,
        fun=median,na.rm=T,filename=paste0("F:/whales/cmip6/ssp245/modelmedian_decadedelta_aug.tif"))


#add delta to avhrr to get downscaled avhrr sst
#for feb
febsst245 <- overlay(avhrrFeb,deltafeb, fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245.tif")
#for aug
augsst245 <- overlay(avhrrAug,deltaaug, fun=sum,filename="F:/whales/cmip6/ssp245/AugSST_ssp245.tif")






###ssp585

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/ssp585"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
febs <- c(2+12*(0:86)) #86 years between 2015-2100 
augs <- c(8+12*(0:86))

for(i in 1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #make single brick with all model years 2015-2100
  r <- raster()
  for(j in 1:length(modelfiles)){
    modelbrick <- brick(modelfiles[j],varname="tos")
    r <- stack(r,modelbrick)
  }
  #extract just the febs. this will be feb temps 2015-2100 (rasterbrick with 86 layers)
  febstack <- stack(r[[febs]])
  #calculate mean decadal temps (e.g. 2016-2025, 2026-2035, ...). will result in a 9 layer stack
  index <- c(1:nlayers(febstack))
  febdecadal <- raster()
  for(k in 0:8){
    ten <- c((k*10+1):(k*10+10))
    window <- intersect(index,ten)
    decade <- stack(febstack[[window]])
    decademean <- calc(decade,fun=mean,na.rm=T)
    febdecadal <- stack(febdecadal,decademean)
  }
  #calculate delta  - difference between model hist base period and decadal mean. save
  febbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb.tif"))
  overlay(febdecadal,febbase, fun=function(x,y,na.rm=T){return(x-y)},
          filename=paste0("F:/whales/cmip6/ssp585/",modelnames[i],"_decadedelta_feb.tif"))
  
  #extract just the augs. this will be aug temps 2015-2100 (rasterbrick with 86 layers)
  augstack <- stack(r[[augs]])
  #calculate mean decadal temps (e.g. 2016-2025, 2026-2035, ...). will result in a 9 layer stack
  index <- c(1:nlayers(augstack))
  augdecadal <- raster()
  for(l in 0:8){
    ten <- c((l*10+1):(l*10+10))
    window <- intersect(index,ten)
    decade <- stack(augstack[[window]])
    decademean <- calc(decade,fun=mean,na.rm=T)
    augdecadal <- stack(augdecadal,decademean)
  }
  #calculate delta  - difference between model hist base period and decadal mean. save
  augbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug.tif"))
  overlay(augdecadal,augbase, fun=function(x,y,na.rm=T){return(x-y)},
          filename=paste0("F:/whales/cmip6/ssp585/",modelnames[i],"_decadedelta_aug.tif"))
}

#calculate model median delta
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/ssp585",pattern="feb",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("feb585.", i)
  r <- brick(febmodelfiles[i])
  assign(object,r)
}
deltafeb <- overlay(feb585.1,feb585.2,feb585.3,feb585.4,feb585.5,feb585.6,feb585.7,feb585.8,feb585.9,feb585.10,
                    feb585.11,feb585.12,feb585.13,feb585.14,feb585.15,feb585.16,feb585.17,feb585.18,feb585.19,feb585.20,
                    fun=median,na.rm=T,filename="F:/whales/cmip6/ssp585/modelmedian_decadedelta_feb.tif")
#for aug
augmodelfiles <- list.files("F:/whales/cmip6/ssp585",pattern="aug",full.names=T)
for(i in 1:length(augmodelfiles)){
  object <- paste0("aug585.", i)
  r <- brick(augmodelfiles[i])
  assign(object,r)
}
deltaaug <- overlay(aug585.1,aug585.2,aug585.3,aug585.4,aug585.5,aug585.6,aug585.7,aug585.8,aug585.9,aug585.10,
                    aug585.11,aug585.12,aug585.13,aug585.14,aug585.15,aug585.16,aug585.17,aug585.18,aug585.19,aug585.20,
                    fun=median,na.rm=T,filename=paste0("F:/whales/cmip6/ssp585/modelmedian_decadedelta_aug.tif"))


#add delta to avhrr to get downscaled avhrr sst
#for feb
febsst585 <- overlay(avhrrfeb,deltafeb, fun=sum,filename="F:/whales/cmip6/ssp585/FebSST_ssp585.tif")
#for aug
augsst585 <- overlay(avhrraug,deltaaug, fun=sum,filename="F:/whales/cmip6/ssp585/AugSST_ssp585.tif")





























###
#old version

setwd("C:/Users/rsett/Documents/whales") 

#import csv files. only keep columns that matter
origcsv <- list.files(path="C:/Users/rsett/Documents/whales/IndivData-CSVfmt",pattern='*.csv$',full.names=T)
origcsvname <- paste("csv",1:9,sep="")
  
for (i in 1:length(origcsvname)) {
  x <- read.csv(file=origcsv[i],row.names=NULL,header=T)
  y <- x[-c(1:4,10:15,19,23:40)] #get rid of columns we don't need
  #calculate decimal degrees
  y$lati <- y$Lat + (y$Mn)/60
  y$latitude <- ifelse(y$X=="S",-(y$lati),y$lati)
  y$long <- y$Lon + (y$Mn.1)/60
  y$longitude <- ifelse(y$X.1=="W",-(y$long),y$long)
  #make all length units in metric (feet/inches are L.u=2 or 4)
  y$length <- ifelse(y$L.u==2 , ((y$Len)/100)/3.281,
                     ifelse(y$L.u==4, (as.numeric(substr(y$Len,1,2))/3.281 + as.numeric(substr(y$Len,3,4))/39.37),y$Len))
  y$length[is.na(y$length)] <- 0 #fix NA values
  y$basin <- substr(origcsv[i],50,51) #add column for ocean basin
  z <- y[-c(3,4,6:12,14)]
  assign(origcsvname[i],z)
}

whaledf <- rbind(csv1,csv2,csv3,csv4,csv5,csv6,csv7,csv8,csv9)
write.csv(whaledf,"C:/Users/rsett/Documents/whales/IndivData-CSVfmt/allwhales.csv", row.names=FALSE)










#####
#sst analysis

#use NOAA MMM file
noaa.mmm <- raster("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

#cmip5 rcp85 2006 sst model median
#cmip <- brick("E:/sst/sstRCP85med_mean5.tif") #this is annual mean
cmip <- brick("E:/sst/SSTRCP85_monthly_temp.nc")
cmip2006 <- cmip[[1]]


#cmip5 2100 sst
cmip2100 <- cmip[[1128]]

#difference between 2100 and 2005 both cmip5 (change in sst)
cmipdiff <- overlay(cmip2100,cmip2006,fun=function(x,y,na.rm=T){return(x-y)})
cmipdiffd <- disaggregate(cmipdiff,fact=2)
cmipdiffa <- crop(cmipdiffd,extent(180,359.5,-90,90))
cmipdiffb <- crop(cmipdiffd,extent(-0.5,180,-90,90))
extent(cmipdiffa) <- c(-180,-0.5,-90,90)
cmipdiffm <- merge(cmipdiffa,cmipdiffb)
cmipdiffmd <- disaggregate(cmipdiffm,fact=10)

#add sst difference to noaa raster to find 2100 temps
temp2100 <- overlay(cmipdiffmd,noaa.mmm,fun=sum,na.rm=T,filename="cmip.rcp85.2100.tif")




###
#this time with cmip6 tos (SST)

#get all cmip6 - historical and ssp5-rcp85

#for each historical model, get mean of february SST for 1995-2015

#for each historical model, get feb means for 1985-1990 + 1993 (for comparison against NOAA). calculate model mean and median

#taylor diagram - justify the model median approach. compare historical model sst to noaa sst for 2015
