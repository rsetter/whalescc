library(raster)
library(gdalUtils)
library(gdal)
library(foreach)
library(doParallel)


setwd("F:/whales")


#ESA SST
#in kelvin


parent.folder <- "F:/whales/ESA_SST"
allfolders <- list.dirs(parent.folder)

#open the feb files
#make all files into a brick 
febfolders <- allfolders[c(3*(1:17))]
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
foreach::foreach(i=1:length(febfolders), .packages=c("raster")) %dopar% {
#for(i in 1:length(febfolders)){
  febyrfiles <- list.files(febfolders[[i]],full.names=T)
  februaries <- raster()
  for(j in 1:length(febyrfiles)){
     r <- raster(febyrfiles[[j]],varname="analysed_sst")
     februaries <- stack(februaries,r)
  }
  #get mean feb SST per year
  calc(februaries,fun=function(x){mean(x,na.rm=T)},filename=paste0("F:/whales/ESA_SST/ESA_feb_mean_",1984+i,".tif"))
}
stopCluster(cl)

#get the 1985-2001 mean february SST 
febmeanyrfiles <- list.files(parent.folder,pattern="feb_mean",full.names=T)
febyears <- raster()
for(i in 1:length(febmeanyrfiles)){
  r <- raster(febmeanyrfiles[[i]])
  febyears <- stack(febyears,r)
}

esafeb <- calc(febyears,fun=function(x){mean(x,na.rm=T)},filename="F:/whales/ESA_SST/ESA_feb_mean.tif")
esafebc <- esafeb-273.15
writeRaster(esafebc, filename="F:/whales/ESA_SST/ESA_feb_meanC.tif")

#open the aug files
#make all files into a brick 
augfolders <- allfolders[c(3*(1:17)+1)]
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
foreach::foreach(i=1:length(augfolders), .packages=c("raster")) %dopar% {
#for(i in 1:length(augfolders)){
  augyrfiles <- list.files(augfolders[[i]],full.names=T)
  augusts <- raster()
  for(j in 1:length(augyrfiles)){
    r <- raster(augyrfiles[[j]],varname="analysed_sst")
    augusts <- stack(augusts,r)
  }
  #get mean aug SST per year
  calc(augusts,fun=function(x){mean(x,na.rm=T)},filename=paste0("F:/whales/ESA_SST/ESA_aug_mean_",1984+i,".tif"))
}
stopCluster(cl)

#get the 1985-2001 mean august SST 
augmeanyrfiles <- list.files(parent.folder,pattern="aug_mean",full.names=T)
augyears <- raster()
for(i in 1:length(augmeanyrfiles)){
  r <- raster(augmeanyrfiles[[i]])
  augyears <- stack(augyears,r)
}

#get the 1985-2001 mean august SST 
esaaug <- calc(augyears,fun=function(x){mean(x,na.rm=T)},filename="F:/whales/ESA_SST/ESA_aug_mean.tif")
esaaugc <- esaaug-273.15
writeRaster(esaaugc, filename="F:/whales/ESA_SST/ESA_aug_meanC.tif")









###
#pre-process all cmip6 datasets - regrid from curvilinear to rectilinear for easy comparison with AVHRR
#use CDO in Terminal https://code.mpimet.mpg.de/projects/cdo/embedded/index.html#x1-6190002.12.1
# see info: cdo sinfon in.nc
# concatenate files: ncrcat 85.nc 86.nc 87.nc 88.nc 89.nc 8589.nc
# remap with bilinear interpolation: cdo remapbil,r360x180 in.nc out.nc

#trouble with IPSL ssp245 & ssp545







###hist 

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/regrid/hist"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
febs <- c(2+12*(0:16)) #17 years between 1985-2001
augs <- c(8+12*(0:16))


for(i in 1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #models are concatenated in one file. open full dataset file
  r <- brick(modelfiles[1],varname="tos")
  #keep just the years 1985-2001
  r8501 <- dropLayer(r,c(1:1620,(nlayers(r)-155):nlayers(r)))
  #extract just the febs
  febstack <- stack(r8501[[febs[1:17]]])
  #calculate mean feb temp for 1985-2001 (base period) for each model. save
  febbase <- calc(febstack,fun=mean)#,filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb.tif"))
  #do focal statistics to fill in coastal gaps
  febbasef <- focal(febbase,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb_focal.tif"))
  
  #extract just the augs
  augstack <- stack(r8501[[augs[1:17]]])
  #calculate mean aug temp for 1985-2001 (base period) for each model. save
  augbase <- calc(augstack,fun=mean)#,filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug.tif"))
  #do focal statistics to fill in coastal gaps
  augbasef <- focal(augbase,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                    filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug_focal.tif"))
}



###ssp245

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/regrid/ssp245"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0", "GFDL-CM4", "GFDL-ESM4", "MIROC6", "MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
#skip "IPSL-CM6A-LR" - could not regrid

febs <- c(2+12*(0:85)) #86 years between 2015-2100 
augs <- c(8+12*(0:85))

#run CAMS-CSM1-0 separately - only goes up to 2099, not 2100. use febstack <- stack(r[[febs[1:85]]]) and augstack <- stack(r[[augs[1:86]]])
for(i in  1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #make single brick with all model years 2015-2100
  r <- brick(modelfiles[1],varname="tos")
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
  febbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb_focal.tif"))
  febdeca <- raster()
  for(a in 1:nlayers(febdecadal)){
    f <- focal(febdecadal[[a]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    febdeca <- stack(febdeca,f)
  }
  overlay(febdeca,febbase, fun=function(x,y,na.rm=T){return(x-y)},
                      filename=paste0("F:/whales/cmip6/ssp245/",modelnames[i],"_decadedelta_feb_focal.tif"))
  
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
  augbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug_focal.tif"))
  augdeca <- raster()
  for(b in 1:nlayers(augdecadal)){
    f <- focal(augdecadal[[b]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    augdeca <- stack(augdeca,f)
  }
  overlay(augdeca,augbase, fun=function(x,y,na.rm=T){return(x-y)},
                      filename=paste0("F:/whales/cmip6/ssp245/",modelnames[i],"_decadedelta_aug_focal.tif"))
}


#calculate model median delta
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/ssp245",pattern="feb_focal",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("feb245.", i)
  r <- brick(febmodelfiles[i],varname="tos")
  assign(object,r)
}
febmodels <- lapply(paste0('feb245.',1:length(febmodelfiles)),get)
#find minimum & maximum delta models
febdataframe <- data.frame()
for(i in 1:19){
  crop <- crop(febmodels[[i]][[9]],extent(180,230,2,37))
  meandelta <- cellStats(crop,stat=mean)
  x <- data.frame(model = i, mean = meandelta)
  febdataframe <- rbind(x,febdataframe)
}
#min = model3
#max = model11

deltafeb <- overlay(feb245.1,feb245.2,feb245.3,feb245.4,feb245.5,feb245.6,feb245.7,feb245.8,feb245.9,feb245.10,
        feb245.11,feb245.12,feb245.13,feb245.14,feb245.15,feb245.16,feb245.17,feb245.18,feb245.19,
        fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_feb.tif")
deltafebd <- disaggregate(deltafeb, fact=20)
i_brick1 <- crop(deltafebd,extent(180,359.5,-90,90))
i_brick2 <- crop(deltafebd,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
deltafebr <- merge(i_brick1,i_brick2)
writeRaster(deltafebr,filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_febr.tif")


#for aug
augmodelfiles <- list.files("F:/whales/cmip6/ssp245",pattern="aug_focal",full.names=T)
for(i in 1:length(augmodelfiles)){
  object <- paste0("aug245.", i)
  r <- brick(augmodelfiles[i])
  assign(object,r)
}
augmodels <- lapply(paste0('aug245.',1:length(augmodelfiles)),get)
#find minimum & maximum delta models
augdataframe <- data.frame()
for(i in 1:19){
  crop <- crop(febmodels[[i]][[9]],extent(180,230,2,37))
  meandelta <- cellStats(crop,stat=mean)
  x <- data.frame(model = i, mean = meandelta)
  augdataframe <- rbind(x,augdataframe)
}
#min = model11
#max = model16

deltaaug <- overlay(aug245.1,aug245.2,aug245.3,aug245.4,aug245.5,aug245.6,aug245.7,aug245.8,aug245.9,aug245.10,
        aug245.11,aug245.12,aug245.13,aug245.14,aug245.15,aug245.16,aug245.17,aug245.18,aug245.19,
        fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_aug.tif")
deltaaugd <- disaggregate(deltaaug, fact=20)
i_brick1 <- crop(deltaaugd,extent(180,359.5,-90,90))
i_brick2 <- crop(deltaaugd,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
deltaaugr <- merge(i_brick1,i_brick2)
writeRaster(deltaaugr,filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_augr.tif")




#add delta to esa to get downscaled esa sst
#for feb
febsst245 <- overlay(esafebc,deltafebr, fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245.tif")
#for aug
augsst245 <- overlay(esaaugc,deltaaugr, fun=sum,filename="F:/whales/cmip6/ssp245/AugSST_ssp245.tif")

#add min model (model3) delta to feb esa 
feb245.3d <- disaggregate(feb245.3, fact=20)
i_brick1 <- crop(feb245.3d,extent(180,359.5,-90,90))
i_brick2 <- crop(feb245.3d,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
feb245.3dm <- merge(i_brick1,i_brick2)
febmin245 <- overlay(esafebc,feb245.3dm,fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245_model3_min.tif")

#add max model (model11) delta to feb esa
feb245.11d <- disaggregate(feb245.11, fact=20)
i_brick1 <- crop(feb245.11d,extent(180,359.5,-90,90))
i_brick2 <- crop(feb245.11d,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
feb245.11dm <- merge(i_brick1,i_brick2)
febmax245 <- overlay(esafebc,feb245.11dm,fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245_model11_max.tif")



###ssp585

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/regrid/ssp585"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
#skip "IPSL-CM6A-LR" - could not regrid

febs <- c(2+12*(0:85)) #86 years between 2015-2100 
augs <- c(8+12*(0:85))

#run CAMS-CSM1-0 separately - only goes up to 2099, not 2100. use febstack <- stack(r[[febs[1:85]]]) and augstack <- stack(r[[augs[1:86]]])
for(i in 1:length(modelnames)){
  modelfiles <- list.files(parent.folder, pattern=modelnames[i],full.names=T)
  #make single brick with all model years 2015-2100
  r <- brick(modelfiles[1],varname="tos")
  #extract just the febs. this will be feb temps 2015-2100 (rasterbrick with 86 layers)
  febstack <- stack(r[[febs[1:86]]])
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
  febbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb_focal.tif"))
  febdeca <- raster()
  for(a in 1:nlayers(febdecadal)){
    f <- focal(febdecadal[[a]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    febdeca<-stack(febdeca,f)
  }
  overlay(febdeca,febbase, fun=function(x,y,na.rm=T){return(x-y)},
          filename=paste0("F:/whales/cmip6/ssp585/",modelnames[i],"_decadedelta_feb_focal.tif"),overwrite=T)
  
  #extract just the augs. this will be aug temps 2015-2100 (rasterbrick with 86 layers)
  augstack <- stack(r[[augs[1:86]]])
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
  augbase <- raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug_focal.tif"))
  augdeca <- raster()
  for(b in 1:nlayers(augdecadal)){
    f <- focal(augdecadal[[b]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    augdeca <- stack(augdeca,f)
  }
  overlay(augdeca,augbase, fun=function(x,y,na.rm=T){return(x-y)},
          filename=paste0("F:/whales/cmip6/ssp585/",modelnames[i],"_decadedelta_aug_focal.tif"),overwrite=T)
}

#calculate model median delta
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/ssp585",pattern="feb_focal",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("feb585.", i)
  r <- brick(febmodelfiles[i])
  assign(object,r)
}
febmodels <- lapply(paste0('feb585.',1:length(febmodelfiles)),get)
#find minimum & maximum delta models
febdataframe <- data.frame()
for(i in 1:18){
  crop <- crop(febmodels[[i]][[9]],extent(180,230,2,37))
  meandelta <- cellStats(crop,stat=mean)
  x <- data.frame(model = i, mean = meandelta)
  febdataframe <- rbind(x,febdataframe)
}
#min = model3
#max = model4

deltafeb <- overlay(feb585.1,feb585.2,feb585.3,feb585.4,feb585.5,feb585.6,feb585.7,feb585.8,feb585.9,feb585.10,
                    feb585.11,feb585.12,feb585.13,feb585.14,feb585.15,feb585.16,feb585.17,feb585.18,
                    fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/ssp585/modelmedian_decadedelta_feb.tif",overwrite=T)
deltafebd <- disaggregate(deltafeb, fact=20)
i_brick1 <- crop(deltafebd,extent(180,359.5,-90,90))
i_brick2 <- crop(deltafebd,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
deltafebr <- merge(i_brick1,i_brick2)
writeRaster(deltafebr,filename="F:/whales/cmip6/ssp585/modelmedian_decadedelta_febr.tif",overwrite=T)

#for aug
augmodelfiles <- list.files("F:/whales/cmip6/ssp585",pattern="aug_focal",full.names=T)
for(i in 1:length(augmodelfiles)){
  object <- paste0("aug585.", i)
  r <- brick(augmodelfiles[i])
  assign(object,r)
}
augmodels <- lapply(paste0('aug585.',1:length(febmodelfiles)),get)
# #find minimum & maximum delta models
# augdataframe <- data.frame()
# for(i in 1:18){
#   crop <- crop(augmodels[[i]][[9]],extent(200,210,17,24))
#   meandelta <- cellStats(crop,stat=mean)
#   x <- data.frame(model = i, mean = meandelta)
#   augdataframe <- rbind(x,augdataframe)
# }
# #min = model
# #max = model
deltaaug <- overlay(aug585.1,aug585.2,aug585.3,aug585.4,aug585.5,aug585.6,aug585.7,aug585.8,aug585.9,aug585.10,
                    aug585.11,aug585.12,aug585.13,aug585.14,aug585.15,aug585.16,aug585.17,aug585.18,
                    fun=function(x){median(x,na.rm=T)},filename=paste0("F:/whales/cmip6/ssp585/modelmedian_decadedelta_aug.tif"),overwrite=T)
deltaaugd <- disaggregate(deltaaug, fact=20)
i_brick1 <- crop(deltaaugd,extent(180,359.5,-90,90))
i_brick2 <- crop(deltaaugd,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
deltaaugr <- merge(i_brick1,i_brick2)
writeRaster(deltaaugr,filename="F:/whales/cmip6/ssp585/modelmedian_decadedelta_augr.tif",overwrite=T)


#add delta to esa to get downscaled esa sst
#for feb
febsst585 <- overlay(esafebc,deltafebr, fun=sum,filename="F:/whales/cmip6/ssp585/FebSST_ssp585.tif",overwrite=T)
#for aug
augsst585 <- overlay(esaaugc,deltaaugr, fun=sum,filename="F:/whales/cmip6/ssp585/AugSST_ssp585.tif",overwrite=T)

#add min model (model3) delta to feb esa 
feb585.3d <- disaggregate(feb585.3, fact=20)
i_brick1 <- crop(feb585.3d,extent(180,359.5,-90,90))
i_brick2 <- crop(feb585.3d,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
feb585.3dm <- merge(i_brick1,i_brick2)
febmin585 <- overlay(esafebc,feb585.3dm,fun=sum,filename="F:/whales/cmip6/ssp585/FebSST_ssp585_model3_min.tif",overwrite=T)

#add max model (model4) delta to feb esa
feb585.4d <- disaggregate(feb585.4, fact=20)
i_brick1 <- crop(feb585.4d,extent(180,359.5,-90,90))
i_brick2 <- crop(feb585.4d,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
feb585.4dm <- merge(i_brick1,i_brick2)
febmax585 <- overlay(esafebc,feb585.4dm,fun=sum,filename="F:/whales/cmip6/ssp585/FebSST_ssp585_model4_max.tif",overwrite=T)



























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



###avhrr


#open avhrr dataset. sst climatology for 1985-2001
#convert HDF files to TIF
# sds <- get_subdatasets("F:/whales/avhrr/month02_combined.hdf")
# gdal_translate(sds[1], dst_dataset = "F:/whales/avhrr/month02_combined_band0.tif")
# 
# avhrrfeb <- raster("F:/whales/avhrr/month02_combined_band5.tif",
# crs="+init=epsg:4267 +proj=longlat +ellps=clrk66 +datum=NAD27
# +no_defs")
# 
# sds <- get_subdatasets("F:/whales/avhrr/month08_combined.hdf")
# gdal_translate(sds[1], dst_dataset = "F:/whales/avhrr/month08_combined_band0.tif")
# 
# avhrraug <- raster("F:/whales/avhrr/month08_combined_band0.tif")


# inpath <- "D:/climatologydata/r"
# setwd(inpath)
# filenames <- list.files(,pattern = ".hdf$", full.names = F)
# for(filename in filenames)
# {
#   sds <- get_subdatasets(filename)
#   gdal_translate(sds[1],dst_dataset = paste0(substr(filename, 1, nchar(filename)-4),".tif"))
# }

##########