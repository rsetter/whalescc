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
  febbase <- calc(febstack,fun=mean)
  #do focal statistics to fill in coastal gaps
  febbasef <- focal(febbase,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb_focal.tif"))
  
  #extract just the augs
  augstack <- stack(r8501[[augs[1:17]]])
  #calculate mean aug temp for 1985-2001 (base period) for each model. save
  augbase <- calc(augstack,fun=mean)
  #do focal statistics to fill in coastal gaps
  augbasef <- focal(augbase,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                    filename=paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_aug_focal.tif"))
}



###ssp245

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/regrid/ssp245"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0", "GFDL-CM4", "GFDL-ESM4", "MIROC6", "MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")

febs <- c(2+12*(0:85)) #86 years between 2015-2100 
augs <- c(8+12*(0:85))

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


#calculate multi-model median delta
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/ssp245",pattern="feb_focal",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("feb245.", i)
  r <- brick(febmodelfiles[i],varname="tos")
  assign(object,r)
}
febmodels <- lapply(paste0('feb245.',1:length(febmodelfiles)),get)

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

deltaaug <- overlay(aug245.1,aug245.2,aug245.3,aug245.4,aug245.5,aug245.6,aug245.7,aug245.8,aug245.9,aug245.10,
        aug245.11,aug245.12,aug245.13,aug245.14,aug245.15,aug245.16,aug245.17,aug245.18,aug245.19,
        fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_aug.tif")
deltaaugd <- disaggregate(deltaaug, fact=20)
i_brick1 <- crop(deltaaugd,extent(180,359.5,-90,90))
i_brick2 <- crop(deltaaugd,extent(-0.5,180,-90,90))
extent(i_brick1) <- c(-180,-0.5,-90,90)
deltaaugr <- merge(i_brick1,i_brick2)
writeRaster(deltaaugr,filename="F:/whales/cmip6/ssp245/modelmedian_decadedelta_augr.tif")



#add ensemble median delta to esa to get downscaled sst
#for feb
febsst245 <- overlay(esafebc,deltafebr, fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245.tif")
#for aug
augsst245 <- overlay(esaaugc,deltaaugr, fun=sum,filename="F:/whales/cmip6/ssp245/AugSST_ssp245.tif")




###ssp585

#open each cmip6 model. get just the feb/august. 
parent.folder <- "F:/cmip6/tos/regrid/ssp585"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")

febs <- c(2+12*(0:85)) #86 years between 2015-2100 
augs <- c(8+12*(0:85))

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


#add each model delta to feb esa
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
foreach::foreach(i=1:length(febmodelfiles), .packages=c("raster")) %dopar% {  modeld <- disaggregate(febmodels[[i]], fact=20)
  i_brick1 <- crop(modeld,extent(180,359.5,-90,90))
  i_brick2 <- crop(modeld,extent(-0.5,180,-90,90))
  extent(i_brick1) <- c(-180,-0.5,-90,90)
  modeldm <- merge(i_brick1,i_brick2)
  overlay(esafebc,modeldm,fun=sum,filename=paste0("F:/whales/cmip6/ssp585/FebSST_ssp585_",modelnames[i],".tif"))
}
stopCluster(cl)




###
### calculate coarse grid model median sst
###
#calculate model median base sst
#for feb
febmodelfiles <- list.files("F:/whales/cmip6/hist",pattern="basesst_feb_focal",full.names=T)
for(i in 1:length(febmodelfiles)){
  object <- paste0("febbase.", i)
  r <- raster(febmodelfiles[i],varname="tos")
  assign(object,r)
}
febbasemed <- overlay(febbase.1,febbase.2,febbase.3,febbase.4,febbase.5,febbase.6,febbase.7,febbase.8,febbase.9,febbase.10,
                      febbase.11,febbase.12,febbase.13,febbase.14,febbase.15,febbase.16,febbase.17,febbase.18,febbase.19,febbase.20,
                      fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/hist/modelmedian_basesst_feb.tif")

#for aug
augmodelfiles <- list.files("F:/whales/cmip6/hist",pattern="basesst_aug_focal",full.names=T)
for(i in 1:length(augmodelfiles)){
  object <- paste0("augbase.", i)
  r <- raster(augmodelfiles[i],varname="tos")
  assign(object,r)
}
augbasemed <- overlay(augbase.1,augbase.2,augbase.3,augbase.4,augbase.5,augbase.6,augbase.7,augbase.8,augbase.9,augbase.10,
                      augbase.11,augbase.12,augbase.13,augbase.14,augbase.15,augbase.16,augbase.17,augbase.18,augbase.19,augbase.20,
                      fun=function(x){median(x,na.rm=T)},filename="F:/whales/cmip6/hist/modelmedian_basesst_aug.tif")

#open deltas
deltaaug585 <- brick("F:/whales/cmip6/ssp585/modelmedian_decadedelta_aug.tif")
deltafeb585 <- brick("F:/whales/cmip6/ssp585/modelmedian_decadedelta_feb.tif")
deltaaug245 <- brick("F:/whales/cmip6/ssp245/modelmedian_decadedelta_aug.tif")
deltafeb245 <- brick("F:/whales/cmip6/ssp245/modelmedian_decadedelta_feb.tif")

#add delta to base sst to get coarse grid projections
#for feb
febsst585c <- overlay(febbasemed,deltafeb585, fun=sum,filename="F:/whales/cmip6/ssp585/FebSST_ssp585_coarse.tif",overwrite=T)
febsst245c <- overlay(febbasemed,deltafeb245, fun=sum,filename="F:/whales/cmip6/ssp245/FebSST_ssp245_coarse.tif",overwrite=T)

#for aug
augsst585c <- overlay(augbasemed,deltaaug585, fun=sum,filename="F:/whales/cmip6/ssp585/AugSST_ssp585_coarse.tif",overwrite=T)
augsst245c <- overlay(augbasemed,deltaaug245, fun=sum,filename="F:/whales/cmip6/ssp245/AugSST_ssp245_coarse.tif",overwrite=T)


#rotate it
febsst585c <- brick("F:/whales/cmip6/ssp585/FebSST_ssp585_coarse.tif")
febsst585r <- rotate(febsst585c,filename="F:/whales/cmip6/ssp585/FebSST_ssp585_coarse180.tif", overwrite=T)

febsst245c <- brick("F:/whales/cmip6/ssp245/FebSST_ssp245_coarse.tif")
febsst245r <- rotate(febsst245c,filename="F:/whales/cmip6/ssp245/FebSST_ssp245_coarse180.tif", overwrite=T)

augsst585c <- brick("F:/whales/cmip6/ssp585/AugSST_ssp585_coarse.tif")
augsst585r <- rotate(augsst585c,filename="F:/whales/cmip6/ssp585/AugSST_ssp585_coarse180.tif", overwrite=T)

augsst245c <- brick("F:/whales/cmip6/ssp245/AugSST_ssp245_coarse.tif")
augsst245r <- rotate(augsst245c,filename="F:/whales/cmip6/ssp245/AugSST_ssp245_coarse180.tif", overwrite=T)






### calculate difference between coarse and downscaled for feb ssp585 year 2100
setwd("/Volumes/BlackSea/Whales with Hannah/whales")
coarse <- brick("/Volumes/BlackSea/Whales with Hannah/whales/cmip6/ssp585/FebSST_ssp585_coarse180.tif")
downscaled <- brick("/Volumes/BlackSea/Whales with Hannah/whales/cmip6/ssp585/FebSST_ssp585.tif")
#make coarse dataset same extent and resolution as downscaled
extent(coarse) <- extent(downscaled)
coarsed <- disaggregate(coarse[[9]],fact=20)

difference <- downscaled[[9]] - coarsed
writeRaster(difference, filename="/Volumes/BlackSea/Whales with Hannah/whales/cmip6/ssp585/FebSST_ssp585_difference.tif")






#### taylor diagram
###
#empirical & hist taylor diagrams
#calculate accuracy: Taylor diagrams

#use overlapping years with historical scenarios (1985 to 2001)
esafebc <- raster("F:/whales/ESA_SST/ESA_feb_meanC.tif")
esafebcdf <- as.data.frame(esafebc,xy=T)



#open each cmip6 model. get just the feb base value (mean value for 1985 to 2001)
parent.folder <- "/Volumes/BlackSea/Whales with Hannah/whales/cmip6/tos/regrid/hist"
modelnames <- c("ACCESS-CM2","ACCESS-ESM1-5","CAMS-CSM1-0","CanESM5","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","FGOALS-f3-L",
                "FGOALS-g3","FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","IPSL-CM6A-LR","MIROC6","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM")
allmodbase <- raster()
for(i in 1:length(modelnames)){
  r <- raster::raster(paste0("F:/whales/cmip6/hist/",modelnames[i], "_basesst_feb_focal.tif"))
  ro <- raster::rotate(r)
  rod <- disaggregate(ro,fact=20)
  extent(rod) <- extent(esafebc)
  rdf <- as.data.frame(rod,xy=T)
  
  a <- paste("mod", i, sep = "")
  assign(a,rdf)
  
  allmodbase <- stack(allmodbase,r)
}

#get the multi-model median 
modmedian <- raster::calc(allmodbase,fun=median,na.rm=T)
ro <- raster::rotate(modmedian)
rod <- disaggregate(ro,fact=20)
extent(rod) <- extent(esafebc)
modmeddf <- as.data.frame(rod,xy=T)

#calculate the multi-model mean
modmean <- raster::calc(allmodbase,fun=mean,na.rm=T)
ro <- raster::rotate(modmean)
rod <- disaggregate(ro,fact=20)
extent(rod) <- extent(esafebc)
modmeandf <- as.data.frame(rod,xy=T)


##join all (by x and y) into one dataframe for comparison
library(dplyr)
alldata<- inner_join(esafebcdf,modmeddf,by=c("x","y"))
alldata<- inner_join(alldata,modmeandf,by=c("x","y"))
alldata<- inner_join(alldata,mod1,by=c("x","y"))
alldata<- inner_join(alldata,mod2,by=c("x","y"))
alldata<- inner_join(alldata,mod3,by=c("x","y"))
alldata<- inner_join(alldata,mod4,by=c("x","y"))
alldata<- inner_join(alldata,mod5,by=c("x","y"))
alldata<- inner_join(alldata,mod6,by=c("x","y"))
alldata<- inner_join(alldata,mod7,by=c("x","y"))
alldata<- inner_join(alldata,mod8,by=c("x","y"))
alldata<- inner_join(alldata,mod9,by=c("x","y"))
alldata<- inner_join(alldata,mod10,by=c("x","y"))
alldata<- inner_join(alldata,mod11,by=c("x","y"))
alldata<- inner_join(alldata,mod12,by=c("x","y"))
alldata<- inner_join(alldata,mod13,by=c("x","y"))
alldata<- inner_join(alldata,mod14,by=c("x","y"))
alldata<- inner_join(alldata,mod15,by=c("x","y"))
alldata<- inner_join(alldata,mod16,by=c("x","y"))
alldata<- inner_join(alldata,mod17,by=c("x","y"))
alldata<- inner_join(alldata,mod18,by=c("x","y"))
alldata<- inner_join(alldata,mod19,by=c("x","y"))
alldata<- inner_join(alldata,mod20,by=c("x","y"))

names(alldata) <- c("x","y","esa","medianmod","meanmod","model1","model2","model3","model4","model5","model6","model7","model8","model9","model10",
                    "model11","model12","model13","model14","model15","model16","model17","model18","model19","model20")
write.csv(alldata,"F:/whales/cmip6/hist/taylordiagram.csv",row.names=F)

#make the diagram
library(plotrix)
oldpar<-taylor.diagram(alldata$esa,alldata$esa,ref.sd=T,normalize=TRUE,sd.arcs=TRUE,pcex=1.3,pch=19,col="blue",gamma.col="blue",
                       main="Taylor Diagram - SST")
oldpar<-taylor.diagram(alldata$esa,alldata$model1,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model2,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model3,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model4,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model5,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model6,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model7,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model8,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model9,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model10,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model11,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model12,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model13,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model14,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model15,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model16,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model17,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model18,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$model20,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1,pch=19,col="black")
oldpar<-taylor.diagram(alldata$esa,alldata$meanmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.3,col="red",pch=19)
oldpar<-taylor.diagram(alldata$esa,alldata$medianmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.3,col="#4cb84c",pch=19)







