##### written by: Sam Murphy
##### purpose: read in VPRM outputs, regridd them to match the wrf outputs (i.e., 
##### convert the 266x354 VPRM NEE layers to 99x99 to match the influence function 
##### grids)

# load in libraries 
library('sp')
library('raster')
library('rgdal')
library('gdalUtils')
library(lubridate)
#library(ggplot2)
library(ncdf4)
library(gtools)

# write the year and month in the same format as the VPRM outputs 
yr_month_vprm <- '202004'
print(yr_month_vprm)

# making a temporary directory (need to do this for raster reprojections)
tmpdir <- paste0('/gpfs/group/nlm136/INFLUX/SLM/temporary_dirs/regrid_vprm_temp_',yr_month_vprm,'/',sep='')
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir

# path to VPRM outputs 
VPRM_path <- paste0('/gpfs/group/nlm136/INFLUX/SLM/VPRM_gridded_run_outputs/',yr_month_vprm,'/')

# get a list of the NEE files for the month - recursive = false that way we dont read in all the individual PFT files 
VPRM_files_month <- list.files(path=VPRM_path, recursive=FALSE, pattern=glob2rx(paste0('*NEE*',yr_month_vprm,'*.tif')), full.names=TRUE)

# set path to write the output 
VPRM_regrid_out <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/vprm_regrid/',yr_month_vprm,'/')
dir.create(VPRM_regrid_out)

# grab an influence function raster file to match the resolution to 
match_res <- raster('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/inf_func_rasters/Tower_09_2018-08/2018-08_hour_1.tif')

# loop over all hours in the month to regrid VPRM to 99x99 
for (i in 1:length(VPRM_files_month)) {
  
  print(paste('regridding layer',i,'of',length(VPRM_files_month),sep=' '))
  
  # read in hour i of the month 
  VPRM_layer <- raster(VPRM_files_month[i])
  
  # remove the factor of 1000 from the VPRM layer 
  VPRM_layer <- VPRM_layer/1000
  
  # aggregate to as close to 99x99 allowed (needs to be a factor of the rows/columns)
  # VPRM_agg will be 118x133 
  VPRM_agg <- aggregate(VPRM_layer,c(3,2),fun=mean,na.rm=TRUE)  
  
  # resample using bilinear interpolation to get the rest of the way to 99x99
  VPRM_disagg <- resample(VPRM_agg, match_res, method = "bilinear") 
  
  # write the 99x99 VPRM raster (still in WG1984 coordinates)
  writeRaster(VPRM_disagg,paste0(VPRM_regrid_out,'VPRM_NEE_99x99_',substr(VPRM_files_month[i],77,86),'.tif'),overwrite=TRUE)
  
  # clean up the loop 
  rm(VPRM_layer,VPRM_disagg, VPRM_agg)
  
  
}

################ check the averages and sums - for testing only ###############
i=1 

# read in hour i of the month 
VPRM_layer <- raster(VPRM_files_month[i])

# remove the factor of 1000 from the VPRM layer 
VPRM_layer <- VPRM_layer/1000

# check mean on input VPRM 
cellStats(VPRM_layer,stat=mean,na.rm=TRUE)
cellStats(VPRM_layer,stat=mean)

# check sum on input VPRM 
cellStats(VPRM_layer,stat=sum,na.rm=TRUE)
cellStats(VPRM_layer,stat=sum)

# aggregate to as close to 99x99 allowed (needs to be a factor of the rows/columns)
# VPRM_agg will be 118x133 
VPRM_agg <- aggregate(VPRM_layer,c(3,2),fun=mean,na.rm=TRUE)  

# check mean on aggregated VPRM 
cellStats(VPRM_agg,stat=mean,na.rm=TRUE)
cellStats(VPRM_agg,stat=mean)

# resample using bilinear interpolation to get the rest of the way to 99x99
VPRM_disagg <- resample(VPRM_agg, match_res, method = "bilinear") 

# check mean on resampled VPRM 
cellStats(VPRM_disagg,stat=mean,na.rm=TRUE)
cellStats(VPRM_disagg,stat=mean)

# check sum on resampled VPRM 
cellStats(VPRM_disagg,stat=sum,na.rm=TRUE)
cellStats(VPRM_disagg,stat=sum)
