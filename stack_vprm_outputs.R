##### written by: Sam Murphy
##### purpose: read in VPRM NEE in the 99x99 grid, stack VPRM files to make 
##### the convolution easier 

# load in libraries 
library('sp')
library('raster')
#library('rgdal')
#library('gdalUtils')
#library('lubridate')
#library(ggplot2)
#library(ncdf4)
library(gtools)


# year and month to run (how the tower influence function dates are formatted)
yr_month <- '2020-04'

# write the year and month in the same format as the VPRM outputs 
yr_month_vprm <- '202004'

# need the month previous to the one you're running for the VPRM to grab the 72 hours before 
yr_month_vprm_pre <- '202003'

print(yr_month)
print(yr_month_vprm)
print(yr_month_vprm_pre)

# making a temporary directory (need to do this for raster reprojections)
tmpdir <- paste0('/gpfs/group/nlm136/INFLUX/SLM/temporary_dirs/vprm_stack_temp_',yr_month,'/',sep='')
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir



############## load in regridded VRPM NEE and stack them ####################

# path to the regridded VPRM NEE (data will now be 99x99) for the month of interest 
VPRM_regrid_month <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/vprm_regrid/',yr_month_vprm,'/')

# path to the regreddied VPRM NEE for the month before the month of interest (needed to be able 
# to grab the 72 hours previous to each hour early in the month)
VPRM_regrid_month_pre <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/vprm_regrid/',yr_month_vprm_pre,'/')

# get a list of filenames for the month of interest 
VPRM_regrid_month_files <- list.files(path=VPRM_regrid_month, recursive=TRUE, pattern=glob2rx(paste0('*VPRM_NEE_99x99_*',yr_month_vprm,'*.tif')), full.names=TRUE)

# get a list of filenames for the month previous 
VPRM_regrid_month_pre_files <- list.files(path=VPRM_regrid_month_pre, recursive=TRUE, pattern=glob2rx(paste0('*VPRM_NEE_99x99_*',yr_month_vprm_pre,'*.tif')), full.names=TRUE)

# grab only the last 72 hours of the previous month (actually the last 70 hours, to account for hour 0 of the month of interest)
# VPRM_regrid_month_pre_files_subset <- VPRM_regrid_month_pre_files[(length(VPRM_regrid_month_pre_files)-71):length(VPRM_regrid_month_pre_files)]
VPRM_regrid_month_pre_files_subset <- VPRM_regrid_month_pre_files[(length(VPRM_regrid_month_pre_files)-70):length(VPRM_regrid_month_pre_files)]

# combine the end (72hrs) of previous month with the rest of the month of interest 
VPRM_regrid <- c(VPRM_regrid_month_pre_files_subset,VPRM_regrid_month_files)

# set and make a directory to the store the 72 hour stacks for the month 
path_72hr_stack <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/vprm_72hr_stacks/',yr_month,'/')
dir.create(path_72hr_stack)

# make an empty stack to store the 72 hour rasters for each hour 
vprm_72hr_stack <- stack()

# loop for stacking the NEE of 72 hours previous to the hour of interest (this 
# is to match the times of the influence function stacks)
n=1 # this is for naming the files (hour 1 to hour # of the end of month)
#for (i in 73:length(VPRM_regrid)){ # do i need one hour into the next month? hour zero of the next month??
for (i in 73:(length(VPRM_regrid)+1)){  
  # note: start at 73, this starts at the first hour of the month of interest 
  
  print(paste('stacking hour number',n,'of',length(73:(length(VPRM_regrid)+1)),sep=' '))
  
  # grab the file names for the 72 hours before hour i 
  vprm_72hr_files <- VPRM_regrid[(i-72):(i-1)]
  #vprm_72hr_files <- VPRM_regrid[(i-73):(i)]
  
  # need to reverse the order of the VPRM files to match that of the influence functions
  # (goes from hour i-1 to hour i-72)
  vprm_72hr_files <- rev(vprm_72hr_files)
  
  # load in the 72 hours as a raster stack 
  vprm_72hr_stack <- stack(vprm_72hr_files)
  
  # write the raster stack with the hour of the month (n) in the name 
  writeRaster(vprm_72hr_stack,filename=paste0(path_72hr_stack,yr_month,'_NEE_stack_hour_',n,'.tif'),overwrite=TRUE)
  
  # clean up variables to use in next loop 
  rm(vprm_72hr_stack,vprm_72hr_files)
  
  # add one to make n the next hour in the month 
  n=n+1
  
}

print('------------------- done :) -------------------')
