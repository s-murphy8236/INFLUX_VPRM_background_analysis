##### written by: Sam Murphy
##### purpose: read in influence functions in the raster format, 
##### read in stacked VPRM NEE in the 99x99 grid, do convolution, and 
##### calculate concentration 
##### Note: terms influence function and footprint are used interchangeably in script

# load in libraries 
library('sp')
library('raster')
library('rgdal')
library('gdalUtils')
library(lubridate)
#library(ggplot2)
library(ncdf4)
library(gtools)

######################### set month and tower to run ###########################

# year and month to run (how the tower influence function dates are formatted)
yr_month <- '2020-04'

# Tower to run (Tower_01H3,Tower_09,Tower_14)
tower <- 'Tower_14'

# write the year and month in the same format as the VPRM outputs 
yr_month_vprm <- '202004'

print(yr_month)
print(tower)
print(yr_month_vprm)

# making a temporary directory (need to do this for raster calculations)
tmpdir <- paste0('/gpfs/group/nlm136/INFLUX/SLM/temporary_dirs/convolution_temp_',tower,'_',yr_month,'/',sep='')
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir


################ read in influence function raster stacks #####################

# directory for influence function rasters 
footprint_raster_path_month <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/inf_func_rasters/',tower,'_',yr_month,'/')

# grab the names of the influence function raster stacks (length = number of hours in the month)
footprint_data_files <- list.files(path=footprint_raster_path_month, recursive=TRUE, pattern=glob2rx(paste0('*','.tif')), full.names=TRUE)
# make sure the files are in numerical order (since files are named 1,2,.. and not 001,002,.. they are pulled out of order)
footprint_data_files <- mixedsort(footprint_data_files)


################ load in prepared data and run convolution ####################

# path to the VPRM NEE raster stacks for the month 
vprm_stack_path <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/vprm_72hr_stacks/',yr_month,'/')
# get list of the VPRM NEE files for the month (length = number of hours in the month)
vprm_stack_files <- list.files(path=vprm_stack_path, recursive=TRUE, pattern=glob2rx(paste0('*',yr_month,'*.tif')), full.names=TRUE)
# force the files to be in numerical order 
vprm_stack_files <- mixedsort(vprm_stack_files)

# set directory for the month's convolution stacks to go 
convolution_stacks_path <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/convolution_stacks/',tower,'_',yr_month,'/')
dir.create(convolution_stacks_path)


# make an empty stack to put the multiplication of the NEE and influence functions
convolution_stack <- stack()

# loop to convert VPRM units and multiply with the influence functions (output is stack)
for (i in 1:length(vprm_stack_files)) { 
  
  print(paste('convolution working on hour',i,sep=' '))
  
  # read in the influence function stack for the hour of the month 
  inf_fun <- stack(footprint_data_files[i])
  
  # read in the VPRM NEE stack for the hour of the month 
  vprm_nee <- stack(vprm_stack_files[i])
  
  # convert vprm units from umol/m^2/s to mol/km^2/hr 
  vprm_nee_convert <- (vprm_nee)*3600
  
  # multiply the NEE and the influence functions, result is in stack (units?)
  convolution_stack <- inf_fun * vprm_nee_convert
  
  # save convolution stacks 
  writeRaster(convolution_stack,filename=paste0(convolution_stacks_path,yr_month,'_convolution_hour_',i,'.tif'),overwrite=TRUE)
  
  # clean up the loop 
  rm(inf_fun,vprm_nee,vprm_nee_convert,convolution_stack)
  
}

# make an empty vector to store the hourly co2 concentration values 
conv_hrly <- rep(NA,length(vprm_stack_files))

# loop to calculate the hourly co2 concentration 
for (i in 1:length(vprm_stack_files)){ 
  
  print(paste('calculating mole frac for hour',i,sep=' '))
  
  # read in the convolution stack for hour i 
  conv_stack_in <- stack(paste0(convolution_stacks_path,yr_month,'_convolution_hour_',i,'.tif'))
  
  # sum over all the raster layers to create one raster layer 
  conv_hr_sum_raster <- sum(conv_stack_in, na.rm=TRUE)
  
  # sum over all cells in the single raster layer that represents the sum 
  conv_hr_sum <- cellStats(conv_hr_sum_raster,stat='sum')
  
  # assign the value of the sum to the vector for hour i of the month (units: ppb)
  conv_hrly[i] <- conv_hr_sum
  
  # clean up 
  rm(conv_hr_sum,conv_hr_sum_raster,conv_stack_in)
  
}


################### unit conversion and saving outputs ########################

# convert the hourly concentration vector to a data frame 
conv_hrly_df <- as.data.frame(conv_hrly)

# convert units from ppb to ppm 
conv_hrly_df <- conv_hrly_df/1000

# save out the data as a csv
write.csv(conv_hrly_df,paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/convolution_molefrac_outputs/biogenic_co2_molefrac_',tower,'_',yr_month,'.csv'))

print('------------------- done :) -------------------')
