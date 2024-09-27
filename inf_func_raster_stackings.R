####### script written by: Sam Murphy
####### purpose: read in the tower influence functions created by Zach from the 
####### netCDF files, turn them into rasters in WGS1984 coordinates (to match
####### the VPRM NEE outputs), and stack them to be used in the VPRM-influence 
####### function convolutions 


# load in libraries 
library('sp')
library('raster')
library('rgdal')
library('gdalUtils')
library(lubridate)
#library(ggplot2)
library(ncdf4)
library(gtools)


# year and month to run (how the tower influence function dates are formatted)
yr_month <- '2020-04'
print(yr_month)

# Tower to run ('Tower_01H3')
tower <- 'Tower_01H3'
print(tower)

# making a temporary directory (need to do this for raster reprojections)
tmpdir <- paste0('/gpfs/group/nlm136/INFLUX/SLM/temporary_dirs/inf_fun_raster_temp_',tower,'_',yr_month,'/',sep='')
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir



################# read in influence function from netCDF file ##################

# path to influence function file 
inf_func_nc_path<- paste0('/gpfs/group/nlm136/INFLUX/LPDM/Footprints/',tower,'/',yr_month,'_d02.nc')

# read in monthly netCDF file 
inf_func_nc <- nc_open(inf_func_nc_path)

# extract the footprint from the netCDF file 
footprint <- ncvar_get(inf_func_nc,varid = 'footprint',verbose = TRUE)

# get the dimensions of the footprint (should be 4d - 99 x 99 x #hours in month x 72)
dim(footprint)



##################### convert footprint to rasters ############################
# note: need to convert footprints to rasters in the same coordinate system as 
# the VPRM files (lat/lon WGS1984) 

# make directory to store the influence function raster stacks
footprint_raster_path_month <- paste0('/gpfs/group/nlm136/INFLUX/SLM/conc_tower_vprm_convolutions/inf_func_rasters/',tower,'_',yr_month,'/')
dir.create(footprint_raster_path_month)

# get the number of hours in the month 
month_hours <- length(footprint[1,1,,1])

# this loop loads in netCDF files, and saves out raster stacks for each of the hours 
# in the month
for (i in 1:month_hours){ # for all hours in the month
  
  # make an empty raster stack 
  footprint_hrly_stack <- stack() 
  
  # make a 3d array with the 72 "layers" for the hour i (dimensions = 99x99x72)
  footprint_hrly_array <- footprint[,,i,] 
  
  # now loop over all 72 hours for hour i 
  for (j in 1:length(footprint_hrly_array[1,1,])){
    
    # make an empty 99x99 raster in WGS1984 coordinates 
    empty_raster_wgs <- raster(nrows=99,ncols=99,crs='+proj=longlat +datum=WGS84')
    
    # make empty raster in WGS1984 coordinates with the wrf grid extents for temperature and radiation
    footprint_geo_hr_ind <- setExtent(empty_raster_wgs,extent(-87.9244,-84.3795,38.5049,41.1683))

    # put the (rotated) values from the temperature matrix in the empty raster (note:
    #  the apply  written rotates the grid into the correct orientation)
    values(footprint_geo_hr_ind) <- apply(t(footprint_hrly_array[,,j]),2,rev)
    
    # add the raster of layer j to the footprint stack (will be 72 layers at the end)
    footprint_hrly_stack <- stack(footprint_hrly_stack,footprint_geo_hr_ind)
    
    # clear the layer j raster for the next loop 
    rm(footprint_geo_hr_ind,empty_raster_wgs)
  }
  
  print(paste('saving hour',i,sep=' '))
  
  # convert to a raster brick before saving (idk why but writeRaster started returning
  # an error without this)
  footprint_hrly_stack <- brick(footprint_hrly_stack)
  
  # write the raster stack (99x99x72) for hour i of the month 
  writeRaster(footprint_hrly_stack,filename=paste0(footprint_raster_path_month,yr_month,'_hour_',i,'.tif'),overwrite=TRUE)
  
  # clear the stack for hour i for the next loop
  rm(footprint_hrly_stack)
  
}

# clear i and j 
rm(i)
rm(j)


print('------------------ done :) ------------------')


