##### script written by: Sam Murphy
##### purpose: extract PFT-specific VPRM NEE results at flux tower sites 

# load in libraries 
library(raster)
library(sp)
#library(rgdal)
library(lubridate)

# making a temporary directory - this is for intermediate files 
#tmpdir <- paste0('/gpfs/group/nlm136/INFLUX/SLM/temp_raster_vprm_extract') 
tmpdir <- paste0('/storage/group/zrb5027/default/INFLUX/smm8236/temporary_dirs/temp_raster_vprm_extract') 
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir

# set the year to run for this 
year = '2018'

# vegetation of the site to run
#pft_type <- 'Corn'
pft_type <- 'Other_Crops'

# set the site name (for saving file)
site_name <- 'A_NWa'
#site_name <- 'A_NWb'
#site_name <- 'A_09a'
#site_name <- 'A_09b'
#site_name <- 'A_14a'

# set the site coordinates 
#site_coords <- cbind(-86.5083,39.8632) # 2020 site 3 (SiteA_NWb)
site_coords <- cbind(-86.5005,39.8645) #  (SiteA_NWa)
#site_coords <- cbind(-85.7005,39.8167) # site A_09a 
#site_coords <- cbind(-85.7066,39.8766) # site A_09b
#site_coords <- cbind(-86.7880,39.9468) # site A_14a 

# set path to my outputs 
vprm_outputs_path <- paste0('/gpfs/group/nlm136/INFLUX/SLM/VPRM_gridded_run_outputs/',year,'*','/',pft_type,'/')

# set path to the output folder
path_out <- '/gpfs/group/nlm136/INFLUX/SLM/flux_comparisons_vprm/vprm_out_point/'

# get the filenames of nee files for the pft type
filenames_nee <- Sys.glob(paste(vprm_outputs_path,'*','NEE*',pft_type,'*',sep='')) 

# get the date and time for all of the files - need to change this for difference veg type names 
#dates.all <- strptime(substring(basename(filenames_nee),19,28),tz='UTC',format='%Y%m%d%H') # this is the number for corn
dates.all <- strptime(substring(basename(filenames_nee),26,35),tz='UTC',format='%Y%m%d%H') # this is the number for other crops

# make an empty vector to store the nee 
nee_vec_at_point <- rep(NA,length(filenames_nee))

# loop over all the rasters in the directory, extract the nee values at the cell where the site is located
for (i in 1:length(filenames_nee)){
  
  # read in a vprm nee output 
  nee_layer <- raster(filenames_nee[i])
  
  # pull the nee at the site location (grid cell the site falls in)
  nee_vec_at_point[i] <- extract(nee_layer,site_coords)
  
  # clear the nee layer 
  rm(nee_layer)
  print(i)
}

# now divide by the factor of 1000
nee_vec_at_point <- nee_vec_at_point/1000

# quick plot of the data 
plot(nee_vec_at_point,type='l',xlab='time',ylab='co2 flux (micromols/m2/s)')


# now let's save the model output nee
flux_matrix <- cbind(as.character(dates.all),nee_vec_at_point)

# write the output file 
write.csv(flux_matrix,file=paste0(path_out,'vprm_nee_Site_',site_name,'_',pft_type,'_',year,'.csv'))







