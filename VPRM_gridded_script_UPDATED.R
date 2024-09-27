# script adapted by sam murphy (originally written by Sharon Gourdji)
# purpose: run Sharon Gourdji's gridded updated VPRM for Indy, with individual PFT saving

# note for 2021** need to adjust the min/max lswi and evi raster file names 

# load in libraries 
library(raster)
library(sp)
library(rgdal) # i dont actually think you need this 
library(lubridate)

# allows for using the shell script
args <- commandArgs(TRUE)   
DATE<- args[1] #Year & month to run
domain<-args[2]

###################### set year and month to run ##############################
# YYYYMM to run 
DATE = '202112' # comment out if using shell script

# year and month separately (need this format for pulling radiatino and temp data)
YEAR <- '2021'  # YYYY
MONTH <- '12'   # MM

# set year of impervious surface area data (not a yearly data set)
impervious_year <- '2019'

# set domain (sharon was using multiple - for me d01 is the 3km wrf grid)
domain='d01'  

# print out dates for the log
print(DATE)
print(YEAR)
print(MONTH)
print(impervious_year)

############################## toggles #######################################
# Note: some toggles are leftover from old versions of script

URBANCORRECTION <- TRUE # if true, the impervious surface area is factored in
NEEOUT <- TRUE          # Output total NEE 
CATSOUT <- TRUE # Output 2 files: RES and GEE with one layer for category
PARAMSFROMFILE <- TRUE  # i think this just means that parameters are read from 
# the csv file instead of input manually (this is a leftover?)

PARAMSPERMONTH <- FALSE # not sure, i guess a leftover 

# unweighted option is left over from previous versions of this script (dont need)
UNWEIGHTED <- FALSE  #veg categories for GEE & RESP weighted by vegfrac, unless this set to T

SAVE_PFT_OUTS <- TRUE #This will save the resp,gpp,and nee of individual PFTs

####################### set paths to inputs and output ########################

# path to my directory on icds 
root.sam <- '/gpfs/group/nlm136/INFLUX/SLM/'  

#path to meteo data (radiations and temp)
path.temp <- paste0(root.sam,'VPRM_gridded_inputs/wrf_var_out/',YEAR,'/')   
path.srad <- paste0(root.sam,'VPRM_gridded_inputs/wrf_var_out/',YEAR,'/')

# path to processed EVI data 
path.EVI <- paste0(root.sam,'VPRM_gridded_inputs/modis_evi_processed_gapfill/') 

# path to processed LSWI data 
path.LSWI <- paste0(root.sam,'VPRM_gridded_inputs/modis_lswi_processed_gapfill/')  

# set path to the land use data (ie pft rasters)
path.landuse <- paste0(root.sam,'VPRM_gridded_inputs/cropland_processed/') 

# path to parameter file 
path.param <- paste0(root.sam,'VPRM_gridded_inputs/landcover_params/')  

# path to save outputs 
path.out <- paste0(root.sam,'VPRM_gridded_run_outputs/',DATE,'/')  
# create a directory to save the specified year and month outputs
dir.create(path.out)

# making a temporary directory 
tmpdir <- paste0(root.sam,'vprm_temp_',DATE,'/') # do this for rasters - otherwise 
# it will mess up some computing stuff 
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir



################# load evi and lswi input filenames and date info ################

# load filenames for LSWI data (sys.glob allows for dealing with wildcard)
filenamesLSWI <- Sys.glob(paste(path.LSWI,'LSWI_','*',sep='')) 

# pull the dates of LSWI data from the filenames 
datesLSWI <- strptime(substring(basename(filenamesLSWI),10,20),tz='UTC',format='%Y_doy%j')

# remove overlapping Julian dates at start of year (?)
ind.na = which(is.na(datesLSWI))  

# loop does something about overlapping dates at the start of year (?)
if (length(ind.na)>0) {
  filenamesLSWI=filenamesLSWI[-ind.na]
  datesLSWI=datesLSWI[-ind.na]}

# load filenames of EVI data 
filenamesEVI <- Sys.glob(paste(path.EVI,'EVI_','*','.tif',sep='')) # evi file paths names

# pull the dates of EVI data from the filenames 
datesEVI <- strptime(substring(basename(filenamesEVI),9,19),tz='UTC',format='%Y_doy%j')   # evi dates 

# remove overlapping Julian dates at start of year (?)
ind.na = which(is.na(datesEVI))  

# loop does something about overlapping dates at the start of year (?)
if (length(ind.na)>0) {
  filenamesEVI=filenamesEVI[-ind.na]
  datesEVI=datesEVI[-ind.na]}


################## load in min/max LSWI and EVI rasters ########################

evi_max <- raster(paste0(root.sam,'VPRM_gridded_inputs/modis_min_max/','max_evi_raster_2021.tif'))  
evi_min <- raster(paste0(root.sam,'VPRM_gridded_inputs/modis_min_max/','min_evi_raster_2021.tif')) 
lswi_max <- raster(paste0(root.sam,'VPRM_gridded_inputs/modis_min_max/','max_lswi_raster_2021.tif')) 
lswi_min <- raster(paste0(root.sam,'VPRM_gridded_inputs/modis_min_max/','min_lswi_raster_2021.tif')) 

########### calculate EVI threshold for phenology scaling factor ###############
evithresh = evi_min + 0.55*(evi_max-evi_min)     # solves for evi threshold

######################## load vegetation parameters #############################

# assign numbers to the PFT categories (according to sharon's category numbers)
cats <- c(83,82,23,90,71,47,45,52,46)   
# labels that go with the category numbers 
cats.label <- c('Corn','Other_Crops','Urban','Wetlands','Herbaceuous_Grassland_Hay_Pasture',
                'Deciduous_Forest','Evergreen_Forest_North','Shrub_Scrub','Evergreen_Forest_South')

#create PFT directories - for saving PFT NEE/GEE/RESP layers 
if(SAVE_PFT_OUTS){
  pftdir.list <- c()
  for (m in 1:length(cats)){
    print('Creating Individual PFT Directories ... ')
    pftdir <- paste0(path.out,cats.label[m],'/')
    pftdir.list[m] <- pftdir
    
    dir.create(pftdir)
  }  
}

# load in csv with parameters for each PFT 
params <- read.csv(paste0(path.param,'vprm_params_20210727','.csv'))

# make a vector for each parameter 
Tmin <- as.numeric(params[1,2:dim(params)[2]])
Topt <- as.numeric(params[2,2:dim(params)[2]])
Tmax <- as.numeric(params[3,2:dim(params)[2]])

Tcrit <- as.numeric(params[4,2:dim(params)[2]])
Tmod <- as.numeric(params[5,2:dim(params)[2]])

lambda <- as.numeric(params[6,2:dim(params)[2]])
PAR0 <- as.numeric(params[7,2:dim(params)[2]])
beta <- as.numeric(params[8,2:dim(params)[2]])
alpha <- as.numeric(params[9,2:dim(params)[2]])
alpha2 <- as.numeric(params[10,2:dim(params)[2]])
gamma <- as.numeric(params[11,2:dim(params)[2]])
theta1 <- as.numeric(params[12,2:dim(params)[2]])
theta2 <- as.numeric(params[13,2:dim(params)[2]])
theta3 <- as.numeric(params[14,2:dim(params)[2]])


######## load in pft fractions for sum and urban imperviousness data ############

# load in the fraction of vegetation cover in each gridcell and compute sum over all layers 
#vegfrac.sum <- stack(paste0(path.landuse,'landcover_percent_',YEAR,'.tif'))  
# sum all vegetation layers (excluding the missing data layer)
#vegfrac.sum <- sum(vegfrac.sum[[1:10]])

#load impervious surface area data 
if (URBANCORRECTION) {
  isa = raster(paste0(root.sam,'VPRM_gridded_inputs/urban_correction_data/','indy_impervious_',impervious_year,'.tif'))
  isa = isa/100  #convert from % to fraction
}


########## load filenames and dates of meteo data #####################
filenames.temp <- Sys.glob(paste(path.temp,'T2_geo__',YEAR,'-',MONTH,'*00:00','*.tif',sep=''))
dates.temp <- strptime(substring(basename(filenames.temp),9,27),tz='UTC',format='%Y-%m-%d_%H:%M:%S')

filenames.srad <- Sys.glob(paste(path.srad,'SRAD_geo__*',YEAR,'-',MONTH,'*00:00','*.tif',sep=''))
dates.srad <- strptime(substring(basename(filenames.srad),11,29),tz='UTC',format='%Y-%m-%d_%H:%M:%S')

#create full date vector (just in case any gaps)
dates.all = seq.POSIXt(dates.temp[1],dates.temp[length(dates.temp)],by=60*60)

#double-check complete hourly coverage for met data
length(dates.temp)
length(dates.srad)
sum(dates.temp==dates.srad)  #dates should match up

#run a specific month
ind.mo = which(format(dates.all,'%Y%m')==DATE)
dates.all = dates.all[ind.mo]


######################## Loop runs VPRM ####################################


#loop through hours, load input data, run VPRM
for (f in 1:length(dates.all)){  #for all loaded dates 
  
  
  tic=Sys.time()  # returns current date 
  date.f <- dates.all[f] # save the date at step f in loop in date.f
  
  # print the date of data in the loop 
  print(paste('Working on date and time: ',date.f,sep=''))
  
  #load temp filenames for the date in the loop 
  # grep function looks for files matching the pattern (pulls the file with the 
  # correct date)
  f.t = grep(format(dates.all[f],'%Y-%m-%d_%H:%M:%S'),filenames.temp)  
  file.temp <- filenames.temp[f.t]
  
  #load srad filenames for the date in the loop 
  f.s = grep(format(dates.all[f],'%Y-%m-%d_%H:%M:%S'),filenames.srad)
  file.srad <- filenames.srad[f.s]
  
  # load gridded temp data 
  Tair <- try(raster(file.temp))
  
  # convert temperature from kelvin to celsius  
  Tair = Tair - 273.15 
  
  # load solar radiation data 
  RAD <- try(raster(file.srad))
  #reset negative values to zero before bias correction
  RAD[RAD<0] = 0  
  
  # convert from short wave downward solar radiation  to PAR
  # divide shortwave by 0.505 to get PAR (Mahadevan), multiply by 4.56 to convert W/m^2 to micromol/m^2/s (Wang)
  PAR = (RAD/0.505)#*4.56 
  # reset negative values to zero (shouldn't be any)
  PAR[PAR<0] = 0         
  
  # sum of all PAR is an indicator for day/night
  sumPAR <- sum(getValues(PAR),na.rm=T) 
  # small values set to zero
  sumPAR[sumPAR<1e3]<-0  
  
  #at start of day, load MODIS data & pre-calculate scaling factors for day
  if (format(dates.all[f],'%H')=='00')  { 
    
    #clean up first (just in case)
    rm(lswi,evi,Wscale1,Wscale2,Pscale1)
    
    #load daily MODIS files 
    lswi <- raster(filenamesLSWI[which(format(datesLSWI,'%Y-%m-%d')==format(date.f,'%Y-%m-%d'))])
    evi <- raster(filenamesEVI[which(format(datesEVI,'%Y-%m-%d')==format(date.f,'%Y-%m-%d'))])
    evi[evi<0] = 0  #reset negative values to 0
    evi[evi>1] = 1  # ensure no values > 1
    
    # calculate water stress scaling factors (one is used in gpp other used in resp)
    Wscale1 <- (1.+lswi)/(1.+lswi_max)
    Wscale2 <- (lswi-lswi_min)/(lswi_max-lswi_min)
    
    
    
    # bound Wscale scaling factors between 0 and 1
    Wscale1[Wscale1<0] = 0
    Wscale1[Wscale1>1] = 1
    Wscale2[Wscale2<0] = 0
    Wscale2[Wscale2>1] = 1
    
    # effect of leaf phenology
    Pscale1 = (1.+lswi)/2.
    
    # get the y and x resolution of lswi (should be 0.01)
    xres <- xres(lswi)
    yres <- yres(lswi)
    
    # get the min and max in the x and y direction 
    xmin <- xmin(lswi)
    xmax <- xmax(lswi)
    ymin <- ymin(lswi)
    ymax <- ymax(lswi)
  } 
  
  
  
  print(paste('Computing common layers at',date()))
  
  # if no solar radiation it is nighttime
  if (sumPAR==0){print('Nighttime so only RESP will be computed')}
  
  #### load the fraction of vegetation type in each cell #####
  # load the stack of rasters (10 layers total)
  load_vegfrac_all <- stack(paste(path.landuse,'landcover_percent_',YEAR,'.tif',sep='')) 
  # note: when the missing data category (layer 1) is included, then the sum of all the layers is 1 
  
  # make a new stack without the first layer (the first layer is the missing data layer) 
  load_vegfrac <- subset(load_vegfrac_all,2:10)
  
  ###### loop through the different pfts ######
  for (m in 1:length(cats)){
    
    print(paste('Working on category: ',cats.label[m],' at ',date(),sep=''))
    
    # pull the fraction of the pft we are at in the loop
    vegfrac <- load_vegfrac[[m]]
    # add in lines to do the saving of the PFT layers 
    PFT_filter <- vegfrac
    # if there is no frac of pft in layer, set to zero
    PFT_filter[PFT_filter<=0] <- 0
    # if there is PFT, set to 1 
    PFT_filter[PFT_filter>0] <- 1
    
    #for daytime, calculate GEE
    if (sumPAR > 0){  # if the PAR is > 0 (ie daylight)
      
      a1 = Tair - Tmin[m]  #degrees above Tmin
      a2 = Tair - Tmax[m]  #degrees above Tmax
      a3 = Tair - Topt[m]  #degrees above Topt
      Tscale = a1*a2/(a1*a2 - a3**2)
      
      # Here a1 or a2 can't be negative
      Tscale = Tscale*(a1>0 & a2<0)
      PARscale = 1./(1. + PAR/PAR0[m])
      
      
      Pscale <- Pscale1
      if (cats[m]==45|cats[m]==46) {  # evergreen needleleaf
        Pscale[] <- 1.
      } else if (cats[m]==71 ) {      # grassland (or savanna but not in NLCD )
      } else {                        # Other vegetation types
        Pscale = (evi >= evithresh) + Pscale*(evi < evithresh)      
      }
      
      #make sure GEE components look good (scaling factors between 0 and 1)
      test.min = min(getValues(Tscale),na.rm=T)>=0 & min(getValues(Pscale),na.rm=T)>=0 & min(getValues(Wscale1),na.rm=T)>=0 & min(getValues(Wscale2),na.rm=T)>=0  & min(getValues(PARscale),na.rm=T)>=0
      test.max = max(getValues(Tscale),na.rm=T)<=1 & max(getValues(Pscale),na.rm=T)<=1 & max(getValues(Wscale1),na.rm=T)<=1 & max(getValues(Wscale2),na.rm=T)<=1 & max(getValues(PARscale),na.rm=T)<=1
      if (test.min==F | test.max==F) {
        ranges = rbind(range(getValues(Tscale),na.rm=T),range(getValues(Pscale),na.rm=T),range(getValues(Wscale1),na.rm=T),range(getValues(Wscale2),na.rm=T),range(getValues(PARscale),na.rm=T))
        rownames(ranges) = c('Tscale','Pscale','Wscale1','Wscale2','PARscale')
        write.csv(ranges,paste0(path.out,'Error_scalingFactors_',date.f,'.csv'))
      }
      
      #calculate GEE
      print(paste('Computing GEE at',date()))
      
      GEE_frac = lambda[m]*Tscale*Pscale*Wscale1*PARscale*evi*PAR*PFT_filter
      
      #write out category GEE values - makes a stack of gee for each PFT
      #save to stack instead (in memory)
      if (m==1) {gee_no_mod = GEE_frac}  else{gee_no_mod = stack(gee_no_mod,GEE_frac)}
      
      GEE_frac = GEE_frac * vegfrac
      
      #write out category GEE values - makes a stack of gee for each PFT
      #save to stack instead (in memory)
      if (m==1) {gee = GEE_frac}  else{gee = stack(gee,GEE_frac)}
      
    } 
    
    #calculate respiration (both day and night)
    print(paste('Computing RESP at',date()))
    ptcrit <- Tcrit[m]
    pTmod <-Tmod[m]
    palpha <- alpha[m]
    palpha2 <- alpha2[m]
    Pbeta <- beta[m]
    Pgamma <- gamma[m]
    Ptheta1 <- theta1[m]
    Ptheta2 <- theta2[m]
    Ptheta3 <- theta3[m]
    
    #Create new variable modifying temperature at low values by PFT for respiration calculation
    Tair.2 = Tair
    
    # note: this equation just has a -1 factored out of the parenthesis compared to the published paper 
    Tair.2[Tair.2<ptcrit] = (Tair.2[Tair.2<ptcrit] - ptcrit)*pTmod + ptcrit
    
    
    RESP_frac =(palpha*Tair.2+Pbeta+Pgamma*evi+palpha2*Tair.2^2+ Ptheta1*Wscale2 + Ptheta2*Wscale2*Tair.2 + Ptheta3*Wscale2*Tair.2^2) * PFT_filter 
    
    
    if (cats[m]==23 & URBANCORRECTION) {   # urban: modify heterotrophic respiration for now, re Hardiman et al
      RESP_frac = RESP_frac * (1-0.5*isa)  # same as resp/2 + resp/2*(1-ISA)
    } 
    
    # make sure there is only positive respiration values 
    RESP_frac = RESP_frac*(RESP_frac>0)  
    
    #write to file or stack (in memory) which ARE NOT modified by veg fraction
    if (m==1) {resp_no_mod = RESP_frac}  else{resp_no_mod = stack(resp_no_mod,RESP_frac)}
    
    RESP_frac = RESP_frac * vegfrac
    
    #write to file or stack (in memory)
    if (m==1) {resp = RESP_frac}  else{resp = stack(resp,RESP_frac)}
    
    removeTmpFiles(h=0)
    rm(Tair.2)
  } # end m categories
  
  removeTmpFiles(h=0)
  
  #Calculate total GEE across categories 
  if (sumPAR>0){
    #calculate weighted GEE across land-use categories
    Tgee = sum(gee)
    # Tgee[vegfrac.sum>1] = Tgee[vegfrac.sum>1]/vegfrac.sum[vegfrac.sum>1]  #divide out by sum of vegfrac categories if greater than one (if <1 because of water)
    
    
    #Tgee = Tgee /vegfrac.sum  #divide out by sum of vegfrac categories to get weighted average on land
    #gee <- stack(gee,Tgee)  #ONLY NEED THIS FOR SAVING OUT INDIVIDUAL PFT'S
  }
  
  #calculate weighted resp across land-use categories
  Tresp = sum(resp)
  
  #Tresp = Tresp/vegfrac.sum
  #resp <- stack(resp,Tresp)  #ONLY NEED THIS FOR SAVING OUT INDIVIDUAL PFT'S
  
  if(SAVE_PFT_OUTS){
    print('Saving Individual PFT Rasters ... ')
    for (m in 1:length(cats)){
      
      vegfrac <- load_vegfrac[[m]]
      
      PFT_filter <- vegfrac
      
      PFT_filter[PFT_filter<=0] <- 0
      
      PFT_filter[PFT_filter>0] <- 1
      
      if (sumPAR>0){
        nee_no_mod = (gee_no_mod[[m]] + resp_no_mod[[m]])*PFT_filter
      } else {
        nee_no_mod = (resp_no_mod[[m]])*PFT_filter
      }
      writeRaster(resp_no_mod[[m]]*1000,filename = paste(pftdir.list[m],'VPRM_RESP_',cats.label[m],'_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                  datatype='INT4S',options=c("COMPRESS=LZW"))
      writeRaster(nee_no_mod*1000,filename = paste(pftdir.list[m],'VPRM_NEE_',cats.label[m],'_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                  datatype='INT4S',options=c("COMPRESS=LZW"))
      if (sumPAR>0){
        writeRaster(gee_no_mod[[m]]*1000,filename = paste(pftdir.list[m],'VPRM_GEE_',cats.label[m],'_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                    datatype='INT4S',options=c("COMPRESS=LZW"))
      } 
    }  
  }
  
  #calculate NEE (GEE + resp)
  if (NEEOUT){
    print(paste('Computing NEE at',date()))
    if (sumPAR>0){
      neeTotal = Tgee + Tresp
    } else {
      neeTotal = Tresp
    }
  }
  
  #remove vegfrac multiplication if unweighted = T, but unweighted is not ever set to TRUE
  if (UNWEIGHTED){
    
    for (m in 1:length(cats)){
      
      vegfrac <- raster(paste(path.landuse,'vprm_cat_frac_',cats[m],'.tif',sep=''))
      
      if (sumPAR>0){
        gee[[m]] = gee[[m]]#/vegfrac
      }
      resp[[m]] = resp[[m]]#/vegfrac
    }
  }
  
  print(paste('Saving output at',date())) 
  
  
  #save out veg categories, GEE + RESP (if CATSOUT=T), only totals for now
  #multiply fluxes by 1000 and save as INT4S (integer) to save space; have to divide later
  if (CATSOUT){
    if (sumPAR>0){
      writeRaster(Tgee*1000,filename = paste(path.out,'VPRM_GEE_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                  datatype='INT4S',options=c("COMPRESS=LZW"))
    }
    writeRaster(Tresp*1000,filename = paste(path.out,'VPRM_RESP_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                datatype='INT4S',options=c("COMPRESS=LZW"))
  }
  #save NEE, if NEEOUT=T
  if (NEEOUT){
    writeRaster(neeTotal*1000,filename = paste(path.out,'VPRM_NEE_',domain,'.',format(date.f,tz='UTC',format='%Y%m%d%H'),'.tif',sep=''),overwrite=T,
                datatype='INT4S',options=c("COMPRESS=LZW"))
  }
  
  
  
  #clean up workspace
  rm(gee,gee_no_mod,resp,resp_no_mod,nee_no_mod,neeTotal,Tgee,Tresp,Tair,Tair.2,RAD,PAR,Tscale,Wscale,Pscale,PARscale,GEE_frac,RESP_frac,vegfrac,PFT_filter)  
  
  removeTmpFiles(h=0)
  print(paste('Finished Date: ',date.f,' at ',date(),sep=''))
  toc=Sys.time()
} # end f dates

print('-------------------------FINISHED VPRM GRIDDED RUN--------------------------------')
