#### script written by: Sam Murphy
#### purpose: look at the landcover in model domain - percent of domain for each PFT

# load in libraries 
library(raster)
library(sp)
#library(rgdal)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(colorspace)
#library(SpaDES.tools)

# making a temporary directory - this is for intermediate files 
tmpdir <- paste0('/storage/group/zrb5027/default/INFLUX/smm8236/temporary_dirs/temp_raster_vprm_in_out_explore') 
dir.create(tmpdir)  # creates a directory names tmpdir at the specified path 
rasterOptions(tmpdir=tmpdir)  # tells r that the temp directory is called tmpdir


##################### read in the vegfrac layers ########################
path_vegfrac <- '/storage/group/zrb5027/default/INFLUX/smm8236/VPRM_gridded_inputs/cropland_processed/'

vegfrac_2021 <- stack(paste0(path_vegfrac,'landcover_percent_2021.tif'))
vegfrac_2020 <- stack(paste0(path_vegfrac,'landcover_percent_2020.tif'))
vegfrac_2019 <- stack(paste0(path_vegfrac,'landcover_percent_2019.tif'))
vegfrac_2018 <- stack(paste0(path_vegfrac,'landcover_percent_2018.tif'))
vegfrac_2017 <- stack(paste0(path_vegfrac,'landcover_percent_2017.tif'))
vegfrac_2016 <- stack(paste0(path_vegfrac,'landcover_percent_2016.tif'))
vegfrac_2015 <- stack(paste0(path_vegfrac,'landcover_percent_2015.tif'))
vegfrac_2014 <- stack(paste0(path_vegfrac,'landcover_percent_2014.tif'))
vegfrac_2013 <- stack(paste0(path_vegfrac,'landcover_percent_2013.tif'))
vegfrac_2012 <- stack(paste0(path_vegfrac,'landcover_percent_2012.tif'))

# layer order - first layer in the stack is the no data / water category 
# ('Corn','Other_Crops','Urban','Wetlands','Herbaceuous_Grassland_Hay_Pasture',
# 'Deciduous_Forest','Evergreen_Forest_North','Shrub_Scrub','Evergreen_Forest_South')

# make data frames for each of the layers for each year for plotting with ggplot
no_data_2021 <- as.data.frame(vegfrac_2021[[1]],xy=TRUE)
no_data_2020 <- as.data.frame(vegfrac_2020[[1]],xy=TRUE)
no_data_2019 <- as.data.frame(vegfrac_2019[[1]],xy=TRUE)
no_data_2018 <- as.data.frame(vegfrac_2018[[1]],xy=TRUE)
no_data_2017 <- as.data.frame(vegfrac_2017[[1]],xy=TRUE)
no_data_2016 <- as.data.frame(vegfrac_2016[[1]],xy=TRUE)
no_data_2015 <- as.data.frame(vegfrac_2015[[1]],xy=TRUE)
no_data_2014 <- as.data.frame(vegfrac_2014[[1]],xy=TRUE)
no_data_2013 <- as.data.frame(vegfrac_2013[[1]],xy=TRUE)
no_data_2012 <- as.data.frame(vegfrac_2012[[1]],xy=TRUE)

corn_2021 <- as.data.frame(vegfrac_2021[[2]],xy=TRUE)
corn_2020 <- as.data.frame(vegfrac_2020[[2]],xy=TRUE)
corn_2019 <- as.data.frame(vegfrac_2019[[2]],xy=TRUE)
corn_2018 <- as.data.frame(vegfrac_2018[[2]],xy=TRUE)
corn_2017 <- as.data.frame(vegfrac_2017[[2]],xy=TRUE)
corn_2016 <- as.data.frame(vegfrac_2016[[2]],xy=TRUE)
corn_2015 <- as.data.frame(vegfrac_2015[[2]],xy=TRUE)
corn_2014 <- as.data.frame(vegfrac_2014[[2]],xy=TRUE)
corn_2013 <- as.data.frame(vegfrac_2013[[2]],xy=TRUE)
corn_2012 <- as.data.frame(vegfrac_2012[[2]],xy=TRUE)

other_crops_2021 <- as.data.frame(vegfrac_2021[[3]],xy=TRUE)
other_crops_2020 <- as.data.frame(vegfrac_2020[[3]],xy=TRUE)
other_crops_2019 <- as.data.frame(vegfrac_2019[[3]],xy=TRUE)
other_crops_2018 <- as.data.frame(vegfrac_2018[[3]],xy=TRUE)
other_crops_2017 <- as.data.frame(vegfrac_2017[[3]],xy=TRUE)
other_crops_2016 <- as.data.frame(vegfrac_2016[[3]],xy=TRUE)
other_crops_2015 <- as.data.frame(vegfrac_2015[[3]],xy=TRUE)
other_crops_2014 <- as.data.frame(vegfrac_2014[[3]],xy=TRUE)
other_crops_2013 <- as.data.frame(vegfrac_2013[[3]],xy=TRUE)
other_crops_2012 <- as.data.frame(vegfrac_2012[[3]],xy=TRUE)

urban_2021 <- as.data.frame(vegfrac_2021[[4]],xy=TRUE)
urban_2020 <- as.data.frame(vegfrac_2020[[4]],xy=TRUE)
urban_2019 <- as.data.frame(vegfrac_2019[[4]],xy=TRUE)
urban_2018 <- as.data.frame(vegfrac_2018[[4]],xy=TRUE)
urban_2017 <- as.data.frame(vegfrac_2017[[4]],xy=TRUE)
urban_2016 <- as.data.frame(vegfrac_2016[[4]],xy=TRUE)
urban_2015 <- as.data.frame(vegfrac_2015[[4]],xy=TRUE)
urban_2014 <- as.data.frame(vegfrac_2014[[4]],xy=TRUE)
urban_2013 <- as.data.frame(vegfrac_2013[[4]],xy=TRUE)
urban_2012 <- as.data.frame(vegfrac_2012[[4]],xy=TRUE)

Wetlands_2021 <- as.data.frame(vegfrac_2021[[5]],xy=TRUE)
Wetlands_2020 <- as.data.frame(vegfrac_2020[[5]],xy=TRUE)
Wetlands_2019 <- as.data.frame(vegfrac_2019[[5]],xy=TRUE)
Wetlands_2018 <- as.data.frame(vegfrac_2018[[5]],xy=TRUE)
Wetlands_2017 <- as.data.frame(vegfrac_2017[[5]],xy=TRUE)
Wetlands_2016 <- as.data.frame(vegfrac_2016[[5]],xy=TRUE)
Wetlands_2015 <- as.data.frame(vegfrac_2015[[5]],xy=TRUE)
Wetlands_2014 <- as.data.frame(vegfrac_2014[[5]],xy=TRUE)
Wetlands_2013 <- as.data.frame(vegfrac_2013[[5]],xy=TRUE)
Wetlands_2012 <- as.data.frame(vegfrac_2012[[5]],xy=TRUE)

Herbaceuous_Grassland_Hay_Pasture_2021 <- as.data.frame(vegfrac_2021[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2020 <- as.data.frame(vegfrac_2020[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2019 <- as.data.frame(vegfrac_2019[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2018 <- as.data.frame(vegfrac_2018[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2017 <- as.data.frame(vegfrac_2017[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2016 <- as.data.frame(vegfrac_2016[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2015 <- as.data.frame(vegfrac_2015[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2014 <- as.data.frame(vegfrac_2014[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2013 <- as.data.frame(vegfrac_2013[[6]],xy=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2012 <- as.data.frame(vegfrac_2012[[6]],xy=TRUE)

Deciduous_Forest_2021 <- as.data.frame(vegfrac_2021[[7]],xy=TRUE)
Deciduous_Forest_2020 <- as.data.frame(vegfrac_2020[[7]],xy=TRUE)
Deciduous_Forest_2019 <- as.data.frame(vegfrac_2019[[7]],xy=TRUE)
Deciduous_Forest_2018 <- as.data.frame(vegfrac_2018[[7]],xy=TRUE)
Deciduous_Forest_2017 <- as.data.frame(vegfrac_2017[[7]],xy=TRUE)
Deciduous_Forest_2016 <- as.data.frame(vegfrac_2016[[7]],xy=TRUE)
Deciduous_Forest_2015 <- as.data.frame(vegfrac_2015[[7]],xy=TRUE)
Deciduous_Forest_2014 <- as.data.frame(vegfrac_2014[[7]],xy=TRUE)
Deciduous_Forest_2013 <- as.data.frame(vegfrac_2013[[7]],xy=TRUE)
Deciduous_Forest_2012 <- as.data.frame(vegfrac_2012[[7]],xy=TRUE)

Evergreen_Forest_North_2021 <- as.data.frame(vegfrac_2021[[8]],xy=TRUE)
Evergreen_Forest_North_2020 <- as.data.frame(vegfrac_2020[[8]],xy=TRUE)
Evergreen_Forest_North_2019 <- as.data.frame(vegfrac_2019[[8]],xy=TRUE)
Evergreen_Forest_North_2018 <- as.data.frame(vegfrac_2018[[8]],xy=TRUE)
Evergreen_Forest_North_2017 <- as.data.frame(vegfrac_2017[[8]],xy=TRUE)
Evergreen_Forest_North_2016 <- as.data.frame(vegfrac_2016[[8]],xy=TRUE)
Evergreen_Forest_North_2015 <- as.data.frame(vegfrac_2015[[8]],xy=TRUE)
Evergreen_Forest_North_2014 <- as.data.frame(vegfrac_2014[[8]],xy=TRUE)
Evergreen_Forest_North_2013 <- as.data.frame(vegfrac_2013[[8]],xy=TRUE)
Evergreen_Forest_North_2012 <- as.data.frame(vegfrac_2012[[8]],xy=TRUE)

Shrub_Scrub_2021 <- as.data.frame(vegfrac_2021[[9]],xy=TRUE)
Shrub_Scrub_2020 <- as.data.frame(vegfrac_2020[[9]],xy=TRUE)
Shrub_Scrub_2019 <- as.data.frame(vegfrac_2019[[9]],xy=TRUE)
Shrub_Scrub_2018 <- as.data.frame(vegfrac_2018[[9]],xy=TRUE)
Shrub_Scrub_2017 <- as.data.frame(vegfrac_2017[[9]],xy=TRUE)
Shrub_Scrub_2016 <- as.data.frame(vegfrac_2016[[9]],xy=TRUE)
Shrub_Scrub_2015 <- as.data.frame(vegfrac_2015[[9]],xy=TRUE)
Shrub_Scrub_2014 <- as.data.frame(vegfrac_2014[[9]],xy=TRUE)
Shrub_Scrub_2013 <- as.data.frame(vegfrac_2013[[9]],xy=TRUE)
Shrub_Scrub_2012 <- as.data.frame(vegfrac_2012[[9]],xy=TRUE)

Evergreen_Forest_South_2021 <- as.data.frame(vegfrac_2021[[10]],xy=TRUE)
Evergreen_Forest_South_2020 <- as.data.frame(vegfrac_2020[[10]],xy=TRUE)
Evergreen_Forest_South_2019 <- as.data.frame(vegfrac_2019[[10]],xy=TRUE)
Evergreen_Forest_South_2018 <- as.data.frame(vegfrac_2018[[10]],xy=TRUE)
Evergreen_Forest_South_2017 <- as.data.frame(vegfrac_2017[[10]],xy=TRUE)
Evergreen_Forest_South_2016 <- as.data.frame(vegfrac_2016[[10]],xy=TRUE)
Evergreen_Forest_South_2015 <- as.data.frame(vegfrac_2015[[10]],xy=TRUE)
Evergreen_Forest_South_2014 <- as.data.frame(vegfrac_2014[[10]],xy=TRUE)
Evergreen_Forest_South_2013 <- as.data.frame(vegfrac_2013[[10]],xy=TRUE)
Evergreen_Forest_South_2012 <- as.data.frame(vegfrac_2012[[10]],xy=TRUE)

# plot the 2020 layers seperatly 
ggplot() + 
  geom_raster(data=no_data_2020, aes(x=x,y=y,fill=landcover_percent_2020_1)) +
  scale_fill_continuous_sequential(palette='Red-Yellow', name='No data/water') + 
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=corn_2020, aes(x=x,y=y,fill=landcover_percent_2020_2)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='corn') + 
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=other_crops_2020, aes(x=x,y=y,fill=landcover_percent_2020_3)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='crops') +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=urban_2020, aes(x=x,y=y,fill=landcover_percent_2020_4)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='urban') +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Wetlands_2020, aes(x=x,y=y,fill=landcover_percent_2020_5)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='wetland',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Herbaceuous_Grassland_Hay_Pasture_2020, aes(x=x,y=y,fill=landcover_percent_2020_6)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='Grasses',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Deciduous_Forest_2020, aes(x=x,y=y,fill=landcover_percent_2020_7)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='DF',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Evergreen_Forest_North_2020, aes(x=x,y=y,fill=landcover_percent_2020_8)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='Ever. N',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Shrub_Scrub_2020, aes(x=x,y=y,fill=landcover_percent_2020_9)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='Shrub',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')

ggplot() + 
  geom_raster(data=Evergreen_Forest_South_2020, aes(x=x,y=y,fill=landcover_percent_2020_10)) +
  scale_fill_continuous_sequential(palette='Red-Yellow',name='Ever. S',limits=c(0,1)) +
  xlab('Lon') + ylab('Lat')


############ total domain fractions ##########
###### 2021 ######
# calculate the sum over all the layers 
vegfrac_layer_sum_2021 <- sum(vegfrac_2021,na.rm=TRUE)
total_sum_2021 <- cellStats(vegfrac_layer_sum_2021,stat='sum')

no_data_2021_sum <- sum(no_data_2021$landcover_percent_2021_1,na.rm=TRUE)  
corn_2021_sum <- sum(corn_2021$landcover_percent_2021_2,na.rm=TRUE)
other_crops_2021_sum <- sum(other_crops_2021$landcover_percent_2021_3,na.rm=TRUE)
urban_2021_sum <- sum(urban_2021$landcover_percent_2021_4,na.rm=TRUE)
Wetlands_2021_sum <- sum(Wetlands_2021$landcover_percent_2021_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2021_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2021$landcover_percent_2021_6,na.rm=TRUE)
Deciduous_Forest_2021_sum <- sum(Deciduous_Forest_2021$landcover_percent_2021_7,na.rm=TRUE)
Evergreen_Forest_North_2021_sum <- sum(Evergreen_Forest_North_2021$landcover_percent_2021_8,na.rm=TRUE)
Shrub_Scrub_2021_sum <- sum(Shrub_Scrub_2021$landcover_percent_2021_9,na.rm=TRUE) 
Evergreen_Forest_South_2021_sum <- sum(Evergreen_Forest_South_2021$landcover_percent_2021_10,na.rm=TRUE) 

no_data_2021_domain_frac <- no_data_2021_sum/total_sum_2021
corn_2021_domain_frac <- corn_2021_sum/total_sum_2021
other_crops_2021_domain_frac <- other_crops_2021_sum/total_sum_2021
urban_2021_domain_frac <- urban_2021_sum/total_sum_2021
Wetlands_2021_domain_frac <- Wetlands_2021_sum/total_sum_2021
Herbaceuous_Grassland_Hay_Pasture_2021_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2021_sum/total_sum_2021
Deciduous_Forest_2021_domain_frac <- Deciduous_Forest_2021_sum/total_sum_2021
Evergreen_Forest_North_2021_domain_frac <- Evergreen_Forest_North_2021_sum/total_sum_2021
Shrub_Scrub_2021_domain_frac <- Shrub_Scrub_2021_sum/total_sum_2021
Evergreen_Forest_South_2021_domain_frac <- Evergreen_Forest_South_2021_sum/total_sum_2021

# look at the fractions 
no_data_2021_domain_frac 
corn_2021_domain_frac 
other_crops_2021_domain_frac 
urban_2021_domain_frac 
Wetlands_2021_domain_frac 
Herbaceuous_Grassland_Hay_Pasture_2021_domain_frac 
Deciduous_Forest_2021_domain_frac 
Evergreen_Forest_North_2021_domain_frac 
Shrub_Scrub_2021_domain_frac 
Evergreen_Forest_South_2021_domain_frac 

# check that the sum is one 
sum(no_data_2021_domain_frac, 
    corn_2021_domain_frac, 
    other_crops_2021_domain_frac, 
    urban_2021_domain_frac,
    Wetlands_2021_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2021_domain_frac, 
    Deciduous_Forest_2021_domain_frac,
    Evergreen_Forest_North_2021_domain_frac, 
    Shrub_Scrub_2021_domain_frac, 
    Evergreen_Forest_South_2021_domain_frac) 


###### 2020 ######
# calculate the sum over all the layers 
vegfrac_layer_sum_2020 <- sum(vegfrac_2020,na.rm=TRUE)
total_sum_2020 <- cellStats(vegfrac_layer_sum_2020,stat='sum')

no_data_2020_sum <- sum(no_data_2020$landcover_percent_2020_1,na.rm=TRUE)  
corn_2020_sum <- sum(corn_2020$landcover_percent_2020_2,na.rm=TRUE)
other_crops_2020_sum <- sum(other_crops_2020$landcover_percent_2020_3,na.rm=TRUE)
urban_2020_sum <- sum(urban_2020$landcover_percent_2020_4,na.rm=TRUE)
Wetlands_2020_sum <- sum(Wetlands_2020$landcover_percent_2020_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2020_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2020$landcover_percent_2020_6,na.rm=TRUE)
Deciduous_Forest_2020_sum <- sum(Deciduous_Forest_2020$landcover_percent_2020_7,na.rm=TRUE)
Evergreen_Forest_North_2020_sum <- sum(Evergreen_Forest_North_2020$landcover_percent_2020_8,na.rm=TRUE)
Shrub_Scrub_2020_sum <- sum(Shrub_Scrub_2020$landcover_percent_2020_9,na.rm=TRUE) 
Evergreen_Forest_South_2020_sum <- sum(Evergreen_Forest_South_2020$landcover_percent_2020_10,na.rm=TRUE) 

no_data_2020_domain_frac <- no_data_2020_sum/total_sum_2020
corn_2020_domain_frac <- corn_2020_sum/total_sum_2020
other_crops_2020_domain_frac <- other_crops_2020_sum/total_sum_2020
urban_2020_domain_frac <- urban_2020_sum/total_sum_2020
Wetlands_2020_domain_frac <- Wetlands_2020_sum/total_sum_2020
Herbaceuous_Grassland_Hay_Pasture_2020_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2020_sum/total_sum_2020
Deciduous_Forest_2020_domain_frac <- Deciduous_Forest_2020_sum/total_sum_2020
Evergreen_Forest_North_2020_domain_frac <- Evergreen_Forest_North_2020_sum/total_sum_2020
Shrub_Scrub_2020_domain_frac <- Shrub_Scrub_2020_sum/total_sum_2020
Evergreen_Forest_South_2020_domain_frac <- Evergreen_Forest_South_2020_sum/total_sum_2020

# look at the fractions 
no_data_2020_domain_frac 
corn_2020_domain_frac 
other_crops_2020_domain_frac 
urban_2020_domain_frac 
Wetlands_2020_domain_frac 
Herbaceuous_Grassland_Hay_Pasture_2020_domain_frac 
Deciduous_Forest_2020_domain_frac 
Evergreen_Forest_North_2020_domain_frac 
Shrub_Scrub_2020_domain_frac 
Evergreen_Forest_South_2020_domain_frac 

# check that the sum is one 
sum(no_data_2020_domain_frac, 
    corn_2020_domain_frac, 
    other_crops_2020_domain_frac, 
    urban_2020_domain_frac,
    Wetlands_2020_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2020_domain_frac, 
    Deciduous_Forest_2020_domain_frac,
    Evergreen_Forest_North_2020_domain_frac, 
    Shrub_Scrub_2020_domain_frac, 
    Evergreen_Forest_South_2020_domain_frac) 

slices_2020 <- c(no_data_2020_domain_frac, 
                 corn_2020_domain_frac, 
                 other_crops_2020_domain_frac, 
                 urban_2020_domain_frac,
                 Wetlands_2020_domain_frac,
                 Herbaceuous_Grassland_Hay_Pasture_2020_domain_frac, 
                 Deciduous_Forest_2020_domain_frac,
                 Evergreen_Forest_North_2020_domain_frac, 
                 Shrub_Scrub_2020_domain_frac, 
                 Evergreen_Forest_South_2020_domain_frac)

pie(slices_2020)

barplot(slices_2020)


###### look at 2019 
# calculate the sum over all the layers 
vegfrac_layer_sum_2019 <- sum(vegfrac_2019,na.rm=TRUE)
total_sum_2019 <- cellStats(vegfrac_layer_sum_2019,stat='sum')

no_data_2019_sum <- sum(no_data_2019$landcover_percent_2019_1,na.rm=TRUE)  
corn_2019_sum <- sum(corn_2019$landcover_percent_2019_2,na.rm=TRUE)
other_crops_2019_sum <- sum(other_crops_2019$landcover_percent_2019_3,na.rm=TRUE)
urban_2019_sum <- sum(urban_2019$landcover_percent_2019_4,na.rm=TRUE)
Wetlands_2019_sum <- sum(Wetlands_2019$landcover_percent_2019_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2019_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2019$landcover_percent_2019_6,na.rm=TRUE)
Deciduous_Forest_2019_sum <- sum(Deciduous_Forest_2019$landcover_percent_2019_7,na.rm=TRUE)
Evergreen_Forest_North_2019_sum <- sum(Evergreen_Forest_North_2019$landcover_percent_2019_8,na.rm=TRUE)
Shrub_Scrub_2019_sum <- sum(Shrub_Scrub_2019$landcover_percent_2019_9,na.rm=TRUE) 
Evergreen_Forest_South_2019_sum <- sum(Evergreen_Forest_South_2019$landcover_percent_2019_10,na.rm=TRUE) 

no_data_2019_domain_frac <- no_data_2019_sum/total_sum_2019
corn_2019_domain_frac <- corn_2019_sum/total_sum_2019
other_crops_2019_domain_frac <- other_crops_2019_sum/total_sum_2019
urban_2019_domain_frac <- urban_2019_sum/total_sum_2019
Wetlands_2019_domain_frac <- Wetlands_2019_sum/total_sum_2019
Herbaceuous_Grassland_Hay_Pasture_2019_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2019_sum/total_sum_2019
Deciduous_Forest_2019_domain_frac <- Deciduous_Forest_2019_sum/total_sum_2019
Evergreen_Forest_North_2019_domain_frac <- Evergreen_Forest_North_2019_sum/total_sum_2019
Shrub_Scrub_2019_domain_frac <- Shrub_Scrub_2019_sum/total_sum_2019
Evergreen_Forest_South_2019_domain_frac <- Evergreen_Forest_South_2019_sum/total_sum_2019

# look at the fractions 
no_data_2019_domain_frac 
corn_2019_domain_frac 
other_crops_2019_domain_frac 
urban_2019_domain_frac 
Wetlands_2019_domain_frac 
Herbaceuous_Grassland_Hay_Pasture_2019_domain_frac 
Deciduous_Forest_2019_domain_frac 
Evergreen_Forest_North_2019_domain_frac 
Shrub_Scrub_2019_domain_frac 
Evergreen_Forest_South_2019_domain_frac 

sum(no_data_2019_domain_frac, 
    corn_2019_domain_frac, 
    other_crops_2019_domain_frac, 
    urban_2019_domain_frac,
    Wetlands_2019_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2019_domain_frac, 
    Deciduous_Forest_2019_domain_frac,
    Evergreen_Forest_North_2019_domain_frac, 
    Shrub_Scrub_2019_domain_frac, 
    Evergreen_Forest_South_2019_domain_frac)


###### look at 2018 
# calculate the sum over all the layers 
vegfrac_layer_sum_2018 <- sum(vegfrac_2018,na.rm=TRUE)
total_sum_2018 <- cellStats(vegfrac_layer_sum_2018,stat='sum')

no_data_2018_sum <- sum(no_data_2018$landcover_percent_2018_1,na.rm=TRUE)  
corn_2018_sum <- sum(corn_2018$landcover_percent_2018_2,na.rm=TRUE)
other_crops_2018_sum <- sum(other_crops_2018$landcover_percent_2018_3,na.rm=TRUE)
urban_2018_sum <- sum(urban_2018$landcover_percent_2018_4,na.rm=TRUE)
Wetlands_2018_sum <- sum(Wetlands_2018$landcover_percent_2018_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2018_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2018$landcover_percent_2018_6,na.rm=TRUE)
Deciduous_Forest_2018_sum <- sum(Deciduous_Forest_2018$landcover_percent_2018_7,na.rm=TRUE)
Evergreen_Forest_North_2018_sum <- sum(Evergreen_Forest_North_2018$landcover_percent_2018_8,na.rm=TRUE)
Shrub_Scrub_2018_sum <- sum(Shrub_Scrub_2018$landcover_percent_2018_9,na.rm=TRUE) 
Evergreen_Forest_South_2018_sum <- sum(Evergreen_Forest_South_2018$landcover_percent_2018_10,na.rm=TRUE) 

no_data_2018_domain_frac <- no_data_2018_sum/total_sum_2018
corn_2018_domain_frac <- corn_2018_sum/total_sum_2018
other_crops_2018_domain_frac <- other_crops_2018_sum/total_sum_2018
urban_2018_domain_frac <- urban_2018_sum/total_sum_2018
Wetlands_2018_domain_frac <- Wetlands_2018_sum/total_sum_2018
Herbaceuous_Grassland_Hay_Pasture_2018_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2018_sum/total_sum_2018
Deciduous_Forest_2018_domain_frac <- Deciduous_Forest_2018_sum/total_sum_2018
Evergreen_Forest_North_2018_domain_frac <- Evergreen_Forest_North_2018_sum/total_sum_2018
Shrub_Scrub_2018_domain_frac <- Shrub_Scrub_2018_sum/total_sum_2018
Evergreen_Forest_South_2018_domain_frac <- Evergreen_Forest_South_2018_sum/total_sum_2018

# look at the fractions 
no_data_2018_domain_frac * 100
corn_2018_domain_frac * 100
other_crops_2018_domain_frac * 100
urban_2018_domain_frac * 100
Wetlands_2018_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2018_domain_frac * 100
Deciduous_Forest_2018_domain_frac * 100
Evergreen_Forest_North_2018_domain_frac * 100
Shrub_Scrub_2018_domain_frac * 100
Evergreen_Forest_South_2018_domain_frac * 100

# check that the sum is one 
sum(no_data_2018_domain_frac, 
    corn_2018_domain_frac, 
    other_crops_2018_domain_frac, 
    urban_2018_domain_frac,
    Wetlands_2018_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2018_domain_frac, 
    Deciduous_Forest_2018_domain_frac,
    Evergreen_Forest_North_2018_domain_frac, 
    Shrub_Scrub_2018_domain_frac, 
    Evergreen_Forest_South_2018_domain_frac)


###### look at 2017 
# calculate the sum over all the layers 
vegfrac_layer_sum_2017 <- sum(vegfrac_2017,na.rm=TRUE)
total_sum_2017 <- cellStats(vegfrac_layer_sum_2017,stat='sum')

no_data_2017_sum <- sum(no_data_2017$landcover_percent_2017_1,na.rm=TRUE)  
corn_2017_sum <- sum(corn_2017$landcover_percent_2017_2,na.rm=TRUE)
other_crops_2017_sum <- sum(other_crops_2017$landcover_percent_2017_3,na.rm=TRUE)
urban_2017_sum <- sum(urban_2017$landcover_percent_2017_4,na.rm=TRUE)
Wetlands_2017_sum <- sum(Wetlands_2017$landcover_percent_2017_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2017_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2017$landcover_percent_2017_6,na.rm=TRUE)
Deciduous_Forest_2017_sum <- sum(Deciduous_Forest_2017$landcover_percent_2017_7,na.rm=TRUE)
Evergreen_Forest_North_2017_sum <- sum(Evergreen_Forest_North_2017$landcover_percent_2017_8,na.rm=TRUE)
Shrub_Scrub_2017_sum <- sum(Shrub_Scrub_2017$landcover_percent_2017_9,na.rm=TRUE) 
Evergreen_Forest_South_2017_sum <- sum(Evergreen_Forest_South_2017$landcover_percent_2017_10,na.rm=TRUE) 

no_data_2017_domain_frac <- no_data_2017_sum/total_sum_2017
corn_2017_domain_frac <- corn_2017_sum/total_sum_2017
other_crops_2017_domain_frac <- other_crops_2017_sum/total_sum_2017
urban_2017_domain_frac <- urban_2017_sum/total_sum_2017
Wetlands_2017_domain_frac <- Wetlands_2017_sum/total_sum_2017
Herbaceuous_Grassland_Hay_Pasture_2017_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2017_sum/total_sum_2017
Deciduous_Forest_2017_domain_frac <- Deciduous_Forest_2017_sum/total_sum_2017
Evergreen_Forest_North_2017_domain_frac <- Evergreen_Forest_North_2017_sum/total_sum_2017
Shrub_Scrub_2017_domain_frac <- Shrub_Scrub_2017_sum/total_sum_2017
Evergreen_Forest_South_2017_domain_frac <- Evergreen_Forest_South_2017_sum/total_sum_2017

# look at the fractions 
no_data_2017_domain_frac * 100
corn_2017_domain_frac * 100
other_crops_2017_domain_frac * 100
urban_2017_domain_frac * 100
Wetlands_2017_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2017_domain_frac * 100
Deciduous_Forest_2017_domain_frac * 100
Evergreen_Forest_North_2017_domain_frac * 100
Shrub_Scrub_2017_domain_frac * 100
Evergreen_Forest_South_2017_domain_frac * 100

# check that the sum is one 
sum(no_data_2017_domain_frac, 
    corn_2017_domain_frac, 
    other_crops_2017_domain_frac, 
    urban_2017_domain_frac,
    Wetlands_2017_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2017_domain_frac, 
    Deciduous_Forest_2017_domain_frac,
    Evergreen_Forest_North_2017_domain_frac, 
    Shrub_Scrub_2017_domain_frac, 
    Evergreen_Forest_South_2017_domain_frac)

###### look at 2016 
# calculate the sum over all the layers 
vegfrac_layer_sum_2016 <- sum(vegfrac_2016,na.rm=TRUE)
total_sum_2016 <- cellStats(vegfrac_layer_sum_2016,stat='sum')

no_data_2016_sum <- sum(no_data_2016$landcover_percent_2016_1,na.rm=TRUE)  
corn_2016_sum <- sum(corn_2016$landcover_percent_2016_2,na.rm=TRUE)
other_crops_2016_sum <- sum(other_crops_2016$landcover_percent_2016_3,na.rm=TRUE)
urban_2016_sum <- sum(urban_2016$landcover_percent_2016_4,na.rm=TRUE)
Wetlands_2016_sum <- sum(Wetlands_2016$landcover_percent_2016_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2016_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2016$landcover_percent_2016_6,na.rm=TRUE)
Deciduous_Forest_2016_sum <- sum(Deciduous_Forest_2016$landcover_percent_2016_7,na.rm=TRUE)
Evergreen_Forest_North_2016_sum <- sum(Evergreen_Forest_North_2016$landcover_percent_2016_8,na.rm=TRUE)
Shrub_Scrub_2016_sum <- sum(Shrub_Scrub_2016$landcover_percent_2016_9,na.rm=TRUE) 
Evergreen_Forest_South_2016_sum <- sum(Evergreen_Forest_South_2016$landcover_percent_2016_10,na.rm=TRUE) 

no_data_2016_domain_frac <- no_data_2016_sum/total_sum_2016
corn_2016_domain_frac <- corn_2016_sum/total_sum_2016
other_crops_2016_domain_frac <- other_crops_2016_sum/total_sum_2016
urban_2016_domain_frac <- urban_2016_sum/total_sum_2016
Wetlands_2016_domain_frac <- Wetlands_2016_sum/total_sum_2016
Herbaceuous_Grassland_Hay_Pasture_2016_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2016_sum/total_sum_2016
Deciduous_Forest_2016_domain_frac <- Deciduous_Forest_2016_sum/total_sum_2016
Evergreen_Forest_North_2016_domain_frac <- Evergreen_Forest_North_2016_sum/total_sum_2016
Shrub_Scrub_2016_domain_frac <- Shrub_Scrub_2016_sum/total_sum_2016
Evergreen_Forest_South_2016_domain_frac <- Evergreen_Forest_South_2016_sum/total_sum_2016

# look at the fractions 
no_data_2016_domain_frac * 100
corn_2016_domain_frac * 100
other_crops_2016_domain_frac * 100
urban_2016_domain_frac * 100
Wetlands_2016_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2016_domain_frac * 100
Deciduous_Forest_2016_domain_frac * 100
Evergreen_Forest_North_2016_domain_frac * 100
Shrub_Scrub_2016_domain_frac * 100
Evergreen_Forest_South_2016_domain_frac * 100

# check that the sum is one 
sum(no_data_2016_domain_frac, 
    corn_2016_domain_frac, 
    other_crops_2016_domain_frac, 
    urban_2016_domain_frac,
    Wetlands_2016_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2016_domain_frac, 
    Deciduous_Forest_2016_domain_frac,
    Evergreen_Forest_North_2016_domain_frac, 
    Shrub_Scrub_2016_domain_frac, 
    Evergreen_Forest_South_2016_domain_frac)


###### look at 2015 
# calculate the sum over all the layers 
vegfrac_layer_sum_2015 <- sum(vegfrac_2015,na.rm=TRUE)
total_sum_2015 <- cellStats(vegfrac_layer_sum_2015,stat='sum')

no_data_2015_sum <- sum(no_data_2015$landcover_percent_2015_1,na.rm=TRUE)  
corn_2015_sum <- sum(corn_2015$landcover_percent_2015_2,na.rm=TRUE)
other_crops_2015_sum <- sum(other_crops_2015$landcover_percent_2015_3,na.rm=TRUE)
urban_2015_sum <- sum(urban_2015$landcover_percent_2015_4,na.rm=TRUE)
Wetlands_2015_sum <- sum(Wetlands_2015$landcover_percent_2015_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2015_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2015$landcover_percent_2015_6,na.rm=TRUE)
Deciduous_Forest_2015_sum <- sum(Deciduous_Forest_2015$landcover_percent_2015_7,na.rm=TRUE)
Evergreen_Forest_North_2015_sum <- sum(Evergreen_Forest_North_2015$landcover_percent_2015_8,na.rm=TRUE)
Shrub_Scrub_2015_sum <- sum(Shrub_Scrub_2015$landcover_percent_2015_9,na.rm=TRUE) 
Evergreen_Forest_South_2015_sum <- sum(Evergreen_Forest_South_2015$landcover_percent_2015_10,na.rm=TRUE) 

no_data_2015_domain_frac <- no_data_2015_sum/total_sum_2015
corn_2015_domain_frac <- corn_2015_sum/total_sum_2015
other_crops_2015_domain_frac <- other_crops_2015_sum/total_sum_2015
urban_2015_domain_frac <- urban_2015_sum/total_sum_2015
Wetlands_2015_domain_frac <- Wetlands_2015_sum/total_sum_2015
Herbaceuous_Grassland_Hay_Pasture_2015_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2015_sum/total_sum_2015
Deciduous_Forest_2015_domain_frac <- Deciduous_Forest_2015_sum/total_sum_2015
Evergreen_Forest_North_2015_domain_frac <- Evergreen_Forest_North_2015_sum/total_sum_2015
Shrub_Scrub_2015_domain_frac <- Shrub_Scrub_2015_sum/total_sum_2015
Evergreen_Forest_South_2015_domain_frac <- Evergreen_Forest_South_2015_sum/total_sum_2015

# look at the fractions 
no_data_2015_domain_frac * 100
corn_2015_domain_frac * 100
other_crops_2015_domain_frac * 100
urban_2015_domain_frac * 100
Wetlands_2015_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2015_domain_frac * 100
Deciduous_Forest_2015_domain_frac * 100
Evergreen_Forest_North_2015_domain_frac * 100
Shrub_Scrub_2015_domain_frac * 100
Evergreen_Forest_South_2015_domain_frac * 100

# check that the sum is one 
sum(no_data_2015_domain_frac, 
    corn_2015_domain_frac, 
    other_crops_2015_domain_frac, 
    urban_2015_domain_frac,
    Wetlands_2015_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2015_domain_frac, 
    Deciduous_Forest_2015_domain_frac,
    Evergreen_Forest_North_2015_domain_frac, 
    Shrub_Scrub_2015_domain_frac, 
    Evergreen_Forest_South_2015_domain_frac)


###### look at 2014 
# calculate the sum over all the layers 
vegfrac_layer_sum_2014 <- sum(vegfrac_2014,na.rm=TRUE)
total_sum_2014 <- cellStats(vegfrac_layer_sum_2014,stat='sum')

no_data_2014_sum <- sum(no_data_2014$landcover_percent_2014_1,na.rm=TRUE)  
corn_2014_sum <- sum(corn_2014$landcover_percent_2014_2,na.rm=TRUE)
other_crops_2014_sum <- sum(other_crops_2014$landcover_percent_2014_3,na.rm=TRUE)
urban_2014_sum <- sum(urban_2014$landcover_percent_2014_4,na.rm=TRUE)
Wetlands_2014_sum <- sum(Wetlands_2014$landcover_percent_2014_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2014_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2014$landcover_percent_2014_6,na.rm=TRUE)
Deciduous_Forest_2014_sum <- sum(Deciduous_Forest_2014$landcover_percent_2014_7,na.rm=TRUE)
Evergreen_Forest_North_2014_sum <- sum(Evergreen_Forest_North_2014$landcover_percent_2014_8,na.rm=TRUE)
Shrub_Scrub_2014_sum <- sum(Shrub_Scrub_2014$landcover_percent_2014_9,na.rm=TRUE) 
Evergreen_Forest_South_2014_sum <- sum(Evergreen_Forest_South_2014$landcover_percent_2014_10,na.rm=TRUE) 

no_data_2014_domain_frac <- no_data_2014_sum/total_sum_2014
corn_2014_domain_frac <- corn_2014_sum/total_sum_2014
other_crops_2014_domain_frac <- other_crops_2014_sum/total_sum_2014
urban_2014_domain_frac <- urban_2014_sum/total_sum_2014
Wetlands_2014_domain_frac <- Wetlands_2014_sum/total_sum_2014
Herbaceuous_Grassland_Hay_Pasture_2014_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2014_sum/total_sum_2014
Deciduous_Forest_2014_domain_frac <- Deciduous_Forest_2014_sum/total_sum_2014
Evergreen_Forest_North_2014_domain_frac <- Evergreen_Forest_North_2014_sum/total_sum_2014
Shrub_Scrub_2014_domain_frac <- Shrub_Scrub_2014_sum/total_sum_2014
Evergreen_Forest_South_2014_domain_frac <- Evergreen_Forest_South_2014_sum/total_sum_2014

# look at the fractions 
no_data_2014_domain_frac * 100
corn_2014_domain_frac * 100
other_crops_2014_domain_frac * 100
urban_2014_domain_frac * 100
Wetlands_2014_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2014_domain_frac * 100
Deciduous_Forest_2014_domain_frac * 100
Evergreen_Forest_North_2014_domain_frac * 100
Shrub_Scrub_2014_domain_frac * 100
Evergreen_Forest_South_2014_domain_frac * 100

# check that the sum is one 
sum(no_data_2014_domain_frac, 
    corn_2014_domain_frac, 
    other_crops_2014_domain_frac, 
    urban_2014_domain_frac,
    Wetlands_2014_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2014_domain_frac, 
    Deciduous_Forest_2014_domain_frac,
    Evergreen_Forest_North_2014_domain_frac, 
    Shrub_Scrub_2014_domain_frac, 
    Evergreen_Forest_South_2014_domain_frac)


###### look at 2013 
# calculate the sum over all the layers 
vegfrac_layer_sum_2013 <- sum(vegfrac_2013,na.rm=TRUE)
total_sum_2013 <- cellStats(vegfrac_layer_sum_2013,stat='sum')

no_data_2013_sum <- sum(no_data_2013$landcover_percent_2013_1,na.rm=TRUE)  
corn_2013_sum <- sum(corn_2013$landcover_percent_2013_2,na.rm=TRUE)
other_crops_2013_sum <- sum(other_crops_2013$landcover_percent_2013_3,na.rm=TRUE)
urban_2013_sum <- sum(urban_2013$landcover_percent_2013_4,na.rm=TRUE)
Wetlands_2013_sum <- sum(Wetlands_2013$landcover_percent_2013_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2013_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2013$landcover_percent_2013_6,na.rm=TRUE)
Deciduous_Forest_2013_sum <- sum(Deciduous_Forest_2013$landcover_percent_2013_7,na.rm=TRUE)
Evergreen_Forest_North_2013_sum <- sum(Evergreen_Forest_North_2013$landcover_percent_2013_8,na.rm=TRUE)
Shrub_Scrub_2013_sum <- sum(Shrub_Scrub_2013$landcover_percent_2013_9,na.rm=TRUE) 
Evergreen_Forest_South_2013_sum <- sum(Evergreen_Forest_South_2013$landcover_percent_2013_10,na.rm=TRUE) 

no_data_2013_domain_frac <- no_data_2013_sum/total_sum_2013
corn_2013_domain_frac <- corn_2013_sum/total_sum_2013
other_crops_2013_domain_frac <- other_crops_2013_sum/total_sum_2013
urban_2013_domain_frac <- urban_2013_sum/total_sum_2013
Wetlands_2013_domain_frac <- Wetlands_2013_sum/total_sum_2013
Herbaceuous_Grassland_Hay_Pasture_2013_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2013_sum/total_sum_2013
Deciduous_Forest_2013_domain_frac <- Deciduous_Forest_2013_sum/total_sum_2013
Evergreen_Forest_North_2013_domain_frac <- Evergreen_Forest_North_2013_sum/total_sum_2013
Shrub_Scrub_2013_domain_frac <- Shrub_Scrub_2013_sum/total_sum_2013
Evergreen_Forest_South_2013_domain_frac <- Evergreen_Forest_South_2013_sum/total_sum_2013

# look at the fractions 
no_data_2013_domain_frac * 100
corn_2013_domain_frac * 100
other_crops_2013_domain_frac * 100
urban_2013_domain_frac * 100
Wetlands_2013_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2013_domain_frac * 100
Deciduous_Forest_2013_domain_frac * 100
Evergreen_Forest_North_2013_domain_frac * 100
Shrub_Scrub_2013_domain_frac * 100
Evergreen_Forest_South_2013_domain_frac * 100

# check that the sum is one 
sum(no_data_2013_domain_frac, 
    corn_2013_domain_frac, 
    other_crops_2013_domain_frac, 
    urban_2013_domain_frac,
    Wetlands_2013_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2013_domain_frac, 
    Deciduous_Forest_2013_domain_frac,
    Evergreen_Forest_North_2013_domain_frac, 
    Shrub_Scrub_2013_domain_frac, 
    Evergreen_Forest_South_2013_domain_frac)


###### look at 2012 
# calculate the sum over all the layers 
vegfrac_layer_sum_2012 <- sum(vegfrac_2012,na.rm=TRUE)
total_sum_2012 <- cellStats(vegfrac_layer_sum_2012,stat='sum')

no_data_2012_sum <- sum(no_data_2012$landcover_percent_2012_1,na.rm=TRUE)  
corn_2012_sum <- sum(corn_2012$landcover_percent_2012_2,na.rm=TRUE)
other_crops_2012_sum <- sum(other_crops_2012$landcover_percent_2012_3,na.rm=TRUE)
urban_2012_sum <- sum(urban_2012$landcover_percent_2012_4,na.rm=TRUE)
Wetlands_2012_sum <- sum(Wetlands_2012$landcover_percent_2012_5,na.rm=TRUE)
Herbaceuous_Grassland_Hay_Pasture_2012_sum <- sum(Herbaceuous_Grassland_Hay_Pasture_2012$landcover_percent_2012_6,na.rm=TRUE)
Deciduous_Forest_2012_sum <- sum(Deciduous_Forest_2012$landcover_percent_2012_7,na.rm=TRUE)
Evergreen_Forest_North_2012_sum <- sum(Evergreen_Forest_North_2012$landcover_percent_2012_8,na.rm=TRUE)
Shrub_Scrub_2012_sum <- sum(Shrub_Scrub_2012$landcover_percent_2012_9,na.rm=TRUE) 
Evergreen_Forest_South_2012_sum <- sum(Evergreen_Forest_South_2012$landcover_percent_2012_10,na.rm=TRUE) 

no_data_2012_domain_frac <- no_data_2012_sum/total_sum_2012
corn_2012_domain_frac <- corn_2012_sum/total_sum_2012
other_crops_2012_domain_frac <- other_crops_2012_sum/total_sum_2012
urban_2012_domain_frac <- urban_2012_sum/total_sum_2012
Wetlands_2012_domain_frac <- Wetlands_2012_sum/total_sum_2012
Herbaceuous_Grassland_Hay_Pasture_2012_domain_frac <- Herbaceuous_Grassland_Hay_Pasture_2012_sum/total_sum_2012
Deciduous_Forest_2012_domain_frac <- Deciduous_Forest_2012_sum/total_sum_2012
Evergreen_Forest_North_2012_domain_frac <- Evergreen_Forest_North_2012_sum/total_sum_2012
Shrub_Scrub_2012_domain_frac <- Shrub_Scrub_2012_sum/total_sum_2012
Evergreen_Forest_South_2012_domain_frac <- Evergreen_Forest_South_2012_sum/total_sum_2012

# look at the fractions 
no_data_2012_domain_frac * 100
corn_2012_domain_frac * 100
other_crops_2012_domain_frac * 100
urban_2012_domain_frac * 100
Wetlands_2012_domain_frac * 100
Herbaceuous_Grassland_Hay_Pasture_2012_domain_frac * 100
Deciduous_Forest_2012_domain_frac * 100
Evergreen_Forest_North_2012_domain_frac * 100
Shrub_Scrub_2012_domain_frac * 100
Evergreen_Forest_South_2012_domain_frac * 100

# check that the sum is one 
sum(no_data_2012_domain_frac, 
    corn_2012_domain_frac, 
    other_crops_2012_domain_frac, 
    urban_2012_domain_frac,
    Wetlands_2012_domain_frac,
    Herbaceuous_Grassland_Hay_Pasture_2012_domain_frac, 
    Deciduous_Forest_2012_domain_frac,
    Evergreen_Forest_North_2012_domain_frac, 
    Shrub_Scrub_2012_domain_frac, 
    Evergreen_Forest_South_2012_domain_frac)

########### make a dataframe to store all the fractions ###########
year <- c(2012,2013,2014,2015,2016,2017,2018,2019,2020,2021)
no_data <- c(no_data_2012_domain_frac,no_data_2013_domain_frac,no_data_2014_domain_frac, 
             no_data_2015_domain_frac,no_data_2016_domain_frac,no_data_2017_domain_frac,
             no_data_2018_domain_frac,no_data_2019_domain_frac,no_data_2020_domain_frac,
             no_data_2021_domain_frac)
corn <- c(corn_2012_domain_frac,corn_2013_domain_frac,corn_2014_domain_frac,
          corn_2015_domain_frac,corn_2016_domain_frac,corn_2017_domain_frac,
          corn_2018_domain_frac,corn_2019_domain_frac,corn_2020_domain_frac,
          corn_2021_domain_frac)
other_crops <- c(other_crops_2012_domain_frac,other_crops_2013_domain_frac,
                 other_crops_2014_domain_frac,other_crops_2015_domain_frac,
                 other_crops_2016_domain_frac,other_crops_2017_domain_frac,
                 other_crops_2018_domain_frac,other_crops_2019_domain_frac,
                 other_crops_2020_domain_frac,other_crops_2021_domain_frac)
urban <- c(urban_2012_domain_frac,urban_2013_domain_frac,urban_2014_domain_frac,
           urban_2015_domain_frac,urban_2016_domain_frac,urban_2017_domain_frac,
           urban_2018_domain_frac,urban_2019_domain_frac,urban_2020_domain_frac,
           urban_2021_domain_frac)
wetlands <- c(Wetlands_2012_domain_frac,Wetlands_2013_domain_frac,
              Wetlands_2014_domain_frac,Wetlands_2015_domain_frac,
              Wetlands_2016_domain_frac,Wetlands_2017_domain_frac,
              Wetlands_2018_domain_frac,Wetlands_2019_domain_frac,
              Wetlands_2020_domain_frac,Wetlands_2021_domain_frac)
Herbaceuous_Grassland_Hay_Pasture <- c(Herbaceuous_Grassland_Hay_Pasture_2012_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2013_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2014_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2015_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2016_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2017_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2018_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2019_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2020_domain_frac,
                                       Herbaceuous_Grassland_Hay_Pasture_2021_domain_frac)
dbf <- c(Deciduous_Forest_2012_domain_frac,Deciduous_Forest_2013_domain_frac,
         Deciduous_Forest_2014_domain_frac,Deciduous_Forest_2015_domain_frac,
         Deciduous_Forest_2016_domain_frac,Deciduous_Forest_2017_domain_frac,
         Deciduous_Forest_2018_domain_frac,Deciduous_Forest_2019_domain_frac,
         Deciduous_Forest_2020_domain_frac,Deciduous_Forest_2021_domain_frac)
Evergreen_North <- c(Evergreen_Forest_North_2012_domain_frac,Evergreen_Forest_North_2013_domain_frac,
                     Evergreen_Forest_North_2014_domain_frac,Evergreen_Forest_North_2015_domain_frac,
                     Evergreen_Forest_North_2016_domain_frac,Evergreen_Forest_North_2017_domain_frac,
                     Evergreen_Forest_North_2018_domain_frac,Evergreen_Forest_North_2019_domain_frac,
                     Evergreen_Forest_North_2020_domain_frac,Evergreen_Forest_North_2021_domain_frac)
Shrub_Scrub <- c(Shrub_Scrub_2012_domain_frac,Shrub_Scrub_2013_domain_frac,
                 Shrub_Scrub_2014_domain_frac,Shrub_Scrub_2015_domain_frac,
                 Shrub_Scrub_2016_domain_frac,Shrub_Scrub_2017_domain_frac,
                 Shrub_Scrub_2018_domain_frac,Shrub_Scrub_2019_domain_frac,
                 Shrub_Scrub_2020_domain_frac,Shrub_Scrub_2021_domain_frac)
Evergreen_South <- c(Evergreen_Forest_South_2012_domain_frac,Evergreen_Forest_South_2013_domain_frac,
                     Evergreen_Forest_South_2014_domain_frac,Evergreen_Forest_South_2015_domain_frac,
                     Evergreen_Forest_South_2016_domain_frac,Evergreen_Forest_South_2017_domain_frac,
                     Evergreen_Forest_South_2018_domain_frac,Evergreen_Forest_South_2019_domain_frac,
                     Evergreen_Forest_South_2020_domain_frac,Evergreen_Forest_South_2021_domain_frac)

PFT_fractions_df <- data.frame(year,no_data,corn,other_crops,urban,wetlands,
                               Herbaceuous_Grassland_Hay_Pasture,dbf,
                               Evergreen_North,Shrub_Scrub,Evergreen_South)

# save out fractions 
write.csv(PFT_fractions_df,'/storage/group/zrb5027/default/INFLUX/smm8236/vprm_gridded_analysis/PFT_yearly_fractions.csv')
