%%% script written by Sam Murphy 
%%% purpose: use the flux footprint model (FFP) from Kljun et al 2015 to
%%% calculate flux footprints for use in ag site convolutions with crop
%%% maps - this is for site ANWa in 2018 when there was the corn and soy split

clear all; 
close all; 
clc; 

%% set input data paths and read in data 

% set path to the flux data set of interest 
path_flux_data = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_data/data_with_flags_removed/'; 

% set the name of the data file for the site/time period of interest 
file_name = 'PROCESSED_SITE_AGa_FLAGS_REMOVED.csv'; 

% set the name of the site (match the input data file)
site = 'site_NWa_US-INd_2018_Corn_Soy_split'; 

% set the path to store the output csvs
path_out = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_data/data_with_flags_removed_and_veg_fraction/'; 

% load in the flux data for the site/time period of interest 
path_flux_data_file = strcat(path_flux_data,file_name); 

% read in the csv file - read in as a matrix 
data_in = readmatrix(path_flux_data_file); 

% read in the csv file as a table 
data_in_table = readtable(path_flux_data_file); 

%% make datetime for data time series %% 

% read date 
dates_time_series(:,1) = data_in_table.index; 

% shift the time series - flux data is marked with the end of the half hour
% averaging period, but we want the timestamp to be at the start 
dates_time_series = dates_time_series - minutes(30); 


%% read in the reanlysis data for the ABL heights 
% also interpolate from an hourly temporal resolution to half hourly 

% path to the directory that stores the abl heights file 
abl_data_path = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_footprint_filtering/ABL_reanalysis/';

% read in the abl heights as a table 
abl_data = readtable(strcat(abl_data_path,'INFLUX_ABL_heights.csv')); 

% grab the column name for the site 
site_abl_heights = abl_data.("ABL_Height_USINd"); 

% grab the times of the data 
%all_time_steps_abl = datetime(string(table2cell(abl_data(:,1))),'InputFormat','MM/dd/yyyy HH:mm');
all_time_steps_abl = abl_data.TmStamp; 

% check that the abl values are hourly (with no gaps)
length(datetime(2012,1,1,0,0,0):hours(1):datetime(2022,12,31,23,0,0))
length(all_time_steps_abl)

% make a timetable with the hourly abl values 
abl_hrly_timetable = timetable(all_time_steps_abl,site_abl_heights); 

% create a half hour datetime vector for same time period as data - add the
% 30 at the end to bring the time series to dec 31 23:30:00 (otherwise it
% is dec 31 23:00:00
time_steps_30min(:,1) = all_time_steps_abl(1):minutes(30):(all_time_steps_abl(end)+minutes(30));

% use retime to interpolate the abl values linearly to a half hourly time step 
abl_30min_timetable = retime(abl_hrly_timetable,time_steps_30min,'linear'); 

% make the new half hourly values into a vector 
abl_30min(:,1) = abl_30min_timetable.site_abl_heights; 

%% grab inputs for the model %%%

% friction velocity (ustar)
ustar = data_in_table.('u_');

% horizontal wind speed 
horz_wind_speed = data_in_table.('u_rot');

% input the measurement height - this might change for some sites so check
% before running - check for 2021 sites*
measure_height = 3; 

% variance of u wind speed - used to calculate the standard dev of lateral
% velocity fluctuations 
uvar = data_in_table.('u_var'); 

% standard dev of lateral velocity fluctuations (standard dev = sqrt of variance of u)
u_std = sqrt(uvar(:)); 

% M-O length - check that this is the right column 
MO_length = data_in_table.('L'); 

% wind direction(in degrees)
wind_dir = data_in_table.('wind_dir'); 

% need to grab the abl heights AT THE TIME THE SITE DATA STARTS 
dates_time_series(1)
dates_time_series(end)

ind_BL_start = find(time_steps_30min==dates_time_series(1)); 
ind_BL_end = find(time_steps_30min==dates_time_series(end)); 

time_steps_30min(ind_BL_start)
time_steps_30min(ind_BL_end)

abl_height = abl_30min(ind_BL_start:ind_BL_end);

%% load in the site map matrix - corn %%

% set path to the site maps in csv format 
site_map_path = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_footprint_filtering/site_map_csvs/';

% set the file name of the site to run
site_map_file = 'site_A_NWa_2018_CORN_matrix.csv'; 

% load in the site map as a table 
site_map = readtable(strcat(site_map_path,site_map_file)); 

% convert the table to an array (double) 
site_map = table2array(site_map); 

% delete row label column 
site_map(:,1) = []; 

%%  run footprint model and multiply with the CORN site map %%

% make an empty matrix to store the map convolution 
map_conv(:,:) = NaN(501,501); 

% make empty vector to store the sum of the footprints for each half hour
footprint_sum(:,1) = NaN(length(data_in),1);

% make an empty vector to store the convolution between the flux footprint
% and the site map 
map_conv_sum(:,1) = NaN(length(data_in),1);

% make an empty vector to store the fraction of the footprint attributed to
% the vegetation of interest 
frac_veg(:,1) = NaN(length(data_in),1);

% make an empty vector to store the footprint error flag
flag_err_all(:,1) = NaN(length(data_in),1);

% loop over all the half hours to calculate footprints, and the fraction of
% hte footprint that is the PFT/vegetation of interest 
for i=1:length(data_in(:,1))
    
    % use climatology function to caclulate the footprint for half hour i 
    [FFP,flag_err] = calc_footprint_FFP_climatology(measure_height,NaN,horz_wind_speed(i),abl_height(i),MO_length(i),u_std(i),ustar(i),wind_dir(i),'dx',1,'dy',1,'domain',[-250 250 -250 250]); 

    % if there is no error in the footprint 
    if flag_err ~= 1 

        % multiply the flipped footprint (need to do this to orient the 
        % matrix correctly) by the site map 
        map_conv(:,:) = flip(FFP(1).fclim_2d) .* site_map;   
    
        % sum over the regular footprint 
        footprint_sum(i) = sum(FFP(1).fclim_2d, 'all'); 

        % sum over the footprint that is multiplied by the site map
        map_conv_sum(i) = sum(map_conv(:,:),'all'); 

    end     

    % calculate the fraction of the footprint attributable to the
    % vegetation of interest 
    frac_veg(i) = map_conv_sum(i)/footprint_sum(i); 

    % save the footprint function error flag here 
    flag_err_all(i) = flag_err; 

    clear FFP flag_err map_conv
end 



%% attach the fraction of corn to the data table 

% convert the vegetation fractions and the flags to a table
frac_veg_table = array2table(horzcat(frac_veg,flag_err_all),"VariableNames",{'Frac Corn','Footprint Flag Corn'}); 

% make a string of dates 
dates_string = string(data_in_table.index,"yyyy-MM-dd HH:mm:ss"); 

dates_string_table = table(dates_string);  

% add the fraction of vegetation and flags to the co2 flux table + add a
% date string 
data_out = [data_in_table dates_string_table frac_veg_table]; 


%% now read in map to run for soy %%

% clear variables before running again 
clear site_map site_map_file flag_veg map_conv_sum map_conv flag_err_all footprint_sum

% set path to the site maps in csv format 
site_map_path = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_footprint_filtering/site_map_csvs/';

% set the file name of the site to run 
site_map_file = 'site_A_NWa_2018_SOY_matrix.csv';

% load in the site map as a table 
site_map = readtable(strcat(site_map_path,site_map_file)); 

% convert the table to an array (double) 
site_map = table2array(site_map); 

% delete row label column 
site_map(:,1) = []; 

%% run footprint model and multiply with the SOY site map %%

% make an empty matrix to store the map convolution 
map_conv(:,:) = NaN(501,501); 

% make empty vector to store the sum of the footprints for each half hour
footprint_sum(:,1) = NaN(length(data_in),1);

% make an empty vector to store the convolution between the flux footprint
% and the site map 
map_conv_sum(:,1) = NaN(length(data_in),1);

% make an empty vector to store the fraction of the footprint attributed to
% the vegetation of interest 
frac_veg(:,1) = NaN(length(data_in),1);

% make an empty vector to store the footprint error flag
flag_err_all(:,1) = NaN(length(data_in),1);

% loop over all the half hours to calculate footprints, and the fraction of
% hte footprint that is the PFT/vegetation of interest 
for i=1:length(data_in(:,1))
    
    % use climatology function to caclulate the footprint for half hour i 
    [FFP,flag_err] = calc_footprint_FFP_climatology(measure_height,NaN,horz_wind_speed(i),abl_height(i),MO_length(i),u_std(i),ustar(i),wind_dir(i),'dx',1,'dy',1,'domain',[-250 250 -250 250]); 

    % if there is no error in the footprint 
    if flag_err ~= 1 

        % multiply the flipped footprint (need to do this to orient the 
        % matrix correctly) by the site map 
        map_conv(:,:) = flip(FFP(1).fclim_2d) .* site_map;   
    
        % sum over the regular footprint 
        footprint_sum(i) = sum(FFP(1).fclim_2d, 'all'); 

        % sum over the footprint that is multiplied by the site map
        map_conv_sum(i) = sum(map_conv(:,:),'all'); 

    end     

    % calculate the fraction of the footprint attributable to the
    % vegetation of interest 
    frac_veg(i) = map_conv_sum(i)/footprint_sum(i); 

    % save the footprint function error flag here 
    flag_err_all(i) = flag_err; 

    clear FFP flag_err map_conv
end 



%% attach the fraction of soy to the data table 

% convert the vegetation fractions and the flags to a table
frac_veg_table = array2table(horzcat(frac_veg,flag_err_all),"VariableNames",{'Frac Soy','Footprint Flag Soy'}); 

% add the fraction of vegetation and flags to the co2 flux table
data_out = [data_out frac_veg_table]; 

%% save file with both fractions 

% write out the data file 
writetable(data_out,strcat(path_out,site,'.csv'));



