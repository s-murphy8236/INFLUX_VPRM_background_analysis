%%%%% written by: Sam Murphy
%%%%% purpose: model-obs comparisons / analysis of model performance for my 
%%%%% modeled biogenic mole fractions at the influx background towers (01, 
%%%%% 09, 14)

% clean up 
clear all 
close all 

%% read in & quick plots of mole frac obs at towers 01, 09, 14

% set path to the mole fraction tower observations (on roar collab)
data_path = ['/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/obspack_data' ...
    '/obspack_multi-species_1_INFLUX_v2.0_2022-03-02/data/nc/']; 

%%%%%%%%%%%%%%%%%%%%%%%%%% load in tower 01 data %%%%%%%%%%%%%%%%%%%%%%%%%%%
% tower 1 co2 concentration at 121 meters agl (highest level on tower)

% set path to tower 1 data file 
co2_t01_121magl_path = horzcat(data_path,'co2_inx01_surface-insitu_60_allhours-121magl.nc'); 

% get info about the data 
ncdisp(co2_t01_121magl_path)

% read in the mole fraction of co2 - this is the variable value, which is
% called  'measured_mole_fraction_of_trace_gas_in_dry_air' in units of
t01_121magl_co2 = ncread(co2_t01_121magl_path,'value'); 

% get the time stamps for the data - this is the variable time, which is 
% called sample_time_in_seconds_since_january_1_1970, in seconds since 
% 1970-01-01T00:00:00Z'. POSIX time. Number of seconds since January 1, 1970 
% in UTC. Time-averaged values are reported at the middle of the averaging 
% interval. For uneven averaging intervals times are rounded down to the nearest second.'
t01_121magl_time_posix = ncread(co2_t01_121magl_path,'time'); 

% convert the dates to datetime format 
t01_121magl_time = datetime(t01_121magl_time_posix,'ConvertFrom','posixtime'); 

% look at the start and end times of the data 
t01_121magl_time(1)
t01_121magl_time(end)

% convert units from mol/mol to ppm 
t01_121magl_co2_ppm = t01_121magl_co2*(1*10^6); 

%%%%%%%%%%%%%%%%%%%%%%%%%% load in tower 09 data %%%%%%%%%%%%%%%%%%%%%%%%%%
% tower 9 co2 concentration at 130 meters agl (highest level on tower)

% set path to tower 14 data file 
co2_t09_130magl_path = horzcat(data_path,'co2_inx09_surface-insitu_60_allhours-130magl.nc'); 

% get info about the data 
ncdisp(co2_t09_130magl_path)

% read in the mole fraction of co2 - this is the variable value, which is
% called  'measured_mole_fraction_of_trace_gas_in_dry_air' in units of
% 'mole per mole of dry air'
t09_130magl_co2 = ncread(co2_t09_130magl_path,'value'); 

% note: _FillValue    = -9.999999790214768e+33 - need to check if there are
% any of these values? 

% get the time stamps for the data - POSIX time
t09_130magl_time_posix = ncread(co2_t09_130magl_path,'time'); 

% convert the dates to datetime format 
t09_130magl_time = datetime(t09_130magl_time_posix,'ConvertFrom','posixtime'); 

% look at the start and end times of the data 
t09_130magl_time(1)
t09_130magl_time(end)

% convert units from mol/mol to ppm 
t09_130magl_co2_ppm = t09_130magl_co2*(1*10^6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%% load in tower 14 data %%%%%%%%%%%%%%%%%%%%%%%%%%
% tower 14 co2 concentration at 76 meters agl (highest level on tower)

% set path to tower 14 data file 
co2_t14_76magl_path = horzcat(data_path,'co2_inx14_surface-insitu_60_allhours-76magl.nc'); 

% get info about the data 
ncdisp(co2_t14_76magl_path)

% read in the mole fraction of co2 - this is the variable value, which is
% called  'measured_mole_fraction_of_trace_gas_in_dry_air' in units of
% 'mole per mole of dry air'
t14_76magl_co2 = ncread(co2_t14_76magl_path,'value'); 

% get the time stamps for the data - POSIX time
t14_76magl_time_posix = ncread(co2_t14_76magl_path,'time'); 

% convert the dates to datetime format 
t14_76magl_time = datetime(t14_76magl_time_posix,'ConvertFrom','posixtime'); 

% look at the start and end times of the data 
t14_76magl_time(1)
t14_76magl_time(end)

% convert units from mol/mol to ppm 
t14_76magl_co2_ppm = t14_76magl_co2*(1*10^6); 


%%%%%%%%%%%%%%%%% quick plots of all data for each site %%%%%%%%%%%%%%%%%%%%% 

% tower 9 with units mol/mol
figure()
t = tiledlayout(3,1); 
nexttile 
plot(t01_121magl_time,t01_121magl_co2,'.'); 
xlabel('time [UTC]'); 
title('Tower 01')
nexttile
plot(t09_130magl_time,t09_130magl_co2,'.'); 
xlabel('time [UTC]'); 
title('Tower 09')
nexttile
plot(t14_76magl_time,t14_76magl_co2,'.'); 
xlabel('time [UTC]'); 
title('Tower 14')
ylabel(t,'CO2 concentration [mol/mol]')


% tower 1 
figure()
plot(t01_121magl_time,t01_121magl_co2_ppm,'.'); 
xlabel('time [UTC]'); 
ylabel('CO2 concentration [ppm]')
title('Tower 1');

% tower 9 with units pmm 
figure()
plot(t09_130magl_time,t09_130magl_co2_ppm,'.'); 
xlabel('time [UTC]'); 
ylabel('CO2 concentration [ppm]')
title('Tower 9');

% tower 14 
figure()
plot(t14_76magl_time,t14_76magl_co2_ppm,'.'); 
xlabel('time [UTC]'); 
ylabel('CO2 concentration [ppm]')
title('Tower 14');

%% "fill" the missing obs times with NaN values
% fill the missing times with Nans because time stamps in the data file are 
% skipped if there are no data available 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tower 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% get the start and end times of the data for tower 1
t01_121magl_time(1)
t01_121magl_time(end)

% make a datetime vector for the entire time length of the data - change
% to be the end of 2020 to make things easier 
%t01_all_time(:,1) = t01_121magl_time(1):hours(1):t01_121magl_time(end);
t01_all_time(:,1) = t01_121magl_time(1):hours(1):datetime(2020,12,31,23,30,0);

% make a NaN vector for the entire time length of the tower 1 data 
t01_co2 = NaN(length(t01_all_time),1); 

% add data values to our vector of nans (which is the correct length for
% the time span of the data set 
for i = 1:length(t01_121magl_time) % for all values in the data file

    % find the corresponding time index in the new vector 
    ind = find(t01_all_time == t01_121magl_time(i)); 

    % save the co2 value to the new hourly co2 vector
    t01_co2(ind) = t01_121magl_co2_ppm(i); 

    clear ind % clean up for the next iteration 
end 
clear i % clean up for later loops 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tower 09 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the start and end times of the tower 9 data 
t09_130magl_time(1)
t09_130magl_time(end)

% make a datetime vector for the entire time length of the data 
t09_all_time(:,1) = t09_130magl_time(1):hours(1):t09_130magl_time(end);

% make a NaN vector for the entire time length of the tower 9 data 
t09_co2 = NaN(length(t09_all_time),1); 

% add data values to our vector of nans (which is the correct length for
% the time span of the data set)
for i = 1:length(t09_130magl_time) % for all values in the data file
    
    % find the corresponding time index in the new vector 
    ind = find(t09_all_time == t09_130magl_time(i));  

    % save the co2 value to the new hourly co2 vector
    t09_co2(ind) = t09_130magl_co2_ppm(i); 

    clear ind % clean up for the next iteration
end 
clear i % clean up for later loops 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tower 14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the start and end times of the data for tower 14
t14_76magl_time(1)
t14_76magl_time(end)

% make a datetime vector for the entire time length of the data
t14_all_time(:,1) = t14_76magl_time(1):hours(1):t14_76magl_time(end);

% make a NaN vector for the entire time length of the tower 14 data 
t14_co2 = NaN(length(t14_all_time),1); 

% add data values to our vector of nans (which is the correct length for
% the time span of the data set)
for i = 1:length(t14_76magl_time)

    % find the corresponding time index in the new vector 
    ind = find(t14_all_time == t14_76magl_time(i)); 

    % save the co2 value to the new hourly co2 vector
    t14_co2(ind) = t14_76magl_co2_ppm(i); 

    clear ind % clean up for the next iteratio
end 
clear i % clean up for later loops 

%% add NaNs to the start of 2017 for Tower 14 %% 
% now, fill the start of 2017 with NaNs for tower 14 - this is just to make
% it easier for analysis later (tower 14 starts in april 2017)

% make datetime vector for the start of the year to the start of the data
% start on the half hour to to match the data 
year_start_t14(:,1) = datetime(2017,1,1,0,30,0):hours(1):t14_76magl_time(1);

% make a NaN vectory for the length of the start of the year before data is
% collected - need to subtract 1 to avoid double counting the first hour of
% data
NaN_t14_year_start = NaN([length(year_start_t14(:,1))-1,1]); 

% make a datetime vector going from the start of 2017 to the dataset end 
% start on the half hour to to match the data 
t14_all_time_ADJUSTED(:,1) = datetime(2017,1,1,0,30,0):hours(1):t14_76magl_time(end);

% now make a new vector of co2 values for tower 14 with the NaNs at the 
% start of year appended to it 
t14_co2_NEW = vertcat(NaN_t14_year_start,t14_co2); 

% clear the "old" t14_co2 time series 
clear t14_co2

% now save the t14 co2 time series with the nans at the start of the year
% as t14_co2 for use 
t14_co2 = t14_co2_NEW; 

% now save the new t14 datetime time series (with start of the 2017
% included) with the old datetime name 
clear t14_all_time
t14_all_time = t14_all_time_ADJUSTED; 

% now let's clean up the vectors we aren't going to use anymore 
clear t14_76magl_time t14_76magl_co2_ppm t01_121magl_time t01_121magl_co2_ppm
clear t09_130magl_time t09_130magl_co2_ppm t14_76magl_time_posix t09_130magl_time_posix
clear t01_121magl_time_posix t14_76magl_co2 t01_121magl_co2 t09_130magl_co2
clear t14_co2_NEW t14_all_time_ADJUSTED NaN_t14_year_start year_start_t14
clear co2_t01_121magl_path co2_t09_130magl_path co2_t14_76magl_path


%% wind direction filtering %% 
% remove data when the towers are in the urban plume (directions removed
% come from Miles et al., 2021)

%%%%%%%%%%%%%%%%%%% load in the wind direction data %%%%%%%%%%%%%%%%%%%%%%%
% set path to the wind direction data for each year 
wind_data_2017_path = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/wind_data_processed/hourly_wind_direction_NEWFILL2017.csv'; 
wind_data_2018_path = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/wind_data_processed/hourly_wind_direction_NEWFILL2018.csv'; 
wind_data_2019_path = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/wind_data_processed/hourly_wind_direction_NEWFILL2019.csv'; 
wind_data_2020_path = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/wind_data_processed/hourly_wind_direction_NEWFILL2020.csv';

% read in the wind data matrix 
wind_data_2017=readmatrix(wind_data_2017_path);
wind_data_2018=readmatrix(wind_data_2018_path); 
wind_data_2019=readmatrix(wind_data_2019_path); 
wind_data_2020=readmatrix(wind_data_2020_path); 

%%%%%%%%%% grab yearly subsets of the co2 obs %%%%%%%%%% 
%%% 2017 %%% 
t14_2017_ind = find(year(t14_all_time)==2017); % get indicies of tower 14 2017 data
t14_co2_2017 = t14_co2(t14_2017_ind);          % grab co2 for 2017 
t14_dates_2017 = t14_all_time(t14_2017_ind);   % get just 2017 dates (mostly for checking indices are correct)

t09_2017_ind = find(year(t09_all_time)==2017); 
t09_co2_2017 = t09_co2(t09_2017_ind); 
t09_dates_2017 = t09_all_time(t09_2017_ind);

t01_2017_ind = find(year(t01_all_time)==2017); 
t01_co2_2017 = t01_co2(t01_2017_ind); 
t01_dates_2017 = t01_all_time(t01_2017_ind);

%%% 2018 %%%  
t14_2018_ind = find(year(t14_all_time)==2018); 
t14_co2_2018 = t14_co2(t14_2018_ind); 
t14_dates_2018 = t14_all_time(t14_2018_ind);

t09_2018_ind = find(year(t09_all_time)==2018); 
t09_co2_2018 = t09_co2(t09_2018_ind); 
t09_dates_2018 = t09_all_time(t09_2018_ind);

t01_2018_ind = find(year(t01_all_time)==2018); 
t01_co2_2018 = t01_co2(t01_2018_ind); 
t01_dates_2018 = t01_all_time(t01_2018_ind);

%%% 2019 %%%  
t14_2019_ind = find(year(t14_all_time)==2019); 
t14_co2_2019 = t14_co2(t14_2019_ind); 
t14_dates_2019 = t14_all_time(t14_2019_ind);

t09_2019_ind = find(year(t09_all_time)==2019); 
t09_co2_2019 = t09_co2(t09_2019_ind); 
t09_dates_2019 = t09_all_time(t09_2019_ind);

t01_2019_ind = find(year(t01_all_time)==2019); 
t01_co2_2019 = t01_co2(t01_2019_ind); 
t01_dates_2019 = t01_all_time(t01_2019_ind);

%%% 2020 %%%  
t14_2020_ind = find(year(t14_all_time)==2020); 
t14_co2_2020 = t14_co2(t14_2020_ind); 
t14_dates_2020 = t14_all_time(t14_2020_ind);

t09_2020_ind = find(year(t09_all_time)==2020); 
t09_co2_2020 = t09_co2(t09_2020_ind); 
t09_dates_2020 = t09_all_time(t09_2020_ind);

t01_2020_ind = find(year(t01_all_time)==2020); 
t01_co2_2020 = t01_co2(t01_2020_ind); 
t01_dates_2020 = t01_all_time(t01_2020_ind);

%%%%%%%%%%%% filter the 2018 co2 data by wind direction %%%%%%%%%%%%%%%%%%%

%%%% tower 01 filtering %%%% 
% make an empty vector to store the points that get removed from the filter
t01_co2_windir_removed_2018 = NaN(length(t01_co2_2018),1); 

% store the unfiltered original co2 for t01 (for reference)
t01_co2_no_winddirfilt_2018 = t01_co2_2018; 

% loop through the tower 1 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 20-65 deg)
for i=1:length(t01_co2_2018) % for all of the 2018 data points 
    % if the wind direction is between 20 and 65 deg 
    if wind_data_2018(i,7) >=20 & wind_data_2018(i,7) <= 65
        % save out the filtered point to the new vector
        t01_co2_windir_removed_2018(i) = t01_co2_2018(i); 

        % nan out the point in the data vector
        t01_co2_2018(i) = NaN; 
    end 
end

% loop through tower 1 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t01_co2_2018)
    if isnan(wind_data_2018(i,7))==true
        % save out the removed point to the new vector
        t01_co2_windir_removed_2018(i) =t01_co2_2018(i); 
        % nan out the point in the data vector 
        t01_co2_2018(i) = NaN; 
    end 
end 


%%%% tower 09 filtering %%%% 

% make an empty vector to store the points that get removed from the filter
t09_co2_windir_removed_2018 = NaN(length(t09_co2_2018),1); 

% store the unfiltered original co2 for t09 (for reference)
t09_co2_no_winddirfilt_2018 = t09_co2_2018; 

% loop through the tower 9 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 235-280 deg)
for i=1:length(t09_co2_2018)
     % if the wind direction is between 235 and 280 deg 
    if wind_data_2018(i,7) >=235 & wind_data_2018(i,7) <= 280 
        % save out the filtered point to the new vector
        t09_co2_windir_removed_2018(i) = t09_co2_2018(i); 
        % nan out the point in the data vector 
        t09_co2_2018(i) = NaN; 
    end 
end 

% loop through tower 9 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t09_co2_2018)
    if isnan(wind_data_2018(i,7))==true
        % save out the removed point to the new vector
        t09_co2_windir_removed_2018(i) =t09_co2_2018(i); 
        % nan out the point in the data vector
        t09_co2_2018(i) = NaN;
    end 
end 

%%%%%%%%%%%% filter the 2019 co2 data by wind direction %%%%%%%%%%%%%%%%%%%

%%%% tower 01 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t01_co2_windir_removed_2019 = NaN(length(t01_co2_2019),1); 

% make a new vector to store the unfiltered (by wind dir) co2 for t01 (for reference)
t01_co2_no_winddirfilt_2019 = t01_co2_2019; 

% loop through the tower 1 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 20-65 deg)
for i=1:length(t01_co2_2019)
    % if the wind direction is between 20 and 65 deg 
    if wind_data_2019(i,7) >=20 & wind_data_2019(i,7) <= 65 
        % save out the filtered point to the new vector
        t01_co2_windir_removed_2019(i) = t01_co2_2019(i); 
        % nan out the removed point 
        t01_co2_2019(i) = NaN; 
    end 
end 

% loop through tower 1 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t01_co2_2019)
    if isnan(wind_data_2019(i,7))==true
        % save out the removed point to the new vector
        t01_co2_windir_removed_2019(i) = t01_co2_2019(i); 
        % nan out the point in the data vector
        t01_co2_2019(i) = NaN; 
    end 
end 


%%%% tower 09 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t09_co2_windir_removed_2019 = NaN(length(t09_co2_2019),1); 

% make a new vector to store the unfiltered (by wind dir) co2 for t01 (for reference)
t09_co2_no_winddirfilt_2019 = t09_co2_2019; 

% loop through the tower 9 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 235-280 deg)
for i=1:length(t09_co2_2019)
    % if the wind direction is between 235 and 280 deg 
    if wind_data_2019(i,7) >=235 & wind_data_2019(i,7) <= 280 
        % save out the filtered point to the new vector
        t09_co2_windir_removed_2019(i) =t09_co2_2019(i); 
        % nan out the point in the data vector 
        t09_co2_2019(i) = NaN; 
    end 
end 

% loop through tower 9 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t09_co2_2019)
    if isnan(wind_data_2019(i,7))==true
        % save out the removed point to the new vector
        t09_co2_windir_removed_2019(i) =t09_co2_2019(i); 
        % nan out the point in the data vector
        t09_co2_2019(i) = NaN; 
    end 
end 


%%%%%%%%%%%% filter the 2017 co2 data by wind direction %%%%%%%%%%%%%%%%%%%

%%%% tower 01 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t01_co2_windir_removed_2017 = NaN(length(t01_co2_2017),1); 

% make a new vector to store the unfiltered co2 for t01 (for reference)
t01_co2_no_winddirfilt_2017 = t01_co2_2017; 

% loop through the tower 1 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 20-65 deg)
for i=1:length(t01_co2_2017)
    % if the wind direction is between 20 and 65 deg
    if wind_data_2017(i,7) >=20 & wind_data_2017(i,7) <= 65 
        % save out the filtered point to the new vector
        t01_co2_windir_removed_2017(i) =t01_co2_2017(i);
        % nan out the filtered point 
        t01_co2_2017(i) = NaN; 
    end 
end 


% loop through tower 1 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t01_co2_2017)
    if isnan(wind_data_2017(i,7))==true
        % save out the removed point to the new vector
        t01_co2_windir_removed_2017(i) = t01_co2_2017(i); 
        % nan out the removed points 
        t01_co2_2017(i) = NaN; 
    end 
end 

%%%% tower 09 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t09_co2_windir_removed_2017 = NaN(length(t09_co2_2017),1); 

% make a new vector to store the unfiltered co2 for t09 (for reference)
t09_co2_no_winddirfilt_2017 = t09_co2_2017; 

% loop through the tower 9 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 235-280 deg)
for i=1:length(t09_co2_2017)
    % if the wind direction is between 235 and 280 deg
    if wind_data_2017(i,7) >=235 & wind_data_2017(i,7) <= 280 
        % save out the filtered point to the new vector
        t09_co2_windir_removed_2017(i) = t09_co2_2017(i); 
        % nan out the point in the data vector 
        t09_co2_2017(i) = NaN; 
    end 
end 

% loop through tower 9 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t09_co2_2017)
    if isnan(wind_data_2017(i,7))==true
        % save out the removed point to the new vector
        t09_co2_windir_removed_2017(i) =t09_co2_2017(i); 
        % nan out the removed points
        t09_co2_2017(i) = NaN; 
    end 
end 

%%%%%%%%%%%% filter the 2020 co2 data by wind direction %%%%%%%%%%%%%%%%%%%

%%%% tower 01 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t01_co2_windir_removed_2020 = NaN(length(t01_co2_2020),1); 

% make a new vector to store the unfiltered co2 for t01 (for reference)
t01_co2_no_winddirfilt_2020 = t01_co2_2020; 

% loop through the tower 1 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 20-65 deg)
for i=1:length(t01_co2_2020)
    % if the wind direction is between 20 and 65 deg
    if wind_data_2020(i,7) >=20 & wind_data_2020(i,7) <= 65 
        % save out the filtered point to the new vector
        t01_co2_windir_removed_2020(i) =t01_co2_2020(i);
        % nan out the filtered point 
        t01_co2_2020(i) = NaN; 
    end 
end 


% loop through tower 1 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t01_co2_2020)
    if isnan(wind_data_2020(i,7))==true
        % save out the removed point to the new vector
        t01_co2_windir_removed_2020(i) = t01_co2_2020(i); 
        % nan out the removed points 
        t01_co2_2020(i) = NaN; 
    end 
end 

%%%% tower 09 filtering %%%%  

% make an empty vector to store the points that get removed from the filter
t09_co2_windir_removed_2020 = NaN(length(t09_co2_2020),1); 

% make a new vector to store the unfiltered co2 for t09 (for reference)
t09_co2_no_winddirfilt_2020 = t09_co2_2020; 

% loop through the tower 9 data and remove points where the tower is in the
% urban plume (according to Miles et al., 2021, that is 235-280 deg)
for i=1:length(t09_co2_2020)
    % if the wind direction is between 235 and 280 deg
    if wind_data_2020(i,7) >=235 & wind_data_2020(i,7) <= 280 
        % save out the filtered point to the new vector
        t09_co2_windir_removed_2020(i) = t09_co2_2020(i); 
        % nan out the point in the data vector 
        t09_co2_2020(i) = NaN; 
    end 
end 

% loop through tower 9 data and remove the points where there is no wind
% direction data (should not be many points)
for i=1:length(t09_co2_2020)
    if isnan(wind_data_2020(i,7))==true
        % save out the removed point to the new vector
        t09_co2_windir_removed_2020(i) =t09_co2_2020(i); 
        % nan out the removed points
        t09_co2_2020(i) = NaN; 
    end 
end 

%% remove points from towers when others are missing data 
% remove points from all towers when any tower is missing data 
% so if one tower is missing data, or if the point was removed through the
% wind direction filtering, we remove that point from all of the towers

%%% for 2017 %%%
% loop through all of the 2017 data, if any of the towers is missing a data
% point, nan out all of the other towers 
for i=1:length(t14_co2_2017) 
    % if any of the towers are missing data 
    if isnan(t14_co2_2017(i))==true | isnan(t09_co2_2017(i))==true | isnan(t01_co2_2017(i))==true  
        % nan out the points in all towers 
        t14_co2_2017(i) = NaN; 
        t01_co2_2017(i) = NaN; 
        t09_co2_2017(i) = NaN; 
    end 
end 

%%% for 2018 %%%
% loop through all of the 2018 data, if any of the towers is missing a data
% point, nan out all of the other towers 
for i=1:length(t14_co2_2018)
    % if any of the towers is missing data 
    if isnan(t14_co2_2018(i))==true | isnan(t09_co2_2018(i))==true | isnan(t01_co2_2018(i))==true  
        % nan out the points in all towers 
        t14_co2_2018(i) = NaN; 
        t01_co2_2018(i) = NaN; 
        t09_co2_2018(i) = NaN; 
    end 
end 


%%% for 2019 %%%
% loop through all of the 2018 data, if any of the towers is missing a data
% point, nan out all of the other towers 
for i=1:length(t14_co2_2019)
    % if any of the towers is missing data 
    if isnan(t14_co2_2019(i))==true | isnan(t09_co2_2019(i))==true | isnan(t01_co2_2019(i))==true  
        % nan out the points in all towers 
        t14_co2_2019(i) = NaN; 
        t01_co2_2019(i) = NaN; 
        t09_co2_2019(i) = NaN; 
    end 
end 


%%% for 2020 %%%
% loop through all of the 2020 data, if any of the towers is missing a data
% point, nan out all of the other towers 
for i=1:length(t14_co2_2020)
    % if any of the towers is missing data 
    if isnan(t14_co2_2020(i))==true | isnan(t09_co2_2020(i))==true | isnan(t01_co2_2020(i))==true  
        % nan out the points in all towers 
        t14_co2_2020(i) = NaN; 
        t01_co2_2020(i) = NaN; 
        t09_co2_2020(i) = NaN; 
    end 
end 


%% read in modeled biogenic co2 mole fracs for towers 01, 09, 14
% read in and format model outputs (estimated impact of biology on mole frac

% path on roar collab 
model_co2_path = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs/'; 

% get the file names for tower 1 2017 
t01_2017_model_files = dir(strcat(model_co2_path,'*Tower_01*2017*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_01*2017*'))

% get the file names for tower 9 2017 
t09_2017_model_files = dir(strcat(model_co2_path,'*Tower_09*2017*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_09*2017*'))

% get the file names for tower 14 2017 
t14_2017_model_files = dir(strcat(model_co2_path,'*Tower_14*2017*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_14*2017*'))


% get the file names for tower 1 2018 
t01_2018_model_files = dir(strcat(model_co2_path,'*Tower_01*2018*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_01*2018*'))

% get the file names for tower 9 2018 
t09_2018_model_files = dir(strcat(model_co2_path,'*Tower_09*2018*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_09*2018*'))

% get the file names for tower 14 2018 
t14_2018_model_files = dir(strcat(model_co2_path,'*Tower_14*2018*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_14*2018*'))


% get the file names for tower 1 2019 
t01_2019_model_files = dir(strcat(model_co2_path,'*Tower_01*2019*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_01*2019*'))

% get the file names for tower 9 2019 
t09_2019_model_files = dir(strcat(model_co2_path,'*Tower_09*2019*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_09*2019*'))

% get the file names for tower 14 2019 
t14_2019_model_files = dir(strcat(model_co2_path,'*Tower_14*2019*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_14*2019*'))


% get the file names for tower 1 2020 
t01_2020_model_files = dir(strcat(model_co2_path,'*Tower_01*2020*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_01*2020*'))

% get the file names for tower 9 2020
t09_2020_model_files = dir(strcat(model_co2_path,'*Tower_09*2020*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_09*2020*'))

% get the file names for tower 14 2020
t14_2020_model_files = dir(strcat(model_co2_path,'*Tower_14*2020*')); 
% check that all the files are there 
dir(strcat(model_co2_path,'*Tower_14*2020*'))


% make an empty vector to store the modeled tower outputs for 2017  
t01_2017_model = []; 
t09_2017_model = []; 
t14_2017_model = []; 

% loop over the files to read in and concatinate the 2017 tower 1 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t01_2017_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t01_2017_model = vertcat(t01_2017_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 


% loop over the files to read in and concatinate the 2017 tower 9 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t09_2017_model_files(i).name)); 

    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t09_2017_model = vertcat(t09_2017_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2017 tower 14 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t14_2017_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t14_2017_model = vertcat(t14_2017_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% make an empty vector to store the modeled tower outputs for 2018  
t01_2018_model = []; 
t09_2018_model = []; 
t14_2018_model = []; 

% loop over the files to read in and concatinate the 2018 tower 1 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t01_2018_model_files(i).name)); 
    % add the current file to the matrix (vertically)
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    t01_2018_model = vertcat(t01_2018_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 


% loop over the files to read in and concatinate the 2018 tower 9 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t09_2018_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t09_2018_model = vertcat(t09_2018_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2018 tower 14 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t14_2018_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t14_2018_model = vertcat(t14_2018_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% make an empty vector to store the modeled tower outputs for 2019  
t01_2019_model = []; 
t09_2019_model = []; 
t14_2019_model = []; 

% loop over the files to read in and concatinate the 2019 tower 1 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t01_2019_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t01_2019_model = vertcat(t01_2019_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2019 tower 9 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t09_2019_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t09_2019_model = vertcat(t09_2019_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2019 tower 14 files 
for i=1:12 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t14_2019_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t14_2019_model = vertcat(t14_2019_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 


% make an empty vector to store the modeled tower outputs for 2020  
t01_2020_model = []; 
t09_2020_model = []; 
t14_2020_model = []; 

% loop over the files to read in and concatinate the 2020 tower 1 files 
for i=1:4 % for all the months in the year - only 4 for 2020 (only year start run)
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t01_2020_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t01_2020_model = vertcat(t01_2020_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2020 tower 9 files 
for i=1:4 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t09_2020_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t09_2020_model = vertcat(t09_2020_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 

% loop over the files to read in and concatinate the 2020 tower 14 files 
for i=1:4 % for all the months in the year 
    % read in the file i of the loop 
    file_in = readmatrix(strcat(model_co2_path,t14_2020_model_files(i).name)); 
    
    % kill the first row of the data (these are nans from reading in the
    % column headers wrong)
    file_in(1,:) = []; 
    
    % add the current file to the matrix (vertically)
    t14_2020_model = vertcat(t14_2020_model,file_in); 
    % clean up to use this in the next loop iteration  
    clear file_in 
end 


%% combine modeled years and combine obs years for each tower %% 

% string together all years of the obs mole fractions 
t01_obs_co2 = vertcat(t01_co2_2017,t01_co2_2018,t01_co2_2019,t01_co2_2020); 
t09_obs_co2 = vertcat(t09_co2_2017,t09_co2_2018,t09_co2_2019,t09_co2_2020); 
t14_obs_co2 = vertcat(t14_co2_2017,t14_co2_2018,t14_co2_2019,t14_co2_2020); 

% make a datetime vector to represent the obs 
all_time_t14_obs = vertcat(t14_dates_2017,t14_dates_2018,t14_dates_2019,t14_dates_2020); 
all_time_t09_obs = vertcat(t09_dates_2017,t09_dates_2018,t09_dates_2019,t09_dates_2020); 
all_time_t01_obs = vertcat(t01_dates_2017,t01_dates_2018,t01_dates_2019,t01_dates_2020); 

% string together all years of the modeled mole fractions - column 2 is the
% mole fraction 
t01_model_co2 = vertcat(t01_2017_model(:,2),t01_2018_model(:,2),t01_2019_model(:,2),t01_2020_model(:,2));  
t09_model_co2 = vertcat(t09_2017_model(:,2),t09_2018_model(:,2),t09_2019_model(:,2),t09_2020_model(:,2));  
t14_model_co2 = vertcat(t14_2017_model(:,2),t14_2018_model(:,2),t14_2019_model(:,2),t14_2020_model(:,2));  

% make datetimes for the times for all the years of model - NOTE: because
% of the influence functions, the data for each year starts on Jan 1, 1
% UTC, and the last hour is Jan 1, 0 UTC of the following year (except
% 2020 for that year i didnt run the full year so it ends May 1, 0 UTC)
model_year_2017(:,1) = datetime(2017,1,1,1,0,0):hours(1):datetime(2018,1,1,0,0,0);
model_year_2018(:,1) = datetime(2018,1,1,1,0,0):hours(1):datetime(2019,1,1,0,0,0);
model_year_2019(:,1) = datetime(2019,1,1,1,0,0):hours(1):datetime(2020,1,1,0,0,0); 
model_year_2020(:,1) = datetime(2020,1,1,1,0,0):hours(1):datetime(2020,5,1,0,0,0); 

% string together the dates of the model 
all_time_model = vertcat(model_year_2017,model_year_2018,model_year_2019,model_year_2020); 

% make just one time series of dates for the obs 
all_time_obs = all_time_t14_obs; 

% now shift the obs time stamps by ONE HALF HOUR (00:30:00 -> 00:00:00) to
% make comparison to model easier 
all_time_obs = all_time_obs - minutes(30); 

% delete all times after May 1, 0UTC 2020 from the obs to match the model
% outputs - wait dont include 0 utc 
model_2020_end = find(all_time_obs==all_time_model(end)); 
all_time_obs((model_2020_end):end) = []; 
t01_obs_co2((model_2020_end):end) = []; 
t09_obs_co2((model_2020_end):end) = []; 
t14_obs_co2((model_2020_end):end) = []; 

% force the obs and data to be the same time frames - first remove the hour
% of 2017 Jan 1 00:00:00 UTC hour in the obs (this is not in the model
% output, plus the data points at that time get filtered out anyway because
% the site 14 data is not available at this time) and I'll cut off the May
% 1 00:00:00 UTC 2020 point in the model 
all_time_obs(1) = []; 
t01_obs_co2(1) = []; 
t09_obs_co2(1) = []; 
t14_obs_co2(1) = []; 

all_time_model(end) = []; 
t01_model_co2(end) = []; 
t09_model_co2(end) = []; 
t14_model_co2(end) = []; 

% check the lengths & start/end times to make sure all vectors are equal 
all_time_obs(1)
all_time_obs(end)

all_time_model(1)
all_time_model(end)

length(all_time_obs)
length(t01_obs_co2)
length(t09_obs_co2)
length(t14_obs_co2)

length(all_time_model)
length(t01_model_co2)
length(t09_model_co2)
length(t14_model_co2)

% remove date variables we don't need anymore 
clear t14_dates_2017 t14_dates_2018 t14_dates_2019 t14_dates_2020
clear t09_dates_2017 t09_dates_2018 t09_dates_2019 t09_dates_2020
clear t01_dates_2017 t01_dates_2018 t01_dates_2019 t01_dates_2020

%% crop time series to remove times with impacts of COVID lockdowns %% 
% remove the times starting at March 10 (march 12 is when all marion
% county schools shut down and march 23 is when lockdown started in
% indiana)

covid_start_date = find(all_time_obs==datetime(2020,03,10,0,0,0)); 
all_time_obs(covid_start_date)
all_time_model(covid_start_date)

all_time_obs((covid_start_date):end) = []; 
t01_obs_co2((covid_start_date):end) = []; 
t09_obs_co2((covid_start_date):end) = []; 
t14_obs_co2((covid_start_date):end) = []; 

all_time_model((covid_start_date):end) = [];
t01_model_co2((covid_start_date):end) = []; 
t09_model_co2((covid_start_date):end) = []; 
t14_model_co2((covid_start_date):end) = []; 

% check the lengths & start/end times to make sure all vectors are equal
all_time_obs(1)
all_time_obs(end)

all_time_model(1)
all_time_model(end)

length(all_time_obs)
length(t01_obs_co2)
length(t09_obs_co2)
length(t14_obs_co2)

length(all_time_model)
length(t01_model_co2)
length(t09_model_co2)
length(t14_model_co2)

%% remove the leap year day (delete feb 29) %% 
% remove from the model and the observations 

% first check that the leap day hours select correctly 
all_time_obs(day(all_time_obs)==29 & month(all_time_obs)==2)
all_time_model(day(all_time_model)==29 & month(all_time_model)==2)

% now remove the leap day from the obs and model (do this before removing
% from the date vectors b/c we use the date vectors to find the leap day
% indeices)
t01_obs_co2(day(all_time_obs)==29 & month(all_time_obs)==2) = []; 
t09_obs_co2(day(all_time_obs)==29 & month(all_time_obs)==2) = [];
t14_obs_co2(day(all_time_obs)==29 & month(all_time_obs)==2) = [];

t01_model_co2(day(all_time_model)==29 & month(all_time_model)==2) = []; 
t09_model_co2(day(all_time_model)==29 & month(all_time_model)==2) = []; 
t14_model_co2(day(all_time_model)==29 & month(all_time_model)==2) = []; 

% now we can remove the dates from the time series 
all_time_obs(day(all_time_obs)==29 & month(all_time_obs)==2) = []; 
all_time_model(day(all_time_model)==29 & month(all_time_model)==2) = []; 

% check once again that everything is the right length 
all_time_obs(1)
all_time_obs(end)

all_time_model(1)
all_time_model(end)

length(all_time_obs)
length(t01_obs_co2)
length(t09_obs_co2)
length(t14_obs_co2)

length(all_time_model)
length(t01_model_co2)
length(t09_model_co2)
length(t14_model_co2)


%% remove model points when missing observations %% 

% check that the obs are filtered correctly (at any time there is a nan,
% all the towers should be nans at that point)
any(isnan(t01_obs_co2) == true & isnan(t09_obs_co2) == false)
any(isnan(t01_obs_co2) == false & isnan(t09_obs_co2) == true)
any(isnan(t09_obs_co2) == false & isnan(t14_obs_co2) == true)  
any(isnan(t09_obs_co2) == true & isnan(t14_obs_co2) == false)  
any(isnan(t01_obs_co2) == false & isnan(t14_obs_co2) == true)  
any(isnan(t01_obs_co2) == true & isnan(t14_obs_co2) == false)

% remove the model data points for times when the data is missing 
for i=1:length(all_time_model)
    if isnan(t01_obs_co2(i)) == true % only using one tower b/c we already filtered all towers to be the same 
        t01_model_co2(i) = NaN; 
        t09_model_co2(i) = NaN;
        t14_model_co2(i) = NaN;
    end     
end 


% quick plot of the model data for all time
figure()
tiledlayout(3,1,'TileSpacing','compact')
nexttile
plot(all_time_model,t01_model_co2,'.','MarkerEdgeColor','black')
ylabel('Modeled CO2 Conc [ppm]')
title('Tower 1')
yline(0)
ylim([-40 50])
grid on 
nexttile
plot(all_time_model,t09_model_co2,'.','MarkerEdgeColor','black')
ylabel('Modeled CO2 Conc [ppm]')
title('Tower 9')
grid on 
yline(0)
ylim([-40 50])
nexttile
plot(all_time_model,t14_model_co2,'.','MarkerEdgeColor','black')
grid on 
xlabel('Time (UTC)')
ylabel('Modeled CO2 Conc [ppm]')
title('Tower 14')
yline(0)
ylim([-40 50])


%% subset afternoon hours for observations and model %% 

% subset obs and model hours seperate (however they should be the same
% indicies)
ind_afternoon_hours_obs = find(hour(all_time_obs) >= 17 & hour(all_time_obs) <= 22); 
ind_afternoon_hours_model = find(hour(all_time_model) >= 17 & hour(all_time_model) <= 22); 

% grab only the afternoon hours for each modeled tower bio conc
t01_model_afternoon_co2(:,1) = t01_model_co2(ind_afternoon_hours_model); 
t09_model_afternoon_co2(:,1) = t09_model_co2(ind_afternoon_hours_model);
t14_model_afternoon_co2(:,1) = t14_model_co2(ind_afternoon_hours_model);

% grab only the afternoon hours for each obs tower conc
t01_obs_afternoon_co2(:,1) = t01_obs_co2(ind_afternoon_hours_obs);
t09_obs_afternoon_co2(:,1) = t09_obs_co2(ind_afternoon_hours_obs);
t14_obs_afternoon_co2(:,1) = t14_obs_co2(ind_afternoon_hours_obs);

% take a look at just the afternoon hours 
afternoon_hrs = all_time_model(ind_afternoon_hours_model); 

% quick plot of the model afternoon hours only 
figure()
t= tiledlayout(3,1,"TileSpacing","compact"); 
nexttile
plot(all_time_model(ind_afternoon_hours_model),t01_model_afternoon_co2,'.')
title('tower 01')
grid on 
nexttile
plot(all_time_model(ind_afternoon_hours_model),t09_model_afternoon_co2,'.')
title('tower 09')
grid on 
nexttile
plot(all_time_model(ind_afternoon_hours_model),t14_model_afternoon_co2,'.')
title('tower 14')
grid on 
ylabel(t,'Model CO2 (ppm)')
title(t,'Modeled CO2 - Afternoon Hours')


% quick plot of the obs afternoon hours only 
figure()
t = tiledlayout(3,1,"TileSpacing","compact"); 
nexttile
plot(all_time_obs(ind_afternoon_hours_obs),t01_obs_afternoon_co2,'.')
title('tower 01')
grid on 
nexttile
plot(all_time_obs(ind_afternoon_hours_obs),t09_obs_afternoon_co2,'.')
title('tower 09')
grid on 
nexttile
plot(all_time_obs(ind_afternoon_hours_obs),t14_obs_afternoon_co2,'.')
title('tower 14')
grid on 
ylabel(t,'Observed CO2 (ppm)')
title(t,'Observed CO2 - Afternoon Hours')

%% calculate afternoon average mole fraction enhancements for model and obs %% 
% NOTE: enhancement = differences between background towers
% NOTE: calculate the difference between afternoon hours at towers
% and then take the average of that value 

% check that the reshaping works as i want it to - yay look good! 
afternoon_hrs_reshape_test = reshape(afternoon_hrs,6,[]);

% calculate the enhancements (tower differences) for the obs hours in afternoon
t01_t09_obs_afternoon_co2_diff(:,1) = t01_obs_afternoon_co2 - t09_obs_afternoon_co2; 
t14_t09_obs_afternoon_co2_diff(:,1) = t14_obs_afternoon_co2 - t09_obs_afternoon_co2; 

% calculate the mean difference between towers for each afternoon - OBS
t01_t09_obs_afternoon_co2_enhc_mean(:,1) =  mean(reshape(t01_t09_obs_afternoon_co2_diff,6,[]),'omitnan'); 
t14_t09_obs_afternoon_co2_enhc_mean(:,1) =  mean(reshape(t14_t09_obs_afternoon_co2_diff,6,[]),'omitnan'); 

% calculate the median difference between towers for each afternoon - OBS
t01_t09_obs_afternoon_co2_enhc_med(:,1) =  median(reshape(t01_t09_obs_afternoon_co2_diff,6,[]),'omitnan'); 
t14_t09_obs_afternoon_co2_enhc_med(:,1) =  median(reshape(t14_t09_obs_afternoon_co2_diff,6,[]),'omitnan');  

% calculate the enhancements (tower differences) for the model in afternoon
t01_t09_model_afternoon_co2_diff(:,1) = t01_model_afternoon_co2 - t09_model_afternoon_co2; 
t14_t09_model_afternoon_co2_diff(:,1) = t14_model_afternoon_co2 - t09_model_afternoon_co2; 

% calculate the mean difference between towers for each afternoon - model
t01_t09_model_afternoon_co2_enhc_mean(:,1) =  mean(reshape(t01_t09_model_afternoon_co2_diff,6,[]),'omitnan'); 
t14_t09_model_afternoon_co2_enhc_mean(:,1) =  mean(reshape(t14_t09_model_afternoon_co2_diff,6,[]),'omitnan'); 

% calculate the median difference between towers for each afternoon - model
t01_t09_model_afternoon_co2_enhc_med(:,1) =  median(reshape(t01_t09_model_afternoon_co2_diff,6,[]),'omitnan'); 
t14_t09_model_afternoon_co2_enhc_med(:,1) =  median(reshape(t14_t09_model_afternoon_co2_diff,6,[]),'omitnan'); 

% make a datetime vector to represent every day (to plot the afternoon
% averages) 
all_time_days(:,1) = datetime(year(afternoon_hrs(1)),month(afternoon_hrs(1)),day(afternoon_hrs(1))):days(1):datetime(year(afternoon_hrs(end)),month(afternoon_hrs(end)),day(afternoon_hrs(end)));

% DELETE THE LEAP DAY FROM THIS DAILY DATETIME VECTOR 
all_time_days(find(month(all_time_days)==2 & day(all_time_days)==29)) =[]; 

% check that the obs and model are the same lenth as the datetime vector 
length(all_time_days)
length(t01_t09_obs_afternoon_co2_enhc_mean)
length(t14_t09_obs_afternoon_co2_enhc_mean)
length(t01_t09_model_afternoon_co2_enhc_mean)
length(t14_t09_model_afternoon_co2_enhc_mean)

% plot the modeled and observed differences between towers - afternoon
% average 
figure()
tiledlayout(2,1,"TileSpacing","compact")
nexttile
plot(all_time_days,t01_t09_obs_afternoon_co2_enhc_mean,'black','LineWidth',1)
hold on 
plot(all_time_days,t01_t09_model_afternoon_co2_enhc_mean,'red','LineWidth',1)
yline(0)
legend('Observation','Model')
ylabel('\Delta CO_2 [ppm]')
title('Mean Afternoon Tower 1 and Tower 9 Difference')
set(gca,'Xticklabel',[])
%set(gca,'FontSize',20)
nexttile
plot(all_time_days,t14_t09_obs_afternoon_co2_enhc_mean,'black','LineWidth',1)
hold on 
plot(all_time_days,t14_t09_model_afternoon_co2_enhc_mean,'red','LineWidth',1)
yline(0)
legend('Observation','Model')
xlabel('Time (UTC)')
ylabel('\Delta CO_2 [ppm]')
title('Mean Afternoon Tower 14 and Tower 9 Difference')
%set(gca,'FontSize',20)


% plot the modeled and observed differences between towers - afternoon
% median 
figure()
tiledlayout(2,1,"TileSpacing","compact")
nexttile
plot(all_time_days,t01_t09_obs_afternoon_co2_enhc_med,'black','LineWidth',1)
hold on 
plot(all_time_days,t01_t09_model_afternoon_co2_enhc_med,'red','LineWidth',1)
yline(0)
legend('Observation','Model')
ylabel('\Delta CO_2 [ppm]')
title('Median Afternoon Tower 1 and Tower 9 Difference')
set(gca,'Xticklabel',[])
%set(gca,'FontSize',20)
nexttile
plot(all_time_days,t14_t09_obs_afternoon_co2_enhc_med,'black','LineWidth',1)
hold on 
plot(all_time_days,t14_t09_model_afternoon_co2_enhc_med,'red','LineWidth',1)
yline(0)
legend('Observation','Model')
xlabel('Time (UTC)')
ylabel('\Delta CO_2 [ppm]')
title('Median Afternoon Tower 14 and Tower 9 Difference')
%set(gca,'FontSize',20)

%% plot observed afternoon average enhancements %% 

figure()
tiledlayout(2,1,"TileSpacing","compact")
nexttile
plot(all_time_days,t01_t09_obs_afternoon_co2_enhc_med,'.','MarkerEdgeColor','black')
grid on 
yline(0)
ylabel('\Delta CO_2 [ppm]')
set(gca,'Xticklabel',[])
%set(gca,'FontSize',20)
nexttile
plot(all_time_days,t14_t09_obs_afternoon_co2_enhc_med,'.','MarkerEdgeColor','black')
grid on 
yline(0)
xlabel('Time (UTC)')
ylabel('\Delta CO_2 [ppm]')
%set(gca,'FontSize',20)


%% calculate the model observation residuals 
% calculate residuals, model - obs, using the afternoon average differences
% (enhancemnts)

% calculate the residuals for each mean afternoon value 
resid_t01_t09_mean_afternoon_enhc = t01_t09_model_afternoon_co2_enhc_mean - t01_t09_obs_afternoon_co2_enhc_mean; 
resid_t14_t09_mean_afternoon_enhc = t14_t09_model_afternoon_co2_enhc_mean - t14_t09_obs_afternoon_co2_enhc_mean; 

%% calculate mean enhancements and residuals for each year %%  
% note that this is for each year, so 2017 and 2020 are shorter than 2018 
% and 2019, the calculation is done by creating vectors of ~365 days later 
% note that this this didn't end up in the paper - just looking 

length(all_time_days)
length(resid_t01_t09_mean_afternoon_enhc)
length(resid_t14_t09_mean_afternoon_enhc)

% calculate the yearly mean residuals and standard error 
t01_mean_resid_2017 = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)),'omitnan'); 
t14_mean_resid_2017 = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)),'omitnan'); 

t01_mean_resid_2018 = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018)),'omitnan');
t14_mean_resid_2018 = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018)),'omitnan');

t01_mean_resid_2019 = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019)),'omitnan');
t14_mean_resid_2019 = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019)),'omitnan');

t01_mean_resid_2020 = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020)),'omitnan');
t14_mean_resid_2020 = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020)),'omitnan');

t01_stderr_resid_2017 = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017))))); 
t14_stderr_resid_2017 = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017))))); 

t01_stderr_resid_2018 = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018)),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018))))); 
t14_stderr_resid_2018 = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018)),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018))))); 

t01_stderr_resid_2019 = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019)),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019))))); 
t14_stderr_resid_2019 = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019)),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019))))); 

t01_stderr_resid_2020 = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020)),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020))))); 
t14_stderr_resid_2020 = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020)),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020))))); 

% make a datetime vector of years 
years_datetime = (datetime(2017,1,1,0,0,0):calyears(1):datetime(2020,1,1,0,0,0))'; 

% quick plot 
figure()
tiledlayout(1,2,"TileSpacing","compact")
nexttile
errorbar(1:4,[t01_mean_resid_2017,t01_mean_resid_2018,t01_mean_resid_2019,t01_mean_resid_2020],[2*t01_stderr_resid_2017,2*t01_stderr_resid_2018,2*t01_stderr_resid_2019,2*t01_stderr_resid_2020],'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
grid on 
yline(0)
xlim([0.5 4.5])
xticks([1 2 3 4])
xticklabels({'2017','2018','2019','2020'})
ylabel('Mean Yearly Residual [ppm]')
%set(gca,'Xticklabel',[])
%set(gca,'FontSize',20)
nexttile
errorbar(1:4,[t14_mean_resid_2017,t14_mean_resid_2018,t14_mean_resid_2019,t14_mean_resid_2020],[2*t14_stderr_resid_2017,2*t14_stderr_resid_2018,2*t14_stderr_resid_2019,2*t14_stderr_resid_2020],'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
grid on 
yline(0)
xlim([0.5 4.5])
xticks([1 2 3 4])
xticklabels({'2017','2018','2019','2020'})
%set(gca,'FontSize',20)


% now calculate the mean yearly enhancents for obs and model 
t01_mean_model_enhc_2017 = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017)),'omitnan'); 
t14_mean_model_enhc_2017 = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017)),'omitnan'); 

t01_mean_obs_enhc_2017 = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017)),'omitnan'); 
t14_mean_obs_enhc_2017 = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017)),'omitnan'); 

t01_mean_model_enhc_2018 = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018)),'omitnan'); 
t14_mean_model_enhc_2018 = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018)),'omitnan'); 

t01_mean_obs_enhc_2018 = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018)),'omitnan'); 
t14_mean_obs_enhc_2018 = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018)),'omitnan'); 

t01_mean_model_enhc_2019 = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019)),'omitnan'); 
t14_mean_model_enhc_2019 = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019)),'omitnan'); 

t01_mean_obs_enhc_2019 = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019)),'omitnan'); 
t14_mean_obs_enhc_2019 = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019)),'omitnan'); 

t01_mean_model_enhc_2020 = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2020)),'omitnan'); 
t14_mean_model_enhc_2020 = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2020)),'omitnan'); 

t01_mean_obs_enhc_2020 = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2020)),'omitnan'); 
t14_mean_obs_enhc_2020 = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2020)),'omitnan'); 

% string together yearly means 
t01_yearly_mean_enhc_model = [t01_mean_model_enhc_2017,t01_mean_model_enhc_2018,t01_mean_model_enhc_2019,t01_mean_model_enhc_2020]; 
t14_yearly_mean_enhc_model = [t14_mean_model_enhc_2017,t14_mean_model_enhc_2018,t14_mean_model_enhc_2019,t14_mean_model_enhc_2020]; 

t01_yearly_mean_enhc_obs = [t01_mean_obs_enhc_2017,t01_mean_obs_enhc_2018,t01_mean_obs_enhc_2019,t01_mean_obs_enhc_2020]; 
t14_yearly_mean_enhc_obs = [t14_mean_obs_enhc_2017,t14_mean_obs_enhc_2018,t14_mean_obs_enhc_2019,t14_mean_obs_enhc_2020]; 

% plot yearly means and residuals 
figure()
tiledlayout(2,2,"TileSpacing","compact")
nexttile
plot(1:4,t01_yearly_mean_enhc_obs,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
hold on 
plot(1:4,t01_yearly_mean_enhc_model,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
title('Tower 01')
yline(0)
xlim([0.5 4.5])
text(0.7,-0.1,'(a)')
xticks([1 2 3 4])
grid on 
ylim([-1 0])
xticklabels([])
ylabel({'Mean Yearly', 'CO_2 Enhancements [ppm]'})
legend({'observations','model'})

nexttile
plot(1:4,t14_yearly_mean_enhc_obs,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
hold on 
plot(1:4,t14_yearly_mean_enhc_model,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
title('Tower 14')
yline(0)
text(0.7,-0.1,'(b)')
xlim([0.5 4.5])
xticks([1 2 3 4])
xticklabels([])
grid on 
ylim([-1 0])

nexttile
errorbar(1:4,[t01_mean_resid_2017,t01_mean_resid_2018,t01_mean_resid_2019,t01_mean_resid_2020],[2*t01_stderr_resid_2017,2*t01_stderr_resid_2018,2*t01_stderr_resid_2019,2*t01_stderr_resid_2020],'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
grid on 
yline(0)
xlim([0.5 4.5])
xticks([1 2 3 4])
ylim([-1 1.25])
text(0.7,1,'(c)')
xticklabels({'2017','2018','2019','2020'})
ylabel('Mean Yearly Residuals [ppm]')
%set(gca,'Xticklabel',[])
%set(gca,'FontSize',20)
nexttile
errorbar(1:4,[t14_mean_resid_2017,t14_mean_resid_2018,t14_mean_resid_2019,t14_mean_resid_2020],[2*t14_stderr_resid_2017,2*t14_stderr_resid_2018,2*t14_stderr_resid_2019,2*t14_stderr_resid_2020],'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
grid on 
yline(0)
xlim([0.5 4.5])
ylim([-1 1.25])
text(0.7,1,'(d)')
xticks([1 2 3 4])
xticklabels({'2017','2018','2019','2020'})
%set(gca,'FontSize',20)

%% calculate the monthly mean residuals %% 
% calculate the monthly mean residuals using the afternoon averages 

% get the years of data available 
all_time_years = unique(year(all_time_days)); 

% make NaN vectors for the number of months in the data set to store
% monthly mean residuals, median residuals, stdev of residuals, and std
% error of the residuals 
resid_t01_t09_monthly_mean = NaN(length(all_time_years)*12,1); 
resid_t14_t09_monthly_mean = NaN(length(all_time_years)*12,1); 

resid_t01_t09_monthly_med = NaN(length(all_time_years)*12,1); 
resid_t14_t09_monthly_med = NaN(length(all_time_years)*12,1); 

resid_t01_t09_monthly_std = NaN(length(all_time_years)*12,1); 
resid_t14_t09_monthly_std = NaN(length(all_time_years)*12,1);

resid_t01_t09_monthly_stderr = NaN(length(all_time_years)*12,1); 
resid_t14_t09_monthly_stderr = NaN(length(all_time_years)*12,1); 

% make NaN vectors to store monthly mean enhancemnts for model and obs 
co2_enhc_model_t01_t09_monthly_mean = NaN(length(all_time_years)*12,1);
co2_enhc_model_t14_t09_monthly_mean = NaN(length(all_time_years)*12,1);
co2_enhc_obs_t01_t09_monthly_mean = NaN(length(all_time_years)*12,1);
co2_enhc_obs_t14_t09_monthly_mean = NaN(length(all_time_years)*12,1); 

% make a vector of all the months in the data set (needs to have a day so
% it defaults to the first of the month)
all_time_months(:,1) = datetime(year(all_time_days(1)),month(all_time_days(1)),day(all_time_days(1))):calmonths(1):datetime(year(all_time_days(end)),month(all_time_days(end)),day(all_time_days(end))); 

% now i need to manually append april 2020 to dec 2020 to the vector time
% series (this is just for plotting the 2020 bar chart)
all_time_months = vertcat(all_time_months,(datetime(2020,4,1):calmonths(1):datetime(2020,12,1))'); 

% loop over all the months of data and calculate the mean monthly residual
% along with the mean monthly enhancement for model and obs 
n=1; 
for i=1:length(unique(year(all_time_days))) % for as many unique years of data 
    
    % grab the year 
    year_i = all_time_years(i); 

    for j=1:12 % for all 12 months in each year 

        % get the indices of the year and month (for the daily spaced
        % datetime vector)
        ind_month_j = find(year(all_time_days)==year_i & month(all_time_days)==j); 

        % calculate the average of the residuals for the month of interest 
        month_mean_t01_t09_j = mean(resid_t01_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 
        month_mean_t14_t09_j = mean(resid_t14_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 

        % calculate the average modeled enhancements for the month of interest 
        month_mean_model_enhc_t01_t09_j = mean(t01_t09_model_afternoon_co2_enhc_mean(ind_month_j),'omitnan'); 
        month_mean_model_enhc_t14_t09_j = mean(t14_t09_model_afternoon_co2_enhc_mean(ind_month_j),'omitnan'); 

        % calculate the average modeled enhancements for the month of interest 
        month_mean_obs_enhc_t01_t09_j = mean(t01_t09_obs_afternoon_co2_enhc_mean(ind_month_j),'omitnan'); 
        month_mean_obs_enhc_t14_t09_j = mean(t14_t09_obs_afternoon_co2_enhc_mean(ind_month_j),'omitnan'); 

        % calculate the median of the residuals for the month of interest 
        month_med_t01_t09_j = median(resid_t01_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 
        month_med_t14_t09_j = median(resid_t14_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 

        % calculate the std of the residuals for the month of interest 
        month_std_t01_t09_j = std(resid_t01_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 
        month_std_t14_t09_j = std(resid_t14_t09_mean_afternoon_enhc(ind_month_j),'omitnan'); 
        
        % calculate the standard error of the residuals for the month of interest 
        month_stderr_t01_t09_j = month_std_t01_t09_j/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(ind_month_j)))); 
        month_stderr_t14_t09_j = month_std_t14_t09_j/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(ind_month_j)))); 
       
        % save out the monthly calculated mean residual 
        resid_t01_t09_monthly_mean(n) = month_mean_t01_t09_j; 
        resid_t14_t09_monthly_mean(n) = month_mean_t14_t09_j; 

        % save out the monthly calculated mean enhancement (model and obs)
        co2_enhc_model_t01_t09_monthly_mean(n) = month_mean_model_enhc_t01_t09_j; 
        co2_enhc_model_t14_t09_monthly_mean(n) = month_mean_model_enhc_t14_t09_j; 
        co2_enhc_obs_t01_t09_monthly_mean(n) = month_mean_obs_enhc_t01_t09_j; 
        co2_enhc_obs_t14_t09_monthly_mean(n) = month_mean_obs_enhc_t14_t09_j; 

        % save out the monthly calculated median residual 
        resid_t01_t09_monthly_med(n) = month_med_t01_t09_j; 
        resid_t14_t09_monthly_med(n) = month_med_t14_t09_j; 

        % save out the monthly calculated standard deviation of the
        % residuals
        resid_t01_t09_monthly_std(n) = month_std_t01_t09_j; 
        resid_t14_t09_monthly_std(n) = month_std_t14_t09_j; 
        
        % save out the monthly calculated standard errors of the
        % residuals
        resid_t01_t09_monthly_stderr(n) = month_stderr_t01_t09_j; 
        resid_t14_t09_monthly_stderr(n) = month_stderr_t14_t09_j; 


        n=n+1; 
        clear ind_month_j month_mean_t01_t09_j month_mean_t14_t09_j 
        clear month_med_t01_t09_j month_med_t14_t09_j 
        clear month_std_t01_t09_j month_std_t14_t09_j 
        clear month_stderr_t01_t09_j month_stderr_t14_t09_j 
        clear month_mean_model_enhc_t01_t09_j month_mean_model_enhc_t14_t09_j
        clear month_mean_obs_enhc_t01_t09_j month_mean_obs_enhc_t14_t09_j

    end     
    clear year_i 
end 

% calculate & look at the yearly mean of monthly residuals 
mean(resid_t01_t09_monthly_mean(find(year(all_time_months)==2017)),'omitnan')
mean(resid_t14_t09_monthly_mean(find(year(all_time_months)==2017)),'omitnan')

mean(resid_t01_t09_monthly_mean(find(year(all_time_months)==2018)),'omitnan')
mean(resid_t14_t09_monthly_mean(find(year(all_time_months)==2018)),'omitnan')

mean(resid_t01_t09_monthly_mean(find(year(all_time_months)==2019)),'omitnan')
mean(resid_t14_t09_monthly_mean(find(year(all_time_months)==2019)),'omitnan')

mean(resid_t01_t09_monthly_mean(find(year(all_time_months)==2020)),'omitnan')
mean(resid_t14_t09_monthly_mean(find(year(all_time_months)==2020)),'omitnan')

% plot monthly mean residuals with standard error 
figure()
t = tiledlayout(3,2,"TileSpacing","compact"); 
nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2017))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2017)),resid_t01_t09_monthly_stderr(find(year(all_time_months)==2017)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2017))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2017)),resid_t14_t09_monthly_stderr(find(year(all_time_months)==2017)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])


nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2018))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2018)),resid_t01_t09_monthly_stderr(find(year(all_time_months)==2018)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2018))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2018)),resid_t14_t09_monthly_stderr(find(year(all_time_months)==2018)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])

nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2019))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2019)),resid_t01_t09_monthly_stderr(find(year(all_time_months)==2019)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2019))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2019)),resid_t14_t09_monthly_stderr(find(year(all_time_months)==2019)),'k.'); 
hold off 
ylim([-4.4 6])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})

ylabel(t,'Monthly Mean Residual [ppm]')%,'FontSize',24)

% plot monthly mean residuals with two times the standard error 
figure()
t = tiledlayout(4,2,"TileSpacing","compact"); 
nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2017))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2017)),2.*resid_t01_t09_monthly_stderr(find(year(all_time_months)==2017)),'k.'); 
text(0.2,-8,'2017','FontSize',18)
text(11,5,'(a)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
title('Tower 01')

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2017))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2017)),2.*resid_t14_t09_monthly_stderr(find(year(all_time_months)==2017)),'k.'); 
text(0.2,-8,'2017','FontSize',18)
text(11,5,'(b)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
title('Tower 14')
set(gca,'Yticklabel',[])

nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2018))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2018)),2.*resid_t01_t09_monthly_stderr(find(year(all_time_months)==2018)),'k.'); 
text(0.2,-8,'2018','FontSize',18)
text(11,5,'(c)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2018))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2018)),2.*resid_t14_t09_monthly_stderr(find(year(all_time_months)==2018)),'k.'); 
text(0.2,-8,'2018','FontSize',18)
text(11,5,'(d)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2019))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2019)),2.*resid_t01_t09_monthly_stderr(find(year(all_time_months)==2019)),'k.'); 
text(0.2,-8,'2019','FontSize',18)
text(11,5,'(e)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2019))); 
b1.FaceColor = '#96BEE6'; 
hold on 
text(0.2,-8,'2019','FontSize',18)
text(11,5,'(f)','FontSize',18)
grid on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2019)),2.*resid_t14_t09_monthly_stderr(find(year(all_time_months)==2019)),'k.'); 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

nexttile 
b1=bar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2020))); 
b1.FaceColor = '#96BEE6'; 
hold on 
er1 = errorbar(1:12,resid_t01_t09_monthly_mean(find(year(all_time_months)==2020)),2.*resid_t01_t09_monthly_stderr(find(year(all_time_months)==2020)),'k.'); 
text(0.2,-8,'2020','FontSize',18)
text(11,5,'(g)','FontSize',18)
grid on 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
xtickangle(360) % i dont know why i have to add this.. something is up with the labels and they were getting rotated weird

nexttile 
b1=bar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2020))); 
b1.FaceColor = '#96BEE6'; 
hold on 
text(0.2,-8,'2020','FontSize',18)
text(11,5,'(h)','FontSize',18)
grid on 
er1 = errorbar(1:12,resid_t14_t09_monthly_mean(find(year(all_time_months)==2020)),2.*resid_t14_t09_monthly_stderr(find(year(all_time_months)==2020)),'k.'); 
hold off 
ylim([-10 6.5])
set(gca,'FontSize',24)
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
xtickangle(360)
ylabel(t,'Mean Monthly Residuals [ppm]','FontSize',16)


% confirm how many months are within two times the standard error (hard to
% tell in the figure for a few months) - should do this for pub 
%resid_t01_t09_monthly_mean
%resid_t01_t09_monthly_stderr

%resid_t14_t09_monthly_mean
%resid_t14_t09_monthly_stderr

%% calculate enhancements and residuals by season %% 
% seasons are defined as AMJ=April,May,Jun; JAS=July,Aug,Sept; OND =
% Oct,Nov,Dec; JFM=Jan,Feb,Mar

% check that all vectors are the same length 
length(all_time_days)
length(resid_t01_t09_mean_afternoon_enhc)
length(resid_t14_t09_mean_afternoon_enhc)

% test that I am finding dates correctly 
all_time_days(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6) ))

%%%%%% calculate mean enhancements - 2017 seasons %%%%%% 
t01_mean_model_enhc_2017_AMJ = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_model_enhc_2017_AMJ = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_mean_obs_enhc_2017_AMJ = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_obs_enhc_2017_AMJ = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 

t01_mean_model_enhc_2017_JAS = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_model_enhc_2017_JAS = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_mean_obs_enhc_2017_JAS = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_obs_enhc_2017_JAS = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 

t01_mean_model_enhc_2017_OND = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_model_enhc_2017_OND = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_mean_obs_enhc_2017_OND = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_obs_enhc_2017_OND = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 

%%%%%% calculate mean enhancements - 2018 seasons %%%%%% 
t01_mean_model_enhc_2018_JFM = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_model_enhc_2018_JFM = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_mean_obs_enhc_2018_JFM = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_obs_enhc_2018_JFM = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 

t01_mean_model_enhc_2018_AMJ = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_model_enhc_2018_AMJ = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_mean_obs_enhc_2018_AMJ = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_obs_enhc_2018_AMJ = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 

t01_mean_model_enhc_2018_JAS = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_model_enhc_2018_JAS = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_mean_obs_enhc_2018_JAS = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_obs_enhc_2018_JAS = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 

t01_mean_model_enhc_2018_OND = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_model_enhc_2018_OND = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_mean_obs_enhc_2018_OND = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_obs_enhc_2018_OND = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 

%%%%%% calculate mean enhancements - 2019 seasons %%%%%% 
t01_mean_model_enhc_2019_JFM = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_model_enhc_2019_JFM = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_mean_obs_enhc_2019_JFM = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_obs_enhc_2019_JFM = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 

t01_mean_model_enhc_2019_AMJ = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_model_enhc_2019_AMJ = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_mean_obs_enhc_2019_AMJ = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_obs_enhc_2019_AMJ = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 

t01_mean_model_enhc_2019_JAS = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_model_enhc_2019_JAS = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_mean_obs_enhc_2019_JAS = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_obs_enhc_2019_JAS = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 

t01_mean_model_enhc_2019_OND = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_model_enhc_2019_OND = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_mean_obs_enhc_2019_OND = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_obs_enhc_2019_OND = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 

%%%%%% calculate mean enhancements - 2020 season (only Jan, Feb, March) %%%%%% 
t01_mean_model_enhc_2020_JFM = mean(t01_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_model_enhc_2020_JFM = mean(t14_t09_model_afternoon_co2_enhc_mean(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_mean_obs_enhc_2020_JFM = mean(t01_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_obs_enhc_2020_JFM = mean(t14_t09_obs_afternoon_co2_enhc_mean(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 


%%%%%% calculate mean & stderr of residuals - 2017 seasons %%%%%% 
t01_mean_resid_2017_AMJ = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_resid_2017_AMJ = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_stderr_resid_2017_AMJ = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 
t14_stderr_resid_2017_AMJ = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 

t01_mean_resid_2017_JAS = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_resid_2017_JAS = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_stderr_resid_2017_JAS = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 
t14_stderr_resid_2017_JAS = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 

t01_mean_resid_2017_OND = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_resid_2017_OND = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_stderr_resid_2017_OND = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 
t14_stderr_resid_2017_OND = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 

%%%%%% calculate mean & stderr of residuals - 2018 seasons %%%%%% 
t01_mean_resid_2018_JFM = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_resid_2018_JFM = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_stderr_resid_2018_JFM = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 
t14_stderr_resid_2018_JFM = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 

t01_mean_resid_2018_AMJ = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_resid_2018_AMJ = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_stderr_resid_2018_AMJ = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 
t14_stderr_resid_2018_AMJ = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 

t01_mean_resid_2018_JAS = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_resid_2018_JAS = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_stderr_resid_2018_JAS = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 
t14_stderr_resid_2018_JAS = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 

t01_mean_resid_2018_OND = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_resid_2018_OND = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_stderr_resid_2018_OND = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 
t14_stderr_resid_2018_OND = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 

%%%%%% calculate mean & stderr of residuals - 2019 seasons %%%%%% 
t01_mean_resid_2019_JFM = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_resid_2019_JFM = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_stderr_resid_2019_JFM = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 
t14_stderr_resid_2019_JFM = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 

t01_mean_resid_2019_AMJ = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t14_mean_resid_2019_AMJ = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan'); 
t01_stderr_resid_2019_AMJ = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 
t14_stderr_resid_2019_AMJ = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6)))))); 

t01_mean_resid_2019_JAS = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t14_mean_resid_2019_JAS = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan'); 
t01_stderr_resid_2019_JAS = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 
t14_stderr_resid_2019_JAS = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9)))))); 

t01_mean_resid_2019_OND = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t14_mean_resid_2019_OND = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan'); 
t01_stderr_resid_2019_OND = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 
t14_stderr_resid_2019_OND = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12)))))); 

%%%%%% calculate mean & stderr of residuals - 2020 season %%%%%% 
t01_mean_resid_2020_JFM = mean(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t14_mean_resid_2020_JFM = mean(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan'); 
t01_stderr_resid_2020_JFM = std(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 
t14_stderr_resid_2020_JFM = std(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))),'omitnan')/sqrt(sum(~isnan(resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3)))))); 

%%%%%% string together seasons %%%%%% 
t01_seasonal_model_enhc = [t01_mean_model_enhc_2017_AMJ,t01_mean_model_enhc_2017_JAS,t01_mean_model_enhc_2017_OND,t01_mean_model_enhc_2018_JFM,t01_mean_model_enhc_2018_AMJ,t01_mean_model_enhc_2018_JAS,t01_mean_model_enhc_2018_OND...
    t01_mean_model_enhc_2019_JFM,t01_mean_model_enhc_2019_AMJ,t01_mean_model_enhc_2019_JAS,t01_mean_model_enhc_2019_OND,t01_mean_model_enhc_2020_JFM]; 
t14_seasonal_model_enhc = [t01_mean_model_enhc_2017_AMJ,t14_mean_model_enhc_2017_JAS,t14_mean_model_enhc_2017_OND,t14_mean_model_enhc_2018_JFM,t14_mean_model_enhc_2018_AMJ,t14_mean_model_enhc_2018_JAS,t14_mean_model_enhc_2018_OND...
    t14_mean_model_enhc_2019_JFM,t14_mean_model_enhc_2019_AMJ,t14_mean_model_enhc_2019_JAS,t14_mean_model_enhc_2019_OND,t14_mean_model_enhc_2020_JFM]; 

t01_seasonal_obs_enhc = [t01_mean_obs_enhc_2017_AMJ,t01_mean_obs_enhc_2017_JAS,t01_mean_obs_enhc_2017_OND,t01_mean_obs_enhc_2018_JFM,t01_mean_obs_enhc_2018_AMJ,t01_mean_obs_enhc_2018_JAS,t01_mean_obs_enhc_2018_OND...
    t01_mean_obs_enhc_2019_JFM,t01_mean_obs_enhc_2019_AMJ,t01_mean_obs_enhc_2019_JAS,t01_mean_obs_enhc_2019_OND,t01_mean_obs_enhc_2020_JFM]; 
t14_seasonal_obs_enhc = [t01_mean_obs_enhc_2017_AMJ,t14_mean_obs_enhc_2017_JAS,t14_mean_obs_enhc_2017_OND,t14_mean_obs_enhc_2018_JFM,t14_mean_obs_enhc_2018_AMJ,t14_mean_obs_enhc_2018_JAS,t14_mean_obs_enhc_2018_OND...
    t14_mean_obs_enhc_2019_JFM,t14_mean_obs_enhc_2019_AMJ,t14_mean_obs_enhc_2019_JAS,t14_mean_obs_enhc_2019_OND,t14_mean_obs_enhc_2020_JFM]; 

t01_seasonal_mean_resid = [t01_mean_resid_2017_AMJ,t01_mean_resid_2017_JAS,t01_mean_resid_2017_OND,t01_mean_resid_2018_JFM,t01_mean_resid_2018_AMJ,t01_mean_resid_2018_JAS,t01_mean_resid_2018_OND...
    t01_mean_resid_2019_JFM,t01_mean_resid_2019_AMJ,t01_mean_resid_2019_JAS,t01_mean_resid_2019_OND,t01_mean_resid_2020_JFM]; 
t14_seasonal_mean_resid = [t14_mean_resid_2017_AMJ,t14_mean_resid_2017_JAS,t14_mean_resid_2017_OND,t14_mean_resid_2018_JFM,t14_mean_resid_2018_AMJ,t14_mean_resid_2018_JAS,t14_mean_resid_2018_OND...
    t14_mean_resid_2019_JFM,t14_mean_resid_2019_AMJ,t14_mean_resid_2019_JAS,t14_mean_resid_2019_OND,t14_mean_resid_2020_JFM]; 

t01_seasonal_stderr_resid = [t01_stderr_resid_2017_AMJ,t01_stderr_resid_2017_JAS,t01_stderr_resid_2017_OND,t01_stderr_resid_2018_JFM,t01_stderr_resid_2018_AMJ,t01_stderr_resid_2018_JAS,t01_stderr_resid_2018_OND...
    t01_stderr_resid_2019_JFM,t01_stderr_resid_2019_AMJ,t01_stderr_resid_2019_JAS,t01_stderr_resid_2019_OND,t01_stderr_resid_2020_JFM]; 
t14_seasonal_stderr_resid = [t14_stderr_resid_2017_AMJ,t14_stderr_resid_2017_JAS,t14_stderr_resid_2017_OND,t14_stderr_resid_2018_JFM,t14_stderr_resid_2018_AMJ,t14_stderr_resid_2018_JAS,t14_stderr_resid_2018_OND...
    t14_stderr_resid_2019_JFM,t14_stderr_resid_2019_AMJ,t14_stderr_resid_2019_JAS,t14_stderr_resid_2019_OND,t14_stderr_resid_2020_JFM]; 

% plot mean seasonal enhancements and mean seasonal residuals 
figure()
t = tiledlayout(2,2,"TileSpacing","compact"); 
nexttile
plot(t01_seasonal_obs_enhc,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
hold on 
plot(t01_seasonal_model_enhc,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
yline(0)
grid on 
ylim([-4 1])
xlim([0 13])
xticklabels([])

nexttile
plot(t14_seasonal_obs_enhc,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
hold on 
plot(t14_seasonal_model_enhc,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
yline(0)
grid on 
xlim([0 13])
ylim([-4 1])
xticklabels([])

nexttile
errorbar(t01_seasonal_mean_resid,2.*t01_seasonal_stderr_resid,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
hold on 
yline(0)
xlim([0 13])
grid on 
ylim([-4 2])
xticks([1 2 3 4 5 6 7 8 9 10 11 12])
xticklabels({'2017 AMJ', '2017 JAS', '2017 OND', '2018 JFM','2018 AMJ','2018 JAS','2018 OND','2019 JFM','2019 AMJ','2019 JAS','2019 OND','2020 JFM'})

nexttile
errorbar(t14_seasonal_mean_resid,2.*t14_seasonal_stderr_resid,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
hold on 
yline(0)
xlim([0 13])
ylim([-4 2])
grid on 
xticks([1 2 3 4 5 6 7 8 9 10 11 12])
xticklabels({'2017 AMJ', '2017 JAS', '2017 OND', '2018 JFM','2018 AMJ','2018 JAS','2018 OND','2019 JFM','2019 AMJ','2019 JAS','2019 OND','2020 JFM'})

%% calculate residuals by year - in 12 month sections %% 
% residuals here are in 12 month sections, not by calendar year 

% first 12 month section is april 2017-march 2018
t01_mean_resid_year1_p1 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)); 
t01_mean_resid_year1_p2 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t01_mean_resid_year1 = mean(vertcat(t01_mean_resid_year1_p1,t01_mean_resid_year1_p2),'omitnan'); 

t14_mean_resid_year1_p1 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2017)); 
t14_mean_resid_year1_p2 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t14_mean_resid_year1 = mean(vertcat(t14_mean_resid_year1_p1,t14_mean_resid_year1_p2),'omitnan'); 

t01_mean_resid_year2_p1 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6 | month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9 | month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))); 
t01_mean_resid_year2_p2 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t01_mean_resid_year2 = mean(vertcat(t01_mean_resid_year2_p1,t01_mean_resid_year2_p2),'omitnan'); 

t14_mean_resid_year2_p1 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2018 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6 | month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9 | month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))); 
t14_mean_resid_year2_p2 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t14_mean_resid_year2 = mean(vertcat(t14_mean_resid_year2_p1,t14_mean_resid_year2_p2),'omitnan'); 

t01_mean_resid_year3_p1 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6 | month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9 | month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))); 
t01_mean_resid_year3_p2 = resid_t01_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t01_mean_resid_year3 = mean(vertcat(t01_mean_resid_year3_p1,t01_mean_resid_year3_p2),'omitnan'); 

t14_mean_resid_year3_p1 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2019 & (month(all_time_days)==4 | month(all_time_days)==5 | month(all_time_days)==6 | month(all_time_days)==7 | month(all_time_days)==8 | month(all_time_days)==9 | month(all_time_days)==10 | month(all_time_days)==11 | month(all_time_days)==12))); 
t14_mean_resid_year3_p2 = resid_t14_t09_mean_afternoon_enhc(find(year(all_time_days)==2020 & (month(all_time_days)==1 | month(all_time_days)==2 | month(all_time_days)==3))); 
t14_mean_resid_year3 = mean(vertcat(t14_mean_resid_year3_p1,t14_mean_resid_year3_p2),'omitnan'); 

%% calculate mean absolute error by day, month, season, year %% 
% Note the MAE for the year is using the residuals from the previous section, 
% which is the residuals in 12 month blocks, not by calendar year 

% calculate MAE at different time scales for each tower seperately 
t01_seasonal_MAE = mean(abs(t01_seasonal_mean_resid),'omitnan'); 
t14_seasonal_MAE = mean(abs(t14_seasonal_mean_resid),'omitnan'); 

t01_monthly_MAE = mean(abs(resid_t01_t09_monthly_mean),'omitnan'); 
t14_monthly_MAE = mean(abs(resid_t14_t09_monthly_mean),'omitnan'); 

t01_daily_MAE = mean(abs(resid_t01_t09_mean_afternoon_enhc),'omitnan'); 
t14_daily_MAE = mean(abs(resid_t14_t09_mean_afternoon_enhc),'omitnan'); 

% calculate MAE at different time scales - each tower data combined
seasonal_MAE = mean(abs([t01_seasonal_mean_resid,t14_seasonal_mean_resid]),'omitnan'); 
monthly_MAE = mean(abs(vertcat(resid_t01_t09_monthly_mean,resid_t14_t09_monthly_mean)),'omitnan'); 
daily_MAE = mean(abs(vertcat(resid_t01_t09_mean_afternoon_enhc,resid_t14_t09_mean_afternoon_enhc)),'omitnan'); 
year_MAE = mean(abs([t01_mean_resid_year1,t01_mean_resid_year2,t01_mean_resid_year3,t14_mean_resid_year1,t14_mean_resid_year2,t14_mean_resid_year3]),'omitnan'); 

%% figure of mean absolute errors %% 

figure()
plot(1:4,[daily_MAE,monthly_MAE,seasonal_MAE,year_MAE],'--o','MarkerEdgeColor','black','MarkerFaceColor','black','Color','black')
grid on 
ylim([0 2.5])
xlim([0.5 4.5])
xticks([1 2 3 4])
axis square
xticklabels({'Day','Month','Season','Year'})
ylabel({'Mean Abosolute', 'Error [ppm]'})
xlabel('Averaging Period')

%% plot the month mean enhancements (model and obs) for all months %% 

% make a datetime vector for all months in the time series 
all_months_series = datetime(2017,1,1):calmonths(1):datetime(2020,12,1); 
all_months_series = all_months_series'; 

figure()
t = tiledlayout(2,2,"TileSpacing","tight"); 
nexttile 
plot(all_months_series,co2_enhc_obs_t01_t09_monthly_mean,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
xlim([datetime(2017,1,1) datetime(2020,4,1)])
hold on 
plot(all_months_series,co2_enhc_model_t01_t09_monthly_mean,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
patch([datetime(2017,7,1),datetime(2017,9,31),datetime(2017,9,31),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,9,31),datetime(2018,9,31),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,9,31),datetime(2019,9,13),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
ylim([-6 6])
text(30,5,'(a)')
xticklabels([])
yline(0)
grid on 
ylabel({'Mean Monthly CO_2','Enhancement [ppm]'})
title('Tower 01')

nexttile 
plot(all_months_series,co2_enhc_obs_t14_t09_monthly_mean,'-o','MarkerFaceColor', 'black','MarkerEdgeColor','black','color','black','MarkerSize',3,'LineWidth',1.7)
xlim([datetime(2017,1,1) datetime(2020,4,1)])
hold on 
plot(all_months_series,co2_enhc_model_t14_t09_monthly_mean,'-o','MarkerFaceColor', 'red','MarkerEdgeColor','red','color','red','MarkerSize',3,'LineWidth',1.7)
patch([datetime(2017,7,1),datetime(2017,9,31),datetime(2017,9,31),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,9,31),datetime(2018,9,31),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,9,31),datetime(2019,9,13),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
ylim([-6 6])
legend('Observation','Model')
text(30,5,'(b)')
xticklabels([])
yline(0)
grid on 
title('Tower 14')

nexttile
errorbar(all_months_series,resid_t01_t09_monthly_mean,2.*resid_t01_t09_monthly_stderr,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
patch([datetime(2017,7,1),datetime(2017,9,31),datetime(2017,9,31),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,9,31),datetime(2018,9,31),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,9,31),datetime(2019,9,13),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
xlim([datetime(2017,1,1) datetime(2020,4,1)])
ylim([-9 6])
text(30,4.5,'(c)')
%ylabel('Mean Monthly Residual [ppm]')
ylabel({'Mean Monthly','Residual [ppm]'})
yline(0)
grid on 
xtickangle(30)

nexttile 
errorbar(all_months_series,resid_t14_t09_monthly_mean,2.*resid_t14_t09_monthly_stderr,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',4,'Color','black')
xlim([datetime(2017,1,1) datetime(2020,4,1)])
patch([datetime(2017,7,1),datetime(2017,9,31),datetime(2017,9,31),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,9,31),datetime(2018,9,31),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,9,31),datetime(2019,9,13),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
ylim([-9 6])
text(30,4.5,'(d)')
yline(0)
grid on 
xtickangle(30)

length(find(resid_t01_t09_monthly_mean<=0.5))
length(find(resid_t14_t09_monthly_mean<=0.5))


length(find(resid_t01_t09_monthly_mean<=1))
length(find(resid_t14_t09_monthly_mean<=1))


sum(~isnan(resid_t01_t09_monthly_mean))
sum(~isnan(resid_t14_t09_monthly_mean))



%% plot afternoon enhancements for model and obs for all time 

% plot the modeled and observed differences between towers - afternoon
% average 
figure()
t = tiledlayout(2,2,"TileSpacing",'tight');%"compact"); 
nexttile
plot(all_time_days,t01_t09_obs_afternoon_co2_enhc_mean,'.','MarkerEdgeColor','black','MarkerSize',5)
hold on 
plot(all_time_days,t01_t09_model_afternoon_co2_enhc_mean,'.','MarkerEdgeColor','red','MarkerSize',5)
yline(0)
%xline(datetime(2017,6,1),'--')
%xline(datetime(2018,6,1),'--')
%xline(datetime(2019,6,1),'--')
%xline(datetime(2017,9,1),'--')
%xline(datetime(2018,9,1),'--')
%xline(datetime(2019,9,1),'--')
patch([datetime(2017,7,1),datetime(2017,10,1),datetime(2017,10,1),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,10,1),datetime(2018,10,1),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,10,1),datetime(2019,10,1),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
%set(gca,'Children',flipud(get(gca,'Children')))
ylim([-30 20])
grid on 
legend('Observation','Model')
ylabel('CO_2 Enhancement [ppm]')
text(18,17,'(a)')
title('Tower 01')
set(gca,'Xticklabel',[])
%set(gca,'FontSize',24)
nexttile
plot(all_time_days,t14_t09_obs_afternoon_co2_enhc_mean,'.','MarkerEdgeColor','black','MarkerSize',5)
hold on 
plot(all_time_days,t14_t09_model_afternoon_co2_enhc_mean,'.','MarkerEdgeColor','red','MarkerSize',5)
patch([datetime(2017,7,1),datetime(2017,10,1),datetime(2017,10,1),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,10,1),datetime(2018,10,1),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,10,1),datetime(2019,10,1),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
yline(0)
ylim([-30 20])
text(18,17,'(b)')
grid on 
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend('Observation','Model','FontSize',16)
%xlabel('Time (UTC)')
%ylabel('CO_2 Enhancement [ppm]')
title('Tower 14')
%set(gca,'FontSize',24)

nexttile
plot(all_time_days,resid_t01_t09_mean_afternoon_enhc,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',3)
hold on 
yline(0)
ylim([-20 20])
text(18,17,'(c)')
grid on 
ylabel('Residual [ppm]')
patch([datetime(2017,7,1),datetime(2017,10,1),datetime(2017,10,1),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,10,1),datetime(2018,10,1),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,10,1),datetime(2019,10,1),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
%title('Tower 1')
xtickangle(30)
%set(gca,'FontSize',24)
nexttile
plot(all_time_days,resid_t14_t09_mean_afternoon_enhc,'^','MarkerEdgeColor','black','MarkerFaceColor','black','MarkerSize',3)
hold on 
yline(0)
grid on 
patch([datetime(2017,7,1),datetime(2017,10,1),datetime(2017,10,1),datetime(2017,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2018,7,1),datetime(2018,10,1),datetime(2018,10,1),datetime(2018,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
patch([datetime(2019,7,1),datetime(2019,10,1),datetime(2019,10,1),datetime(2019,7,1)],[-30,-30,20,20],[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeAlpha',0)
ylim([-20 20])
text(18,17,'(d)')%,'FontSize',16)
set(gca,'Yticklabel',[])
xtickangle(30)
%xlabel('Time (UTC)')
%ylabel('CO_2 Enhancement [ppm]')
%title('Tower 14')
%set(gca,'FontSize',24)


%% calculate & plot Depletion of CO2 mole frac by PFT (vegfrac included) %% 
% code in this section was originally from a script I wrote called
% "vegfrac_breakdown_plots.m"
% this section is the one done for the NEE with the vegfrac included 

% this analysis is only done for july 2018 

% load in the model outputs 
% biogenic co2 mole frac from corn 
co2_CORN_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Corn_Tower_01H3_2018-07.csv'); 
co2_CORN_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Corn_Tower_09_2018-07.csv'); 
co2_CORN_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Corn_Tower_14_2018-07.csv'); 
% biogenic co2 mole frac from other crops 
co2_OTHER_CROPS_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Other_Crops_Tower_01H3_2018-07.csv'); 
co2_OTHER_CROPS_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Other_Crops_Tower_09_2018-07.csv'); 
co2_OTHER_CROPS_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Other_Crops_Tower_14_2018-07.csv'); 
% biogenic co2 mole frac from deciduous broadleaf forest 
co2_DBF_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Deciduous_Forest_Tower_01H3_2018-07.csv'); 
co2_DBF_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Deciduous_Forest_Tower_09_2018-07.csv'); 
co2_DBF_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs/biogenic_co2_molefrac_Deciduous_Forest_Tower_14_2018-07.csv'); 

% crop the first row off all model outputs (this is just the headers)
co2_CORN_t01(1,:) = []; 
co2_CORN_t09(1,:) = [];
co2_CORN_t14(1,:) = [];
co2_OTHER_CROPS_t01(1,:) = []; 
co2_OTHER_CROPS_t09(1,:) = [];
co2_OTHER_CROPS_t14(1,:) = [];
co2_DBF_t01(1,:) = []; 
co2_DBF_t09(1,:) = [];
co2_DBF_t14(1,:) = [];

% set the dates - this should match the influence function format 
hours_jul2018(:,1) = datetime(2018,7,1,1,0,0):hours(1):datetime(2018,8,1,0,0,0); 

% quick plot of the model outputs 
figure()
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(hours_jul2018,co2_CORN_t01(:,2)); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,co2_CORN_t09(:,2)); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,co2_CORN_t14(:,2)); 
set(gca,'FontSize',14)
xlabel('Time (UTC)')
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)

% grab just the july 2018 data from one of the hourly obs/model time series
length(all_time_model)
length(t01_obs_co2)
length(t01_model_co2)

% get index for july 2018 in the time series with all hours available 
ind_jul_2018 = find(year(all_time_model)==2018 & month(all_time_model)==7); 

% check that the times grabbed make sense and assign to a new vector 
all_time_july2018 = all_time_model(ind_jul_2018); 

% grab the obs and the total model for the month of july 2018 
t01_obs_co2_july2018 = t01_obs_co2(ind_jul_2018); 
t09_obs_co2_july2018 = t09_obs_co2(ind_jul_2018); 
t14_obs_co2_july2018 = t14_obs_co2(ind_jul_2018); 

t01_model_co2_july2018 = t01_model_co2(ind_jul_2018); 
t09_model_co2_july2018 = t09_model_co2(ind_jul_2018); 
t14_model_co2_july2018 = t14_model_co2(ind_jul_2018); 

% kill the first hour (0 UTC) to match the NEE breakdown model outputs (no
% hour 0 UTC in that data .. would need to run all of june to get it)
all_time_july2018(1) = []; 

t01_obs_co2_july2018(1) = [];  
t09_obs_co2_july2018(1) = [];  
t14_obs_co2_july2018(1) = [];  

t01_model_co2_july2018(1) = []; 
t09_model_co2_july2018(1) = []; 
t14_model_co2_july2018(1) = []; 

% kill the last hour of the NEE breakdown model outputs (this is 0 UTC of
% august 2018)
hours_jul2018(end) = []; 
length(hours_jul2018)

co2_CORN_t01(end,:) = []; 
co2_CORN_t09(end,:) = [];
co2_CORN_t14(end,:) = [];
co2_OTHER_CROPS_t01(end,:) = []; 
co2_OTHER_CROPS_t09(end,:) = [];
co2_OTHER_CROPS_t14(end,:) = [];
co2_DBF_t01(end,:) = []; 
co2_DBF_t09(end,:) = [];
co2_DBF_t14(end,:) = [];

% check that the length is right 
length(co2_CORN_t01)

% now kill the first column of the monthly breakdowns (this is just the
% point number column - dont need anymore)
co2_CORN_t01(:,1) = []; 
co2_CORN_t09(:,1) = [];
co2_CORN_t14(:,1) = [];
co2_OTHER_CROPS_t01(:,1) = []; 
co2_OTHER_CROPS_t09(:,1) = [];
co2_OTHER_CROPS_t14(:,1) = [];
co2_DBF_t01(:,1) = []; 
co2_DBF_t09(:,1) = [];
co2_DBF_t14(:,1) = [];

% need to filter out the points that are removed (for wind direction or
% missing nans) from the obs and the total model - filter the NEE breakdown
% for july 2018 
for i=1:length(t01_obs_co2_july2018)   
    if isnan(t01_obs_co2_july2018(i)) == true % only using one tower b/c we already filtered all towers to be the same 
        
        co2_CORN_t01(i) = NaN; 
        co2_CORN_t09(i) = NaN;
        co2_CORN_t14(i) = NaN;
        co2_OTHER_CROPS_t01(i) = NaN; 
        co2_OTHER_CROPS_t09(i) = NaN;
        co2_OTHER_CROPS_t14(i) = NaN;
        co2_DBF_t01(i) = NaN; 
        co2_DBF_t09(i) = NaN;
        co2_DBF_t14(i) = NaN;

    end     
end 

%%%%%% calculate the afternoon average biogenic co2 for each PFT %%%%%%
% find the afternoon hours 
ind_afternoon_july2018 = find(hour(all_time_july2018)>= 17 & hour(all_time_july2018) <= 22);

% grab the datetime format hours for the afternoon 
afternoon_hours_july2018 = all_time_july2018(ind_afternoon_july2018); 

% reshape the afternoon hours to make sure everything looks good 
reshape_afternoon_hours_july2018 = reshape(afternoon_hours_july2018,6,[]); 

% grab the PFT specific co2 concentrations for july 2018 afternoon 
co2_CORN_t01_afternoon = co2_CORN_t01(ind_afternoon_july2018);
co2_CORN_t09_afternoon = co2_CORN_t09(ind_afternoon_july2018);
co2_CORN_t14_afternoon = co2_CORN_t14(ind_afternoon_july2018);
co2_OTHER_CROPS_t01_afternoon = co2_OTHER_CROPS_t01(ind_afternoon_july2018);
co2_OTHER_CROPS_t09_afternoon = co2_OTHER_CROPS_t09(ind_afternoon_july2018);
co2_OTHER_CROPS_t14_afternoon = co2_OTHER_CROPS_t14(ind_afternoon_july2018);
co2_DBF_t01_afternoon = co2_DBF_t01(ind_afternoon_july2018);
co2_DBF_t09_afternoon = co2_DBF_t09(ind_afternoon_july2018);
co2_DBF_t14_afternoon = co2_DBF_t14(ind_afternoon_july2018);

% now calculate the afternoon averages for the PFT co2 breakdown 
co2_CORN_t01_afternoon_avg(:,1) = mean(reshape(co2_CORN_t01_afternoon,6,[]),'omitnan'); 
co2_CORN_t09_afternoon_avg(:,1) = mean(reshape(co2_CORN_t09_afternoon,6,[]),'omitnan'); 
co2_CORN_t14_afternoon_avg(:,1) = mean(reshape(co2_CORN_t14_afternoon,6,[]),'omitnan'); 

co2_OTHER_CROPS_t01_afternoon_avg(:,1) = mean(reshape(co2_OTHER_CROPS_t01_afternoon,6,[]),'omitnan'); 
co2_OTHER_CROPS_t09_afternoon_avg(:,1) = mean(reshape(co2_OTHER_CROPS_t09_afternoon,6,[]),'omitnan'); 
co2_OTHER_CROPS_t14_afternoon_avg(:,1) = mean(reshape(co2_OTHER_CROPS_t14_afternoon,6,[]),'omitnan'); 

co2_DBF_t01_afternoon_avg(:,1) = mean(reshape(co2_DBF_t01_afternoon,6,[]),'omitnan'); 
co2_DBF_t09_afternoon_avg(:,1) = mean(reshape(co2_DBF_t09_afternoon,6,[]),'omitnan'); 
co2_DBF_t14_afternoon_avg(:,1) = mean(reshape(co2_DBF_t14_afternoon,6,[]),'omitnan'); 

% make a datetime vector of days in July 2018 (1 afternoon avg per day)
% make the hour = 1, hour doesn't really matter just do it 
days_jul2018(:,1) = datetime(2018,7,1):days(1):datetime(2018,7,31); 

% now, need to calculate the afternoon average CO2 for each tower
% seperately (for total model and obs) for the month of july 2018 - before
% i calculated the afternoon average of the differences between towers
% (enhancements)

t01_obs_co2_july2018_afternoon = t01_obs_co2_july2018(ind_afternoon_july2018); 
t09_obs_co2_july2018_afternoon = t09_obs_co2_july2018(ind_afternoon_july2018); 
t14_obs_co2_july2018_afternoon = t14_obs_co2_july2018(ind_afternoon_july2018); 

t01_model_co2_july2018_afternoon = t01_model_co2_july2018(ind_afternoon_july2018); 
t09_model_co2_july2018_afternoon = t09_model_co2_july2018(ind_afternoon_july2018); 
t14_model_co2_july2018_afternoon = t14_model_co2_july2018(ind_afternoon_july2018); 

t01_obs_co2_july2018_afternoon_avg = mean(reshape(t01_obs_co2_july2018_afternoon,6,[]),'omitnan'); 
t09_obs_co2_july2018_afternoon_avg = mean(reshape(t09_obs_co2_july2018_afternoon,6,[]),'omitnan'); 
t14_obs_co2_july2018_afternoon_avg = mean(reshape(t14_obs_co2_july2018_afternoon,6,[]),'omitnan'); 

t01_model_co2_july2018_afternoon_avg = mean(reshape(t01_model_co2_july2018_afternoon,6,[]),'omitnan'); 
t09_model_co2_july2018_afternoon_avg = mean(reshape(t09_model_co2_july2018_afternoon,6,[]),'omitnan');
t14_model_co2_july2018_afternoon_avg = mean(reshape(t14_model_co2_july2018_afternoon,6,[]),'omitnan'); 

% now calculate monthly mean of the afternoon average values for the obs,
% total model, and PFT breakdown model for all three towers 
CORN_t01_afternoon_avg_MONTHMEAN = mean(co2_CORN_t01_afternoon_avg,'omitnan'); 
CORN_t09_afternoon_avg_MONTHMEAN = mean(co2_CORN_t09_afternoon_avg,'omitnan'); 
CORN_t14_afternoon_avg_MONTHMEAN = mean(co2_CORN_t14_afternoon_avg,'omitnan'); 

OTHER_CROPS_t01_afternoon_MONTHMEAN = mean(co2_OTHER_CROPS_t01_afternoon_avg,'omitnan'); 
OTHER_CROPS_t09_afternoon_MONTHMEAN = mean(co2_OTHER_CROPS_t09_afternoon_avg,'omitnan'); 
OTHER_CROPS_t14_afternoon_MONTHMEAN = mean(co2_OTHER_CROPS_t14_afternoon_avg,'omitnan'); 

DBF_t01_afternoon_avg_MONTHMEAN = mean(co2_DBF_t01_afternoon_avg,'omitnan'); 
DBF_t09_afternoon_avg_MONTHMEAN = mean(co2_DBF_t09_afternoon_avg,'omitnan'); 
DBF_t14_afternoon_avg_MONTHMEAN = mean(co2_DBF_t14_afternoon_avg,'omitnan'); 

t01_obs_co2_july2018_afternoon_avg_MONTHMEAN = mean(t01_obs_co2_july2018_afternoon_avg,'omitnan'); 
t09_obs_co2_july2018_afternoon_avg_MONTHMEAN = mean(t09_obs_co2_july2018_afternoon_avg,'omitnan'); 
t14_obs_co2_july2018_afternoon_avg_MONTHMEAN = mean(t14_obs_co2_july2018_afternoon_avg,'omitnan'); 

t01_model_co2_july2018_afternoon_avg_MONTHMEAN = mean(t01_model_co2_july2018_afternoon_avg,'omitnan'); 
t09_model_co2_july2018_afternoon_avg_MONTHMEAN = mean(t09_model_co2_july2018_afternoon_avg,'omitnan');
t14_model_co2_july2018_afternoon_avg_MONTHMEAN = mean(t14_model_co2_july2018_afternoon_avg,'omitnan'); 

% let's try to make all the bars on one plot but color code things
bar_grouping = [t01_model_co2_july2018_afternoon_avg_MONTHMEAN ...
    t09_model_co2_july2018_afternoon_avg_MONTHMEAN ... 
    t14_model_co2_july2018_afternoon_avg_MONTHMEAN; ...
    CORN_t01_afternoon_avg_MONTHMEAN CORN_t09_afternoon_avg_MONTHMEAN ...
    CORN_t14_afternoon_avg_MONTHMEAN ; ... 
    OTHER_CROPS_t01_afternoon_MONTHMEAN OTHER_CROPS_t09_afternoon_MONTHMEAN ...
    OTHER_CROPS_t14_afternoon_MONTHMEAN; ... 
    DBF_t01_afternoon_avg_MONTHMEAN DBF_t09_afternoon_avg_MONTHMEAN ...
    DBF_t14_afternoon_avg_MONTHMEAN]; 

% new colors options 
figure()
b = bar(bar_grouping); 
grid on 
b(1).FaceColor = '#00345d'; 
b(2).FaceColor = '#b2614c'; 
b(3).FaceColor = '#fef7a9'; 
xticklabels({'Model Total','Corn','Other Crops','DBF'})
ylabel({'Depletion of CO_2 [ppm]'},'FontSize',12)
legend('Tower 01','Tower 09','Tower 14','FontSize',12)
axis square
set(gca,'FontSize',12)

%% calculate & plot Depletion of CO2 mole frac by PFT (vegfrac NOT included) %% 
% this section is the one done for the NEE with the vegfrac NOT included 

% load in the model outputs 
% biogenic co2 mole frac from corn 
CORN_nofrac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Corn_NOVEGFRAC_Tower_01H3_2018-07.csv'); 
CORN_nofrac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Corn_NOVEGFRAC_Tower_09_2018-07.csv'); 
CORN_nofrac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Corn_NOVEGFRAC_Tower_14_2018-07.csv'); 
% biogenic co2 mole frac from other crops 
OTHER_CROPS_nofrac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Other_Crops_NOVEGFRAC_Tower_01H3_2018-07.csv'); 
OTHER_CROPS_nofrac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Other_Crops_NOVEGFRAC_Tower_09_2018-07.csv'); 
OTHER_CROPS_nofrac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Other_Crops_NOVEGFRAC_Tower_14_2018-07.csv'); 
% biogenic co2 mole frac from deciduous broadleaf forest 
DBF_nofrac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Deciduous_Forest_NOVEGFRAC_Tower_01H3_2018-07.csv'); 
DBF_nofrac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Deciduous_Forest_NOVEGFRAC_Tower_09_2018-07.csv'); 
DBF_nofrac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_molefrac_outputs_PFTs_NOVEGFRAC/biogenic_co2_molefrac_Deciduous_Forest_NOVEGFRAC_Tower_14_2018-07.csv'); 

% crop the first row off all model outputs (this is just the headers)
CORN_nofrac_t01(1,:) = []; 
CORN_nofrac_t09(1,:) = [];
CORN_nofrac_t14(1,:) = [];
OTHER_CROPS_nofrac_t01(1,:) = []; 
OTHER_CROPS_nofrac_t09(1,:) = [];
OTHER_CROPS_nofrac_t14(1,:) = [];
DBF_nofrac_t01(1,:) = []; 
DBF_nofrac_t09(1,:) = [];
DBF_nofrac_t14(1,:) = [];

% now let's make the model and obs agree in length (kill the last hour of
% the model (this hour is august 1, 0 UTC, 2018)
length(hours_jul2018)

CORN_nofrac_t01(end,:) = []; 
CORN_nofrac_t09(end,:) = [];
CORN_nofrac_t14(end,:) = [];
OTHER_CROPS_nofrac_t01(end,:) = []; 
OTHER_CROPS_nofrac_t09(end,:) = [];
OTHER_CROPS_nofrac_t14(end,:) = [];
DBF_nofrac_t01(end,:) = []; 
DBF_nofrac_t09(end,:) = [];
DBF_nofrac_t14(end,:) = [];

% check that the length is right 
length(CORN_nofrac_t01)

% now kill the first column of the monthly breakdowns (this is just the
% point number column / hr of month - dont need anymore)
CORN_nofrac_t01(:,1) = []; 
CORN_nofrac_t09(:,1) = [];
CORN_nofrac_t14(:,1) = [];
OTHER_CROPS_nofrac_t01(:,1) = []; 
OTHER_CROPS_nofrac_t09(:,1) = [];
OTHER_CROPS_nofrac_t14(:,1) = [];
DBF_nofrac_t01(:,1) = []; 
DBF_nofrac_t09(:,1) = [];
DBF_nofrac_t14(:,1) = [];

% make a quick plot of the outputs 
figure()
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(hours_jul2018,CORN_nofrac_t01); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,CORN_nofrac_t09); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,CORN_nofrac_t14); 
set(gca,'FontSize',14)
xlabel('Time (UTC)')
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)


figure()
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(hours_jul2018,OTHER_CROPS_nofrac_t01); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,OTHER_CROPS_nofrac_t09); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,OTHER_CROPS_nofrac_t14); 
set(gca,'FontSize',14)
xlabel('Time (UTC)')
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)

figure()
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(hours_jul2018,DBF_nofrac_t01); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,DBF_nofrac_t09); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
nexttile
plot(hours_jul2018,DBF_nofrac_t14); 
set(gca,'FontSize',14)
xlabel('Time (UTC)')
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)


% need to filter out the points that are removed (for wind direction or
% missing nans) from the obs and the total model - filter the PFT novegfrac 
% for july 2018 
for i=1:length(t01_obs_co2_july2018)   
    if isnan(t01_obs_co2_july2018(i)) == true % only using one tower b/c we already filtered all towers to be the same 
        
        CORN_nofrac_t01(i) = NaN; 
        CORN_nofrac_t09(i) = NaN;
        CORN_nofrac_t14(i) = NaN;
        OTHER_CROPS_nofrac_t01(i) = NaN; 
        OTHER_CROPS_nofrac_t09(i) = NaN;
        OTHER_CROPS_nofrac_t14(i) = NaN;
        DBF_nofrac_t01(i) = NaN; 
        DBF_nofrac_t09(i) = NaN;
        DBF_nofrac_t14(i) = NaN;

    end     
end 


% grab the PFT specific co2 concentrations for july 2018 afternoon 
CORN_nofrac_t01_afternoon = CORN_nofrac_t01(ind_afternoon_july2018);
CORN_nofrac_t09_afternoon = CORN_nofrac_t09(ind_afternoon_july2018);
CORN_nofrac_t14_afternoon = CORN_nofrac_t14(ind_afternoon_july2018);
OTHER_CROPS_nofrac_t01_afternoon = OTHER_CROPS_nofrac_t01(ind_afternoon_july2018);
OTHER_CROPS_nofrac_t09_afternoon = OTHER_CROPS_nofrac_t09(ind_afternoon_july2018);
OTHER_CROPS_nofrac_t14_afternoon = OTHER_CROPS_nofrac_t14(ind_afternoon_july2018);
DBF_nofrac_t01_afternoon = DBF_nofrac_t01(ind_afternoon_july2018);
DBF_nofrac_t09_afternoon = DBF_nofrac_t09(ind_afternoon_july2018);
DBF_nofrac_t14_afternoon = DBF_nofrac_t14(ind_afternoon_july2018);

% now calculate the afternoon averages for the PFT co2 breakdown 
CORN_nofrac_t01_afternoon_avg(:,1) = mean(reshape(CORN_nofrac_t01_afternoon,6,[]),'omitnan'); 
CORN_nofrac_t09_afternoon_avg(:,1) = mean(reshape(CORN_nofrac_t09_afternoon,6,[]),'omitnan'); 
CORN_nofrac_t14_afternoon_avg(:,1) = mean(reshape(CORN_nofrac_t14_afternoon,6,[]),'omitnan'); 

OTHER_CROPS_nofrac_t01_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_nofrac_t01_afternoon,6,[]),'omitnan'); 
OTHER_CROPS_nofrac_t09_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_nofrac_t09_afternoon,6,[]),'omitnan'); 
OTHER_CROPS_nofrac_t14_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_nofrac_t14_afternoon,6,[]),'omitnan'); 

DBF_nofrac_t01_afternoon_avg(:,1) = mean(reshape(DBF_nofrac_t01_afternoon,6,[]),'omitnan'); 
DBF_nofrac_t09_afternoon_avg(:,1) = mean(reshape(DBF_nofrac_t09_afternoon,6,[]),'omitnan'); 
DBF_nofrac_t14_afternoon_avg(:,1) = mean(reshape(DBF_nofrac_t14_afternoon,6,[]),'omitnan'); 

% quick plots of the afternoon averages 
figure()
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(days_jul2018,CORN_nofrac_t01_afternoon_avg); 
hold on 
plot(days_jul2018,CORN_nofrac_t09_afternoon_avg);
plot(days_jul2018,CORN_nofrac_t14_afternoon_avg);
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
legend('01','09','14')
nexttile
plot(days_jul2018,OTHER_CROPS_nofrac_t01_afternoon_avg); 
hold on 
plot(days_jul2018,OTHER_CROPS_nofrac_t09_afternoon_avg); 
plot(days_jul2018,OTHER_CROPS_nofrac_t14_afternoon_avg); 
set(gca,'Xticklabel',[])
set(gca,'FontSize',14)
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
legend('01','09','14')
nexttile
plot(days_jul2018,DBF_nofrac_t01_afternoon_avg); 
hold on 
plot(days_jul2018,DBF_nofrac_t09_afternoon_avg); 
plot(days_jul2018,DBF_nofrac_t14_afternoon_avg); 
set(gca,'FontSize',14)
xlabel('Time (UTC)')
ylabel('CO_2 [ppm]')
%ylim([-20 20])
yline(0)
legend('01','09','14')

% now caluclate the monthly means for each PFT (with no vegfrac)
CORN_nofrac_t01_afternoon_avg_MONTHMEAN = mean(CORN_nofrac_t01_afternoon_avg,'omitnan'); 
CORN_nofrac_t09_afternoon_avg_MONTHMEAN = mean(CORN_nofrac_t09_afternoon_avg,'omitnan'); 
CORN_nofrac_t14_afternoon_avg_MONTHMEAN = mean(CORN_nofrac_t14_afternoon_avg,'omitnan'); 

OTHER_CROPS_nofrac_t01_afternoon_MONTHMEAN = mean(OTHER_CROPS_nofrac_t01_afternoon_avg,'omitnan'); 
OTHER_CROPS_nofrac_t09_afternoon_MONTHMEAN = mean(OTHER_CROPS_nofrac_t09_afternoon_avg,'omitnan'); 
OTHER_CROPS_nofrac_t14_afternoon_MONTHMEAN = mean(OTHER_CROPS_nofrac_t14_afternoon_avg,'omitnan'); 

DBF_nofrac_t01_afternoon_avg_MONTHMEAN = mean(DBF_nofrac_t01_afternoon_avg,'omitnan'); 
DBF_nofrac_t09_afternoon_avg_MONTHMEAN = mean(DBF_nofrac_t09_afternoon_avg,'omitnan'); 
DBF_nofrac_t14_afternoon_avg_MONTHMEAN = mean(DBF_nofrac_t14_afternoon_avg,'omitnan'); 


% let's try to make all the bars on one plot but color code things
bar_grouping_PFT_nofrac = [CORN_nofrac_t01_afternoon_avg_MONTHMEAN CORN_nofrac_t09_afternoon_avg_MONTHMEAN ...
    CORN_nofrac_t14_afternoon_avg_MONTHMEAN ; ... 
    OTHER_CROPS_nofrac_t01_afternoon_MONTHMEAN OTHER_CROPS_nofrac_t09_afternoon_MONTHMEAN ...
    OTHER_CROPS_nofrac_t14_afternoon_MONTHMEAN; ... 
    DBF_nofrac_t01_afternoon_avg_MONTHMEAN DBF_nofrac_t09_afternoon_avg_MONTHMEAN ...
    DBF_nofrac_t14_afternoon_avg_MONTHMEAN]; 

figure()
b = bar(bar_grouping_PFT_nofrac); 
grid on 
% b(1).FaceColor = 
xticklabels({'Corn','Other Crops','DBF'})
ylabel('Unweighted Average July 2018 Afternoon Biogenic CO_2 [ppm]','FontSize',12)
legend('Tower 01','Tower 09','Tower 14')



%% Vegfrac and influence function convolution analysis 
% MAKE SURE TO CHANGE THE NAMES HERE 
CORN_frac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Corn_Tower_01H3_2018-07.csv'); 
CORN_frac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Corn_Tower_09_2018-07.csv'); 
CORN_frac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Corn_Tower_14_2018-07.csv'); 

OTHER_CROPS_frac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Other_Crops_Tower_01H3_2018-07.csv'); 
OTHER_CROPS_frac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Other_Crops_Tower_09_2018-07.csv'); 
OTHER_CROPS_frac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Other_Crops_Tower_14_2018-07.csv'); 
 
DBF_frac_t01 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Deciduous_Forest_Tower_01H3_2018-07.csv'); 
DBF_frac_t09 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Deciduous_Forest_Tower_09_2018-07.csv'); 
DBF_frac_t14 = readmatrix('/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/convolution_outputs_VEGFRAC/conv_vegfrac_Deciduous_Forest_Tower_14_2018-07.csv'); 


% crop the first row off all model outputs (this is just the headers)
CORN_frac_t01(1,:) = []; 
CORN_frac_t09(1,:) = [];
CORN_frac_t14(1,:) = [];
OTHER_CROPS_frac_t01(1,:) = []; 
OTHER_CROPS_frac_t09(1,:) = [];
OTHER_CROPS_frac_t14(1,:) = [];
DBF_frac_t01(1,:) = []; 
DBF_frac_t09(1,:) = [];
DBF_frac_t14(1,:) = [];


% now kill the first column of the monthly vegfrac convolutions (this is just the
% point number column - dont need anymore)
CORN_frac_t01(:,1) = []; 
CORN_frac_t09(:,1) = [];
CORN_frac_t14(:,1) = [];
OTHER_CROPS_frac_t01(:,1) = []; 
OTHER_CROPS_frac_t09(:,1) = [];
OTHER_CROPS_frac_t14(:,1) = [];
DBF_frac_t01(:,1) = []; 
DBF_frac_t09(:,1) = [];
DBF_frac_t14(:,1) = [];

% kill the last hour of the vegfrac conv model outputs (this is 0 UTC of
% august 2018)
CORN_frac_t01(end,:) = []; 
CORN_frac_t09(end,:) = [];
CORN_frac_t14(end,:) = [];
OTHER_CROPS_frac_t01(end,:) = []; 
OTHER_CROPS_frac_t09(end,:) = [];
OTHER_CROPS_frac_t14(end,:) = [];
DBF_frac_t01(end,:) = []; 
DBF_frac_t09(end,:) = [];
DBF_frac_t14(end,:) = [];

% check that the length is right 
length(CORN_frac_t01)
length(hours_jul2018)

% need to filter out the points that are removed (for wind direction or
% missing nans) from the obs and the total model - filter the vegfrac
% breakdown 
for i=1:length(t01_obs_co2_july2018)   
    if isnan(t01_obs_co2_july2018(i)) == true % only using one tower b/c we already filtered all towers to be the same 
        
        CORN_frac_t01(i) = NaN; 
        CORN_frac_t09(i) = NaN;
        CORN_frac_t14(i) = NaN;
        OTHER_CROPS_frac_t01(i) = NaN; 
        OTHER_CROPS_frac_t09(i) = NaN;
        OTHER_CROPS_frac_t14(i) = NaN;
        DBF_frac_t01(i) = NaN; 
        DBF_frac_t09(i) = NaN;
        DBF_frac_t14(i) = NaN;

    end     
end 

% okay now I need to grab the afternoons (this is the same line as in the
% above section - just pasted here incase it is run seperately)
ind_afternoon_july2018 = find(hour(all_time_july2018)>= 17 & hour(all_time_july2018) <= 22);

% grab the datetime format hours for the afternoon (this is the same line as in the
% above section - just pasted here incase it is run seperately)
afternoon_hours_july2018 = all_time_july2018(ind_afternoon_july2018); 

% grab the vegfrac conv model outputs for july 2018 afternoon 
CORN_frac_t01_afternoon = CORN_frac_t01(ind_afternoon_july2018);
CORN_frac_t09_afternoon = CORN_frac_t09(ind_afternoon_july2018);
CORN_frac_t14_afternoon = CORN_frac_t14(ind_afternoon_july2018);
OTHER_CROPS_frac_t01_afternoon = OTHER_CROPS_frac_t01(ind_afternoon_july2018);
OTHER_CROPS_frac_t09_afternoon = OTHER_CROPS_frac_t09(ind_afternoon_july2018);
OTHER_CROPS_frac_t14_afternoon = OTHER_CROPS_frac_t14(ind_afternoon_july2018);
DBF_frac_t01_afternoon = DBF_frac_t01(ind_afternoon_july2018);
DBF_frac_t09_afternoon = DBF_frac_t09(ind_afternoon_july2018);
DBF_frac_t14_afternoon = DBF_frac_t14(ind_afternoon_july2018);

% now I think instead of averaging I want the sum here ... (i did in
% average in my old seperate code)
CORN_frac_t01_afternoon_sum = sum(CORN_frac_t01_afternoon,'omitnan');
CORN_frac_t09_afternoon_sum = sum(CORN_frac_t09_afternoon,'omitnan');
CORN_frac_t14_afternoon_sum = sum(CORN_frac_t14_afternoon,'omitnan');
OTHER_CROPS_frac_t01_afternoon_sum = sum(OTHER_CROPS_frac_t01_afternoon,'omitnan');
OTHER_CROPS_frac_t09_afternoon_sum = sum(OTHER_CROPS_frac_t09_afternoon,'omitnan');
OTHER_CROPS_frac_t14_afternoon_sum = sum(OTHER_CROPS_frac_t14_afternoon,'omitnan');
DBF_frac_t01_afternoon_sum = sum(DBF_frac_t01_afternoon,'omitnan');
DBF_frac_t09_afternoon_sum = sum(DBF_frac_t09_afternoon,'omitnan');
DBF_frac_t14_afternoon_sum = sum(DBF_frac_t14_afternoon,'omitnan');
 


% let's try to make all the bars on one plot but color code things
bar_grouping_vegfrac = [CORN_frac_t01_afternoon_sum ...
    CORN_frac_t09_afternoon_sum CORN_frac_t14_afternoon_sum ; ...
    OTHER_CROPS_frac_t01_afternoon_sum OTHER_CROPS_frac_t09_afternoon_sum ...
    OTHER_CROPS_frac_t14_afternoon_sum ; DBF_frac_t01_afternoon_sum ...
    DBF_frac_t09_afternoon_sum DBF_frac_t14_afternoon_sum ]; 

figure()
b = bar(bar_grouping_vegfrac); 
grid on 
b(1).FaceColor = '#B7E6A5'; 
b(2).FaceColor = '#089099'; 
b(3).FaceColor = '#045275'; 
xticklabels({'Corn','Other Crops','DBF'})
ylabel('Contribution to influence [ppb mol^{-1} km^2 hr]','FontSize',24) 
legend('Tower 01','Tower 09','Tower 14','FontSize',16)
axis square
set(gca,'FontSize',24)


%% "total" influence function analysis 
% set path to the total influence function csv files 
path_in_total_inf = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/inf_fun_sum_outputs/';

% read in the vector of total hourly influence function values for each
% hour of july 2018 
hr_sum_t01 = readmatrix(strcat(path_in_total_inf,'inf_func_sum_Tower_01H3_2018-07.csv')); 
hr_sum_t09 = readmatrix(strcat(path_in_total_inf,'inf_func_sum_Tower_09_2018-07.csv')); 
hr_sum_t14 = readmatrix(strcat(path_in_total_inf,'inf_func_sum_Tower_14_2018-07.csv')); 

% delete the first row of the inputs, these are column headers 
hr_sum_t01(1,:) = [];  
hr_sum_t09(1,:) = [];  
hr_sum_t14(1,:) = [];  

% now delete the last row of the inputs (this is hour 0 utc of aug 1 2018
% b/c that is how the influence functions are set up)
hr_sum_t01(end,:) = [];  
hr_sum_t09(end,:) = [];  
hr_sum_t14(end,:) = []; 

% take a quick look at the values for the whole month 

figure()
tiledlayout(3,1,'TileSpacing','compact')
nexttile
plot(hours_jul2018,hr_sum_t01(:,2))
%ylabel('influence function [ppb mol^{-1} km^2 hr]')
title('Tower 1')
ylim([0 2])
nexttile
plot(hours_jul2018,hr_sum_t09(:,2))
ylabel('influence function [ppb mol^{-1} km^2 hr]')
title('Tower 9')
ylim([0 2])
nexttile
plot(hours_jul2018,hr_sum_t14(:,2))
xlabel('Time (UTC)')
%ylabel('influence function [ppb mol^{-1} km^2 hr]')
title('Tower 14')
ylim([0 2])

% now sum over the whole month 
monthly_sum_t01 = sum(hr_sum_t01(:,2)); 
monthly_sum_t09 = sum(hr_sum_t09(:,2)); 
monthly_sum_t14 = sum(hr_sum_t14(:,2)); 

% now we need to filter out the influence function sum values for each hour
% when we have either no obs point / the obs point is filtered out 
for i=1:length(t01_obs_co2_july2018)   
    if isnan(t01_obs_co2_july2018(i)) == true % only using one tower b/c we already filtered all towers to be the same 
        hr_sum_t01(i,2) = NaN; 
        hr_sum_t09(i,2) = NaN;
        hr_sum_t14(i,2) = NaN;
    end     
end 

% let's pull out only the afternoon hours - find index  
afternoon_hours_sum_ind = find(hour(hours_jul2018) >= 17 & hour(hours_jul2018) <= 22); 

% grab just the afternoon hours of the sum values 
hr_sum_t01_afternoon = hr_sum_t01(afternoon_hours_sum_ind,2); 
hr_sum_t09_afternoon = hr_sum_t09(afternoon_hours_sum_ind,2); 
hr_sum_t14_afternoon = hr_sum_t14(afternoon_hours_sum_ind,2); 

% grab just the afternoon hours of the datetime vector 
hrs_jul2018_afternoon = hours_jul2018(afternoon_hours_sum_ind); 

% take a look at the just afternoon filtered values 
figure()
tiledlayout(3,1,'TileSpacing','compact')
nexttile
plot(hrs_jul2018_afternoon,hr_sum_t01_afternoon(:))
title('Tower 1')
ylim([0 2])
nexttile
plot(hrs_jul2018_afternoon,hr_sum_t09_afternoon(:))
ylabel('influence function [ppb mol^{-1} km^2 hr]')
title('Tower 9')
ylim([0 2])
nexttile
plot(hrs_jul2018_afternoon,hr_sum_t14_afternoon(:))
xlabel('Time (UTC)')
title('Tower 14')
ylim([0 2])


% now let's check the afternoon sum 
monthly_sum_t01_afternoon = sum(hr_sum_t01_afternoon(:),'omitnan'); 
monthly_sum_t09_afternoon = sum(hr_sum_t09_afternoon(:),'omitnan'); 
monthly_sum_t14_afternoon = sum(hr_sum_t14_afternoon(:),'omitnan'); 


% let's try to make the PFT fractional breakdown plot again but add the
% totals 

% let's try to make all the bars on one plot but color code things
bar_grouping_vegfrac_tot = [monthly_sum_t01_afternoon monthly_sum_t09_afternoon monthly_sum_t14_afternoon; ...
    CORN_frac_t01_afternoon_sum ...
    CORN_frac_t09_afternoon_sum CORN_frac_t14_afternoon_sum ; ...
    OTHER_CROPS_frac_t01_afternoon_sum OTHER_CROPS_frac_t09_afternoon_sum ...
    OTHER_CROPS_frac_t14_afternoon_sum ; DBF_frac_t01_afternoon_sum ...
    DBF_frac_t09_afternoon_sum DBF_frac_t14_afternoon_sum ;]; 

figure()
b = bar(bar_grouping_vegfrac_tot); 
grid on 
b(1).FaceColor = '#B7E6A5'; 
b(2).FaceColor = '#089099'; 
b(3).FaceColor = '#045275'; 
xticklabels({'Total','Corn','Other Crops','DBF'})
ylabel({'Contribution to influence function', '[ppb mol^{-1} km^2 hr]'},'FontSize',16) 
legend('Tower 01','Tower 09','Tower 14')
axis square
set(gca,'FontSize',12)

% figure with new colors
figure()
b = bar(bar_grouping_vegfrac_tot); 
grid on 
b(1).FaceColor = '#00345d'; 
b(2).FaceColor = '#b2614c'; 
b(3).FaceColor = '#fef7a9'; 
xticklabels({'Total','Corn','Other Crops','DBF'})
ylabel({'Contribution to influence', 'function [ppb mol^{-1} km^2 hr]'})%,'FontSize',16) 
legend('Tower 01','Tower 09','Tower 14')
axis square
%set(gca,'FontSize',12)

% now let's try to make the vegfrac breakdown but relative to the total
% influence function value 
bar_grouping_vegfrac_percent = [CORN_frac_t01_afternoon_sum/monthly_sum_t01_afternoon ...
    CORN_frac_t09_afternoon_sum/monthly_sum_t09_afternoon CORN_frac_t14_afternoon_sum/monthly_sum_t14_afternoon ; ...
    OTHER_CROPS_frac_t01_afternoon_sum/monthly_sum_t01_afternoon OTHER_CROPS_frac_t09_afternoon_sum/monthly_sum_t09_afternoon ...
    OTHER_CROPS_frac_t14_afternoon_sum/monthly_sum_t14_afternoon ; DBF_frac_t01_afternoon_sum/monthly_sum_t01_afternoon ...
    DBF_frac_t09_afternoon_sum/monthly_sum_t09_afternoon DBF_frac_t14_afternoon_sum/monthly_sum_t14_afternoon ;]; 

figure()
b = bar(bar_grouping_vegfrac_percent.*100); 
grid on 
% b(1).FaceColor = 
xticklabels({'Corn','Other Crops','DBF'})
ylabel('Contribution to influence function (%)','FontSize',12) 
legend('Tower 01','Tower 09','Tower 14')

%% look at the NEE productivity in each tower influence area 
% for each hour, divide out the impact of the inlufnece functions 

CORN_NEE_t01 = CORN_nofrac_t01./(hr_sum_t01(:,2).*(1/1000)); 
CORN_NEE_t09 = CORN_nofrac_t09./(hr_sum_t09(:,2).*(1/1000)); 
CORN_NEE_t14 = CORN_nofrac_t14./(hr_sum_t14(:,2)*(1/1000)); 

OTHER_CROPS_NEE_t01 = OTHER_CROPS_nofrac_t01./(hr_sum_t01(:,2).*(1/1000)); 
OTHER_CROPS_NEE_t09 = OTHER_CROPS_nofrac_t09./(hr_sum_t09(:,2).*(1/1000)); 
OTHER_CROPS_NEE_t14 = OTHER_CROPS_nofrac_t14./(hr_sum_t14(:,2).*(1/1000)); 

DBF_NEE_t01 = DBF_nofrac_t01./(hr_sum_t01(:,2).*(1/1000)); 
DBF_NEE_t09 = DBF_nofrac_t09./(hr_sum_t09(:,2).*(1/1000)); 
DBF_NEE_t14 = DBF_nofrac_t14./(hr_sum_t14(:,2).*(1/1000)); 

nee_convert_factor = (1e6)/(1000*1000*3600); 

CORN_NEE_t01 = CORN_NEE_t01.*nee_convert_factor; 
CORN_NEE_t09 = CORN_NEE_t09.*nee_convert_factor; 
CORN_NEE_t14 = CORN_NEE_t14.*nee_convert_factor; 

OTHER_CROPS_NEE_t01 = OTHER_CROPS_NEE_t01.*nee_convert_factor; 
OTHER_CROPS_NEE_t09 = OTHER_CROPS_NEE_t09.*nee_convert_factor; 
OTHER_CROPS_NEE_t14 = OTHER_CROPS_NEE_t14.*nee_convert_factor; 

DBF_NEE_t01 = DBF_NEE_t01.*nee_convert_factor; 
DBF_NEE_t09 = DBF_NEE_t09.*nee_convert_factor; 
DBF_NEE_t14 = DBF_NEE_t14.*nee_convert_factor; 

plot(CORN_NEE_t01,'.')
hold on 
plot(CORN_NEE_t09,'.')
plot(CORN_NEE_t14,'.')


CORN_NEE_t01_afternoon_avg(:,1) = mean(reshape(CORN_NEE_t01(ind_afternoon_july2018),6,[]),'omitnan'); 
CORN_NEE_t09_afternoon_avg(:,1) = mean(reshape(CORN_NEE_t09(ind_afternoon_july2018),6,[]),'omitnan'); 
CORN_NEE_t14_afternoon_avg(:,1) = mean(reshape(CORN_NEE_t14(ind_afternoon_july2018),6,[]),'omitnan'); 

CORN_NEE_t01_mean = mean(CORN_NEE_t01_afternoon_avg,'omitnan'); 
CORN_NEE_t09_mean = mean(CORN_NEE_t09_afternoon_avg,'omitnan'); 
CORN_NEE_t14_mean = mean(CORN_NEE_t14_afternoon_avg,'omitnan'); 

OTHER_CROPS_NEE_t01_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_NEE_t01(ind_afternoon_july2018),6,[]),'omitnan'); 
OTHER_CROPS_NEE_t09_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_NEE_t09(ind_afternoon_july2018),6,[]),'omitnan'); 
OTHER_CROPS_NEE_t14_afternoon_avg(:,1) = mean(reshape(OTHER_CROPS_NEE_t14(ind_afternoon_july2018),6,[]),'omitnan'); 

OTHER_CROPS_NEE_t01_mean = mean(OTHER_CROPS_NEE_t01_afternoon_avg,'omitnan'); 
OTHER_CROPS_NEE_t09_mean = mean(OTHER_CROPS_NEE_t09_afternoon_avg,'omitnan'); 
OTHER_CROPS_NEE_t14_mean = mean(OTHER_CROPS_NEE_t14_afternoon_avg,'omitnan'); 

DBF_NEE_t01_afternoon_avg(:,1) = mean(reshape(DBF_NEE_t01(ind_afternoon_july2018),6,[]),'omitnan'); 
DBF_NEE_t09_afternoon_avg(:,1) = mean(reshape(DBF_NEE_t09(ind_afternoon_july2018),6,[]),'omitnan'); 
DBF_NEE_t14_afternoon_avg(:,1) = mean(reshape(DBF_NEE_t14(ind_afternoon_july2018),6,[]),'omitnan'); 

DBF_NEE_t01_mean = mean(DBF_NEE_t01_afternoon_avg,'omitnan'); 
DBF_NEE_t09_mean = mean(DBF_NEE_t09_afternoon_avg,'omitnan'); 
DBF_NEE_t14_mean = mean(DBF_NEE_t14_afternoon_avg,'omitnan'); 

% let's try to make all the bars on one plot but color code things
bar_grouping_nee = [CORN_NEE_t01_mean CORN_NEE_t09_mean CORN_NEE_t14_mean;...
                    OTHER_CROPS_NEE_t01_mean OTHER_CROPS_NEE_t09_mean OTHER_CROPS_NEE_t14_mean;...
                    DBF_NEE_t01_mean DBF_NEE_t09_mean DBF_NEE_t14_mean]; 


figure()
b = bar(bar_grouping_nee); 
grid on 
b(1).FaceColor = '#B7E6A5'; 
b(2).FaceColor = '#089099'; 
b(3).FaceColor = '#045275'; 
xticklabels({'Corn','Other Crops','DBF'})
%ylabel({'Mean July 2018 Afternoon'; 'NEE [\mu mol m^{-2} s^{-1} ]'},'FontSize',16) 
ylabel({'Mean July 2018 Afternoon NEE [\mu mol m^{-2} s^{-1} ]'},'FontSize',16) 
legend('Tower 01','Tower 09','Tower 14')
axis square
set(gca,'FontSize',12)

% update colors 
figure()
b = bar(bar_grouping_nee); 
grid on 
b(1).FaceColor = '#00345d'; 
b(2).FaceColor = '#b2614c'; 
b(3).FaceColor = '#fef7a9'; 
xticklabels({'Corn','Other Crops','DBF'})
%ylabel({'Mean July 2018 Afternoon'; 'NEE [\mu mol m^{-2} s^{-1} ]'},'FontSize',16) 
ylabel({'NEE [\mumol m^{-2} s^{-1} ]'})%,'FontSize',16) 
legend('Tower 01','Tower 09','Tower 14')
axis square
%set(gca,'FontSize',12)

%% evaluate nighttime impact %% 
% determine if the influence functions are strongly influenced by nighttime 
% hours during the afternoon 

% set path to nighttime outputs
path_night = '/storage/group/zrb5027/default/INFLUX/smm8236/convolution_result_analysis/night_impact_outputs/'; 

% read in the files for the nighttime impact for each tower as tables 
night_t01_data = readtable(strcat(path_night,'night_impact_Tower_01H3_2018-07.csv'));
night_t09_data = readtable(strcat(path_night,'night_impact_Tower_09_2018-07.csv'));
night_t14_data = readtable(strcat(path_night,'night_impact_Tower_14_2018-07.csv'));

% kill the last row of the data - this is the first hour of aug 1 2018 
night_t01_data(end,:) = []; 
night_t09_data(end,:) = []; 
night_t14_data(end,:) = []; 

% let's just plot the data 
figure()
plot(hours_jul2018,night_t01_data{:,3})
hold on 
plot(hours_jul2018,night_t09_data{:,3})
plot(hours_jul2018,night_t14_data{:,3})
legend('T01','T09','T14')
title('Fraction of influence function impacted by night')

% now let's get just the afternoon hours 
hours_jul2018(ind_afternoon_july2018)

% grab the fractions in the afternoon 
night_t01_frac_afternoon = night_t01_data{ind_afternoon_july2018,3}; 
night_t09_frac_afternoon = night_t09_data{ind_afternoon_july2018,3}; 
night_t14_frac_afternoon = night_t14_data{ind_afternoon_july2018,3}; 

% grab the total night influence for the afternoon hours
night_t01_afternoon = night_t01_data{ind_afternoon_july2018,4}; 
night_t09_afternoon = night_t09_data{ind_afternoon_july2018,4}; 
night_t14_afternoon = night_t14_data{ind_afternoon_july2018,4}; 

% plot the afternoon influence 
figure()
tiledlayout(2,1,'TileSpacing','compact')
nexttile
plot(hours_jul2018(ind_afternoon_july2018),night_t01_afternoon)
title('summed nighttime influence')
hold on 
plot(hours_jul2018(ind_afternoon_july2018),night_t09_afternoon)
plot(hours_jul2018(ind_afternoon_july2018),night_t14_afternoon)
legend('T01','T09','T14')
nexttile
plot(hours_jul2018(ind_afternoon_july2018),night_t01_frac_afternoon)
title('night influence fraction')
hold on 
plot(hours_jul2018(ind_afternoon_july2018),night_t09_frac_afternoon)
plot(hours_jul2018(ind_afternoon_july2018),night_t14_frac_afternoon)

% get the avergae influence from nighttime during the afternoon 
night_t01_frac_afternoon_mean = mean(night_t01_frac_afternoon); 
night_t09_frac_afternoon_mean = mean(night_t09_frac_afternoon); 
night_t14_frac_afternoon_mean = mean(night_t14_frac_afternoon); 

% let's caculate the average value for each hour of the afternoon 
mean(reshape(night_t01_frac_afternoon,6,[]),2,'omitnan')
mean(reshape(night_t09_frac_afternoon,6,[]),2,'omitnan')
mean(reshape(night_t14_frac_afternoon,6,[]),2,'omitnan')
test_hours = reshape(hours_jul2018(ind_afternoon_july2018),6,[]); 
