%%%%%%%% written by: Sam Murphy
%%%%%%%% purpose: compare vprm modeled nee and observed co2 fluxes at ag
%%%%%%%% sites (in INFLUX)

clear all 
close all

%% set paths %% 

% set path to flux data files 
data_in_path = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_data/data_with_flags_removed_and_veg_fraction/'; 

%% load in and format ANWa or US-INd - Corn 2017 %%

% load in as table 
ANWa_table = readtable(strcat(data_in_path,'site_NWa_US-INd_2017.csv')); 

% convert the cell dates to datetime format
ANWa_time_stamps_IN = datetime(table2array(ANWa_table(:,1)),'InputFormat','yyyy-MM-dd HH:mm:ss');

% SHIFT TIME STAMP ***
ANWa_time_stamps = ANWa_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear ANWa_time_stamps_IN

% check if there are an -9999 in the data 
%find(ANWa_table{:,3:126}==-9999)
%find(ANWa_table{:,128:205}==-9999)

% quick plot - check that this generally looks resonable 
figure()
plot(ANWa_time_stamps,ANWa_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('ANWa or US-INd - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% subset 2017 data at this site (corn)
ANWa_time_stamps_2017 = ANWa_time_stamps(year(ANWa_time_stamps)==2017); 
ANWa_2017_ind = find(year(ANWa_time_stamps)==2017); 
ANWa_time_stamps(ANWa_2017_ind)

ANWa_co2_flux = ANWa_table.('co2_flux'); 
ANWa_corn_frac = ANWa_table.('FracVeg'); 

ANWa_co2_flux_2017 = ANWa_co2_flux(ANWa_2017_ind); 
ANWa_corn_frac_2017 = ANWa_corn_frac(ANWa_2017_ind); 

% clean up - dont want to grab these variables on accident 
clear ANWa_corn_frac ANWa_co2_flux ANWa_time_stamps ANWa_2017_ind 

find(isnan(ANWa_co2_flux_2017)==false & isnan(ANWa_corn_frac_2017)==true)

%% load in and format ANWa or US-INd - Corn/Soy split in 2018 %%

% load in as table 
ANWa_2018_split_table = readtable(strcat(data_in_path,'site_NWa_US-INd_2018_Corn_Soy_split.csv')); 

% convert the cell dates to datetime format
ANWa_2018_split_time_stamps_IN = datetime(table2array(ANWa_2018_split_table(:,1)),'InputFormat','yyyy-MM-dd HH:mm:ss');

% SHIFT TIME STAMP ***
ANWa_2018_split_time_stamps = ANWa_2018_split_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear ANWa_2018_split_time_stamps_IN

% check if there are an -9999 in the data 
%find(ANWa_2018_split_table{:,3:126}==-9999)
%find(ANWa_2018_split_table{:,128:205}==-9999)

% quick plot - check that this generally looks resonable - not yet split
% into years (so this is the same as the plot above)
figure()
plot(ANWa_2018_split_time_stamps,ANWa_2018_split_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('ANWa or US-INd - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% grab 2018 only and make vectors for the corn and soy side 
ANWa_time_stamps_2018 = ANWa_2018_split_time_stamps(year(ANWa_2018_split_time_stamps)==2018); 
ANWa_2018_ind = find(year(ANWa_2018_split_time_stamps)==2018); 
ANWa_2018_split_time_stamps(ANWa_2018_ind)

ANWa_co2_flux = ANWa_2018_split_table.('co2_flux'); 
ANWa_corn_frac = ANWa_2018_split_table.('FracCorn'); 
ANWa_soy_frac = ANWa_2018_split_table.('FracSoy'); 

ANWa_co2_flux_2018 = ANWa_co2_flux(ANWa_2018_ind); 
ANWa_corn_frac_2018 = ANWa_corn_frac(ANWa_2018_ind); 
ANWa_soy_frac_2018 = ANWa_soy_frac(ANWa_2018_ind); 

% clean up - dont want to grab these variables on accident 
clear ANWa_corn_frac ANWa_co2_flux ANWa_soy_frac ANWa_2018_split_time_stamps ANWa_2018_ind

find(isnan(ANWa_co2_flux_2018)==false & isnan(ANWa_corn_frac_2018)==true)

%% load in and format ANWb or US-INe  %%
% soybean in 2017, corn in 2018, soy in 2019, corn in 2020 

% load in as table 
ANWb_table = readtable(strcat(data_in_path,'site_NWb_US-INe.csv')); 

% convert the cell dates to datetime format
ANWb_time_stamps_IN = datetime(table2array(ANWb_table(:,1)),'InputFormat','MM/dd/yyyy HH:mm');

% SHIFT TIME STAMP ***
ANWb_time_stamps = ANWb_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear ANWb_time_stamps_IN

% quick plot - check that this generally looks resonable 
figure()
plot(ANWb_time_stamps,ANWb_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('ANWb or US-INe - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% grab the co2 fluxes and the veg fraction for all times 
ANWb_co2_flux = ANWb_table.('co2_flux'); 
ANWb_veg_frac = ANWb_table.('FracVeg'); 

% now seperate out the fluxes by year - 2017 
ANWb_time_stamps_2017 = ANWb_time_stamps(year(ANWb_time_stamps)==2017); 
ANWb_2017_ind = find(year(ANWb_time_stamps)==2017); 
ANWb_time_stamps(ANWb_2017_ind)

ANWb_co2_flux_2017 = ANWb_co2_flux(ANWb_2017_ind); 
ANWb_soy_frac_2017 = ANWb_veg_frac(ANWb_2017_ind); 

% now seperate out the fluxes by year - 2018
ANWb_time_stamps_2018 = ANWb_time_stamps(year(ANWb_time_stamps)==2018); 
ANWb_2018_ind = find(year(ANWb_time_stamps)==2018); 
ANWb_time_stamps(ANWb_2018_ind)

ANWb_co2_flux_2018 = ANWb_co2_flux(ANWb_2018_ind); 
ANWb_corn_frac_2018 = ANWb_veg_frac(ANWb_2018_ind); 

% now seperate out the fluxes by year - 2019 
ANWb_time_stamps_2019 = ANWb_time_stamps(year(ANWb_time_stamps)==2019); 
ANWb_2019_ind = find(year(ANWb_time_stamps)==2019); 
ANWb_time_stamps(ANWb_2019_ind)

ANWb_co2_flux_2019 = ANWb_co2_flux(ANWb_2019_ind); 
ANWb_soy_frac_2019 = ANWb_veg_frac(ANWb_2019_ind);

% now seperate out the fluxes by year - 2020
ANWb_time_stamps_2020 = ANWb_time_stamps(year(ANWb_time_stamps)==2020); 
ANWb_2020_ind = find(year(ANWb_time_stamps)==2020); 
ANWb_time_stamps(ANWb_2020_ind)

ANWb_co2_flux_2020 = ANWb_co2_flux(ANWb_2020_ind); 
ANWb_corn_frac_2020 = ANWb_veg_frac(ANWb_2020_ind);

find(isnan(ANWb_co2_flux_2020)==false & isnan(ANWb_corn_frac_2020)==true)

% clean up 
clear ANWb_co2_flux ANWb_veg_frac ANWb_time_stamps

%% load in and format A09a or US-INi  %%
% soybean in 2019, corn in 2021

% load in as table 
A09a_table = readtable(strcat(data_in_path,'site_09a_US-INi.csv')); 

% convert the cell dates to datetime format
A09a_time_stamps_IN = datetime(table2array(A09a_table(:,1)),'InputFormat','yyyy-MM-dd HH:mm:ss');

% SHIFT TIME STAMP ***
A09a_time_stamps = A09a_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear A09a_time_stamps_IN

% quick plot - check that this generally looks resonable 
figure()
plot(A09a_time_stamps,A09a_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('A09a or US-INi - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% grab the co2 fluxes and the veg fraction for all times 
A09a_co2_flux = A09a_table.('co2_flux'); 
A09a_veg_frac = A09a_table.('FracVeg'); 

% now seperate out the fluxes by year - 2019 
A09a_time_stamps_2019 = A09a_time_stamps(year(A09a_time_stamps)==2019); 
A09a_2019_ind = find(year(A09a_time_stamps)==2019); 
A09a_time_stamps(A09a_2019_ind)

A09a_co2_flux_2019 = A09a_co2_flux(A09a_2019_ind); 
A09a_soy_frac_2019 = A09a_veg_frac(A09a_2019_ind);

% now seperate out the fluxes by year - 2021
A09a_time_stamps_2021 = A09a_time_stamps(year(A09a_time_stamps)==2021); 
A09a_2021_ind = find(year(A09a_time_stamps)==2021); 
A09a_time_stamps(A09a_2021_ind)

A09a_co2_flux_2021 = A09a_co2_flux(A09a_2021_ind); 
A09a_corn_frac_2021 = A09a_veg_frac(A09a_2021_ind);

%% load in and format A09b or US-INj  %%
% corn in 2020, also corn (i think) in 2022 but im not gonna use that year
% for this comparison since i dont have vprm for those times 

% load in as table 
A09b_table = readtable(strcat(data_in_path,'site_09b_US-INj.csv')); 

% convert the cell dates to datetime format
A09b_time_stamps_IN = datetime(table2array(A09b_table(:,1)),'InputFormat','yyyy-MM-dd HH:mm:ss');

% SHIFT TIME STAMP ***
A09b_time_stamps = A09b_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear A09a_time_stamps_IN

% quick plot - check that this generally looks resonable 
figure()
plot(A09b_time_stamps,A09b_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('A09b or US-INj - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% grab the co2 fluxes and the veg fraction for all times 
A09b_co2_flux = A09b_table.('co2_flux'); 
A09b_veg_frac = A09b_table.('FracVeg'); 

% now seperate out the fluxes by year - 2020 
A09b_time_stamps_2020 = A09b_time_stamps(year(A09b_time_stamps)==2020); 
A09b_2020_ind = find(year(A09b_time_stamps)==2020); 
A09b_time_stamps(A09b_2020_ind)

A09b_co2_flux_2020 = A09b_co2_flux(A09b_2020_ind); 
A09b_corn_frac_2020 = A09b_veg_frac(A09b_2020_ind);

%% load in and format A14a or US-INn %%
% corn in 2019 and 2021 

% load in as table 
A14a_table = readtable(strcat(data_in_path,'site_14a_US-INn.csv')); 

% convert the cell dates to datetime format
A14a_time_stamps_IN = datetime(table2array(A14a_table(:,1)),'InputFormat','yyyy-MM-dd HH:mm:ss');

% SHIFT TIME STAMP ***
A14a_time_stamps = A14a_time_stamps_IN - minutes(30); 

% clear original time stamp 
clear A09a_time_stamps_IN

% quick plot - check that this generally looks resonable 
figure()
plot(A14a_time_stamps,A14a_table.('co2_flux'),'.','MarkerEdgeColor','black'); 
title('A14a or US-INn - QC flags removed')
ylabel('co2 flux (umol/m2s)')
grid on 

% grab the co2 fluxes and the veg fraction for all times 
A14a_co2_flux = A14a_table.('co2_flux'); 
A14a_veg_frac = A14a_table.('FracVeg'); 

% now seperate out the fluxes by year - 2019 
A14a_time_stamps_2019 = A14a_time_stamps(year(A14a_time_stamps)==2019); 
A14a_2019_ind = find(year(A14a_time_stamps)==2019); 
A14a_time_stamps(A14a_2019_ind)

A14a_co2_flux_2019 = A14a_co2_flux(A14a_2019_ind); 
A14a_corn_frac_2019 = A14a_veg_frac(A14a_2019_ind);

% now seperate out the fluxes by year - 2021
A14a_time_stamps_2021 = A14a_time_stamps(year(A14a_time_stamps)==2021); 
A14a_2021_ind = find(year(A14a_time_stamps)==2021); 
A14a_time_stamps(A14a_2021_ind)

A14a_co2_flux_2021 = A14a_co2_flux(A14a_2021_ind); 
A14a_corn_frac_2021 = A14a_veg_frac(A14a_2021_ind);

%% double check that u* filtering worked %% 

length(find(A14a_table.("u_star_filt_flag")==1 & isnan(A14a_table.('co2_flux'))==false))
A14a_table.('co2_flux')
A14a_co2_flux(A14a_table.("u_star_filt_flag")==1)

%ANWb_co2_flux(find(ANWb_table.("u_star_filt_flag")==1))
length(find(ANWb_table.("u_star_filt_flag")==1 & isnan(ANWb_table.('co2_flux'))==false))

%% filter the 30 min fluxes by vegfrac %% 
% remove times where the ag fraction is less than 90% 

%%%%%%%%%%%%% ANWa or US-INd for 2017 (corn) %%%%%%%%%%%%% 
ANWa_co2_flux_2017_UNFILT = ANWa_co2_flux_2017; 
ANWa_co2_flux_2017(ANWa_corn_frac_2017<0.9) = NaN; 
ANWa_co2_flux_2017(isnan(ANWa_corn_frac_2017)==true) = NaN; 

% quick plot 
figure()
plot(ANWa_time_stamps_2017,ANWa_co2_flux_2017_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWa_time_stamps_2017,ANWa_co2_flux_2017,'.','MarkerEdgeColor','black')
grid on 

% see how many points were removed 
sum(~isnan(ANWa_co2_flux_2017))
sum(~isnan(ANWa_co2_flux_2017_UNFILT))

%%%%%%%%%%%%% ANWa or US-INd for 2018 (corn/soy split) %%%%%%%%%%%%% 
ANWa_co2_flux_2018_corn_UNFILT = ANWa_co2_flux_2018;
ANWa_co2_flux_2018_soy_UNFILT = ANWa_co2_flux_2018;
ANWa_co2_flux_2018_corn = ANWa_co2_flux_2018; 
ANWa_co2_flux_2018_soy = ANWa_co2_flux_2018; 

ANWa_co2_flux_2018_corn(ANWa_corn_frac_2018<0.9) = NaN; 
ANWa_co2_flux_2018_soy(ANWa_soy_frac_2018<0.9) = NaN; 
ANWa_co2_flux_2018_corn(isnan(ANWa_corn_frac_2018)==true) = NaN; 
ANWa_co2_flux_2018_soy(isnan(ANWa_soy_frac_2018)==true) = NaN; 

figure()
plot(ANWa_time_stamps_2018,ANWa_co2_flux_2018_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWa_time_stamps_2018,ANWa_co2_flux_2018_corn,'.','MarkerEdgeColor','black')

figure()
plot(ANWa_time_stamps_2018,ANWa_co2_flux_2018_soy_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWa_time_stamps_2018,ANWa_co2_flux_2018_soy,'.','MarkerEdgeColor','black')

sum(~isnan(ANWa_co2_flux_2018_corn)) % not enough points for an analysis with corn 
sum(~isnan(ANWa_co2_flux_2018_soy))

%%%%%%%%%%%%% ANWb or US-INe for all times (one map for all) %%%%%%%%%%%%%
ANWb_co2_flux_2020_corn_UNFILT = ANWb_co2_flux_2020; 
ANWb_co2_flux_2020_corn = ANWb_co2_flux_2020; 
ANWb_co2_flux_2020_corn(ANWb_corn_frac_2020<0.9) = NaN; 
ANWb_co2_flux_2020_corn(isnan(ANWb_corn_frac_2020)==true) = NaN; 
sum(~isnan(ANWb_co2_flux_2020_corn))
sum(~isnan(ANWb_co2_flux_2020_corn_UNFILT))

figure()
plot(ANWb_time_stamps_2020,ANWb_co2_flux_2020_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWb_time_stamps_2020,ANWb_co2_flux_2020_corn,'.','MarkerEdgeColor','black')

ANWb_co2_flux_2019_soy_UNFILT = ANWb_co2_flux_2019; 
ANWb_co2_flux_2019_soy = ANWb_co2_flux_2019; 
ANWb_co2_flux_2019_soy(ANWb_soy_frac_2019<0.9) = NaN; 
ANWb_co2_flux_2019_soy(isnan(ANWb_soy_frac_2019)==true) = NaN; 
sum(~isnan(ANWb_co2_flux_2019_soy))
sum(~isnan(ANWb_co2_flux_2019_soy_UNFILT))

figure()
plot(ANWb_time_stamps_2019,ANWb_co2_flux_2019_soy_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWb_time_stamps_2019,ANWb_co2_flux_2019_soy,'.','MarkerEdgeColor','black')

ANWb_co2_flux_2018_corn_UNFILT = ANWb_co2_flux_2018; 
ANWb_co2_flux_2018_corn = ANWb_co2_flux_2018; 
ANWb_co2_flux_2018_corn(ANWb_corn_frac_2018<0.9) = NaN; 
ANWb_co2_flux_2018_corn(isnan(ANWb_corn_frac_2018)==true) = NaN; 
sum(~isnan(ANWb_co2_flux_2018_corn))
sum(~isnan(ANWb_co2_flux_2018_corn_UNFILT))

figure()
plot(ANWb_time_stamps_2018,ANWb_co2_flux_2018_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWb_time_stamps_2018,ANWb_co2_flux_2018_corn,'.','MarkerEdgeColor','black')

ANWb_co2_flux_2017_soy_UNFILT = ANWb_co2_flux_2017; 
ANWb_co2_flux_2017_soy = ANWb_co2_flux_2017; 
ANWb_co2_flux_2017_soy(ANWb_soy_frac_2017<0.9) = NaN; 
ANWb_co2_flux_2017_soy(isnan(ANWb_soy_frac_2017)==true) = NaN;
sum(~isnan(ANWb_co2_flux_2017_soy))
sum(~isnan(ANWb_co2_flux_2017_soy_UNFILT))

figure()
plot(ANWb_time_stamps_2017,ANWb_co2_flux_2017_soy_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(ANWb_time_stamps_2017,ANWb_co2_flux_2017_soy,'.','MarkerEdgeColor','black')


%%%%%%%%%%%%% A09a or US_INi %%%%%%%%%%%%%
A09a_co2_flux_2019_soy_UNFILT = A09a_co2_flux_2019; 
A09a_co2_flux_2019_soy = A09a_co2_flux_2019; 
A09a_co2_flux_2019_soy(A09a_soy_frac_2019<0.9) = NaN; 
A09a_co2_flux_2019_soy(isnan(A09a_soy_frac_2019)==true) = NaN; 
sum(~isnan(A09a_co2_flux_2019_soy))
sum(~isnan(A09a_co2_flux_2019_soy_UNFILT))

figure()
plot(A09a_time_stamps_2019,A09a_co2_flux_2019_soy_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(A09a_time_stamps_2019,A09a_co2_flux_2019_soy,'.','MarkerEdgeColor','black')

A09a_co2_flux_2021_corn_UNFILT = A09a_co2_flux_2021; 
A09a_co2_flux_2021_corn = A09a_co2_flux_2021; 
A09a_co2_flux_2021_corn(A09a_corn_frac_2021<0.9) = NaN; 
A09a_co2_flux_2021_corn(isnan(A09a_corn_frac_2021)==true) = NaN; 
sum(~isnan(A09a_co2_flux_2021_corn))
sum(~isnan(A09a_co2_flux_2021_corn_UNFILT))

figure()
plot(A09a_time_stamps_2021,A09a_co2_flux_2021_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(A09a_time_stamps_2021,A09a_co2_flux_2021_corn,'.','MarkerEdgeColor','black')

%%%%%%%%%%%%% A09b or US-INj %%%%%%%%%%%%%
A09b_co2_flux_2020_corn_UNFILT = A09b_co2_flux_2020; 
A09b_co2_flux_2020_corn = A09b_co2_flux_2020; 
A09b_co2_flux_2020_corn(A09b_corn_frac_2020<0.9) = NaN; 
A09b_co2_flux_2020_corn(isnan(A09b_corn_frac_2020)==true) = NaN; 
sum(~isnan(A09b_co2_flux_2020_corn))
sum(~isnan(A09b_co2_flux_2020_corn_UNFILT))

figure()
plot(A09b_time_stamps_2020,A09b_co2_flux_2020_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(A09b_time_stamps_2020,A09b_co2_flux_2020_corn,'.','MarkerEdgeColor','black')

%%%%%%%%%%%%% A14a or US-INn %%%%%%%%%%%%%
A14a_co2_flux_2021_corn_UNFILT = A14a_co2_flux_2021; 
A14a_co2_flux_2021_corn = A14a_co2_flux_2021; 
A14a_co2_flux_2021_corn(A14a_corn_frac_2021<0.9) = NaN; 
A14a_co2_flux_2021_corn(isnan(A14a_corn_frac_2021)==true) = NaN; 
sum(~isnan(A14a_co2_flux_2021_corn))
sum(~isnan(A14a_co2_flux_2021_corn_UNFILT))

figure()
plot(A14a_time_stamps_2021,A14a_co2_flux_2021_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(A14a_time_stamps_2021,A14a_co2_flux_2021_corn,'.','MarkerEdgeColor','black')

A14a_co2_flux_2019_corn_UNFILT = A14a_co2_flux_2019; 
A14a_co2_flux_2019_corn = A14a_co2_flux_2019; 
A14a_co2_flux_2019_corn(A14a_corn_frac_2019<0.9) = NaN; 
A14a_co2_flux_2019_corn(isnan(A14a_corn_frac_2019)==true) = NaN; 
sum(~isnan(A14a_co2_flux_2019_corn))
sum(~isnan(A14a_co2_flux_2019_corn_UNFILT))

figure()
plot(A14a_time_stamps_2019,A14a_co2_flux_2019_corn_UNFILT,'.','MarkerEdgeColor','blue')
hold on 
plot(A14a_time_stamps_2019,A14a_co2_flux_2019_corn,'.','MarkerEdgeColor','black')


%% calculate hourly average of half hour flux data %% 

%%%%%%%%%%%%% ANWa or US-INd - 2017 %%%%%%%%%%%%%
% make a time table with the flux data (in order to do average)
ANWa_2017_tt = timetable(ANWa_time_stamps_2017,ANWa_co2_flux_2017); 

% the mean in this function omits nans (according to help page)
% i checked the first few half hours manually and this seems to work :)
ANWa_2017_hrly = retime(ANWa_2017_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(ANWa_2017_hrly.ANWa_time_stamps_2017,ANWa_2017_hrly.ANWa_co2_flux_2017,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2017 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% ANWa or US-INd - 2018 - corn and soy %%%%%%%%%%%%%
% make a time table with the flux data (in order to do average)
ANWa_2018_corn_tt = timetable(ANWa_time_stamps_2018,ANWa_co2_flux_2018_corn); 

% calculate the mean for each hour with each half hour 
ANWa_2018_corn_hrly = retime(ANWa_2018_corn_tt ,'hourly','mean'); 

% make a time table with the flux data (in order to do average)
ANWa_2018_soy_tt = timetable(ANWa_time_stamps_2018,ANWa_co2_flux_2018_soy); 

% calculate the mean for each hour with each half hour 
ANWa_2018_soy_hrly = retime(ANWa_2018_soy_tt,'hourly','mean'); 

% quick plots of the data 
figure()
plot(ANWa_2018_corn_hrly.ANWa_time_stamps_2018,ANWa_2018_corn_hrly.ANWa_co2_flux_2018_corn,'.')
title('ANWa or US-INd - 2018 - Corn Side')
ylabel('co2 flux (umol/m2s)')
grid on 

figure()
plot(ANWa_2018_soy_hrly.ANWa_time_stamps_2018,ANWa_2018_soy_hrly.ANWa_co2_flux_2018_soy,'.')
title('ANWa or US-INd - 2018 - Soy Side')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn, 2017 soy 

% 2017 - make a time table with the flux data (in order to do average)
ANWb_2017_soy_tt = timetable(ANWb_time_stamps_2017,ANWb_co2_flux_2017_soy); 

% calculate the mean for each hour with each half hour 
ANWb_2017_soy_hrly = retime(ANWb_2017_soy_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(ANWb_2017_soy_hrly.ANWb_time_stamps_2017,ANWb_2017_soy_hrly.ANWb_co2_flux_2017_soy,'.')
title('ANWb or US-INe - 2017')
ylabel('co2 flux (umol/m2s)')
grid on

% 2018 - make a time table with the flux data (in order to do average)
ANWb_2018_corn_tt = timetable(ANWb_time_stamps_2018,ANWb_co2_flux_2018_corn); 

% calculate the mean for each hour with each half hour 
ANWb_2018_corn_hrly = retime(ANWb_2018_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(ANWb_2018_corn_hrly.ANWb_time_stamps_2018, ANWb_2018_corn_hrly.ANWb_co2_flux_2018_corn,'.')
title('ANWb or US-INe - 2018')
ylabel('co2 flux (umol/m2s)')
grid on

% 2019 - make a time table with the flux data (in order to do average) 
ANWb_2019_soy_tt = timetable(ANWb_time_stamps_2019,ANWb_co2_flux_2019_soy); 

% calculate the mean for each hour with each half hour 
ANWb_2019_soy_hrly = retime(ANWb_2019_soy_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(ANWb_2019_soy_hrly.ANWb_time_stamps_2019,ANWb_2019_soy_hrly.ANWb_co2_flux_2019_soy,'.')
title('ANWb or US-INe - 2019')
ylabel('co2 flux (umol/m2s)')
grid on

% 2020 - make a time table with the flux data (in order to do average)  
ANWb_2020_corn_tt = timetable(ANWb_time_stamps_2020,ANWb_co2_flux_2020_corn); 

% calculate the mean for each hour with each half hour 
ANWb_2020_corn_hrly = retime(ANWb_2020_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(ANWb_2020_corn_hrly.ANWb_time_stamps_2020,ANWb_2020_corn_hrly.ANWb_co2_flux_2020_corn,'.')
title('ANWb or US-INe - 2020')
ylabel('co2 flux (umol/m2s)')
grid on

%%%%%%%%%%%%% A09a or US_INi %%%%%%%%%%%%%
% 2019 soy, 2021 corn 

% 2019 - make a time table with the flux data 
A09a_2019_soy_tt = timetable(A09a_time_stamps_2019,A09a_co2_flux_2019_soy); 

% calculate the mean for each hour with each half hour 
A09a_2019_soy_hrly = retime(A09a_2019_soy_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(A09a_2019_soy_hrly.A09a_time_stamps_2019,A09a_2019_soy_hrly.A09a_co2_flux_2019_soy,'.')
title('A09a or US-INi - 2019')
ylabel('co2 flux (umol/m2s)')
grid on

% 2021 - make a time table with the flux data 
A09a_2021_corn_tt = timetable(A09a_time_stamps_2021,A09a_co2_flux_2021_corn); 

% calculate the mean for each hour with each half hour 
A09a_2021_corn_hrly = retime(A09a_2021_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(A09a_2021_corn_hrly.A09a_time_stamps_2021,A09a_2021_corn_hrly.A09a_co2_flux_2021_corn,'.')
title('A09a or US-INi - 2021')
ylabel('co2 flux (umol/m2s)')
grid on

%%%%%%%%%%%%% A09b or US-INj %%%%%%%%%%%%%
% 2020 corn

% 2020 - make a time table with the flux data 
A09b_2020_corn_tt = timetable(A09b_time_stamps_2020,A09b_co2_flux_2020_corn); 

% calculate the mean for each hour with each half hour 
A09b_2020_corn_hrly = retime(A09b_2020_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(A09b_2020_corn_hrly.A09b_time_stamps_2020,A09b_2020_corn_hrly.A09b_co2_flux_2020_corn,'.')
title('A09b or US-INj - 2020')
ylabel('co2 flux (umol/m2s)')
grid on

%%%%%%%%%%%%% A14a or US-INn %%%%%%%%%%%%%
% 2021 corn, 2019 corn 

% 2021 - make a time table with the flux data 
A14a_2021_corn_tt = timetable(A14a_time_stamps_2021,A14a_co2_flux_2021_corn); 

% calculate the mean for each hour with each half hour 
A14a_2021_corn_hrly = retime(A14a_2021_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(A14a_2021_corn_hrly.A14a_time_stamps_2021,A14a_2021_corn_hrly.A14a_co2_flux_2021_corn,'.')
title('A14a or US-INn - 2021')
ylabel('co2 flux (umol/m2s)')
grid on

% 2019 - make a time table with the flux data 
A14a_2019_corn_tt = timetable(A14a_time_stamps_2019,A14a_co2_flux_2019_corn); 

% calculate the mean for each hour with each half hour 
A14a_2019_corn_hrly = retime(A14a_2019_corn_tt,'hourly','mean'); 

% quick plot of the data 
figure()
plot(A14a_2019_corn_hrly.A14a_time_stamps_2019,A14a_2019_corn_hrly.A14a_co2_flux_2019_corn,'.')
title('A14a or US-INn - 2019')
ylabel('co2 flux (umol/m2s)')
grid on

%% select only the month of JULY for each year of observations %%

%%%%%%%%%%%%%  ANWa 2017 - JULY selection %%%%%%%%%%%%% 
ANWa_JUL_2017_ind = find(month(ANWa_2017_hrly.ANWa_time_stamps_2017)==7); 
ANWa_JUL_2017 = ANWa_2017_hrly.ANWa_time_stamps_2017(ANWa_JUL_2017_ind); 
ANWa_co2_flux_JUL_2017 = ANWa_2017_hrly.ANWa_co2_flux_2017(ANWa_JUL_2017_ind);

figure()
plot(ANWa_JUL_2017,ANWa_co2_flux_JUL_2017,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2017 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWa_co2_flux_JUL_2017))

%%%%%%%%%%%%% ANWa 2018 - JULY selection (both corn and soy) %%%%%%%%%%%%% 
% soy 
ANWa_JUL_2018_soy_ind = find(month(ANWa_2018_soy_hrly.ANWa_time_stamps_2018)==7); 
ANWa_JUL_2018_soy = ANWa_2018_soy_hrly.ANWa_time_stamps_2018(ANWa_JUL_2018_soy_ind); 
ANWa_co2_flux_JUL_2018_soy =ANWa_2018_soy_hrly.ANWa_co2_flux_2018_soy(ANWa_JUL_2018_soy_ind);

% corn
ANWa_JUL_2018_corn_ind = find(month(ANWa_2018_corn_hrly.ANWa_time_stamps_2018)==7); 
ANWa_JUL_2018_corn = ANWa_2018_corn_hrly.ANWa_time_stamps_2018(ANWa_JUL_2018_corn_ind); 
ANWa_co2_flux_JUL_2018_corn =ANWa_2018_corn_hrly.ANWa_co2_flux_2018_corn(ANWa_JUL_2018_corn_ind);

figure()
plot(ANWa_JUL_2018_soy,ANWa_co2_flux_JUL_2018_soy,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2018 - Soy Side')
ylabel('co2 flux (umol/m2s)')
grid on 

figure()
plot(ANWa_JUL_2018_corn,ANWa_co2_flux_JUL_2018_corn,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2018 - Corn Side')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWa_co2_flux_JUL_2018_corn))
sum(~isnan(ANWa_co2_flux_JUL_2018_soy))

%%%%%%%%%%%%% ANWb or US-INe - July selection %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn, 2017 soy 

% 2017 soy 
ANWb_JUL_2017_ind = find(month(ANWb_2017_soy_hrly.ANWb_time_stamps_2017)==7); 
ANWb_JUL_2017 = ANWb_2017_soy_hrly.ANWb_time_stamps_2017(ANWb_JUL_2017_ind); 
ANWb_co2_flux_JUL_2017_soy = ANWb_2017_soy_hrly.ANWb_co2_flux_2017_soy(ANWb_JUL_2017_ind);

figure()
plot(ANWb_JUL_2017,ANWb_co2_flux_JUL_2017_soy,'.')
title('ANWb or US-INe - 2017 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2018 corn 
ANWb_JUL_2018_ind = find(month(ANWb_2018_corn_hrly.ANWb_time_stamps_2018)==7); 
ANWb_JUL_2018 = ANWb_2018_corn_hrly.ANWb_time_stamps_2018(ANWb_JUL_2018_ind); 
ANWb_co2_flux_JUL_2018_corn =ANWb_2018_corn_hrly.ANWb_co2_flux_2018_corn(ANWb_JUL_2018_ind);

figure()
plot(ANWb_JUL_2018,ANWb_co2_flux_JUL_2018_corn,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2018 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_JUL_2018_corn))

% 2019 soy 
ANWb_JUL_2019_ind = find(month(ANWb_2019_soy_hrly.ANWb_time_stamps_2019)==7); 
ANWb_JUL_2019 = ANWb_2019_soy_hrly.ANWb_time_stamps_2019(ANWb_JUL_2019_ind); 
ANWb_co2_flux_JUL_2019_soy =ANWb_2019_soy_hrly.ANWb_co2_flux_2019_soy(ANWb_JUL_2019_ind);

figure()
plot(ANWb_JUL_2019,ANWb_co2_flux_JUL_2019_soy,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_JUL_2019_soy))

% 2020 corn 
ANWb_JUL_2020_ind = find(month(ANWb_2020_corn_hrly.ANWb_time_stamps_2020)==7); 
ANWb_JUL_2020 = ANWb_2020_corn_hrly.ANWb_time_stamps_2020(ANWb_JUL_2020_ind); 
ANWb_co2_flux_JUL_2020_corn =ANWb_2020_corn_hrly.ANWb_co2_flux_2020_corn(ANWb_JUL_2020_ind);

figure()
plot(ANWb_JUL_2020,ANWb_co2_flux_JUL_2020_corn,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_JUL_2020_corn))

%%%%%%%%%%%%% A09a or US-INi %%%%%%%%%%%%%
% 2019 soy, 2021 corn 

% 2019 soy 
A09a_JUL_2019_ind = find(month(A09a_2019_soy_hrly.A09a_time_stamps_2019)==7); 
A09a_JUL_2019 = A09a_2019_soy_hrly.A09a_time_stamps_2019(A09a_JUL_2019_ind); 
A09a_co2_flux_JUL_2019_soy = A09a_2019_soy_hrly.A09a_co2_flux_2019_soy(A09a_JUL_2019_ind);

figure()
plot(A09a_JUL_2019,A09a_co2_flux_JUL_2019_soy,'.','MarkerEdgeColor','black')
title('A09a or US-INi - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09a_co2_flux_JUL_2019_soy))

% 2021 corn 
A09a_JUL_2021_ind = find(month(A09a_2021_corn_hrly.A09a_time_stamps_2021)==7); 
A09a_JUL_2021 = A09a_2021_corn_hrly.A09a_time_stamps_2021(A09a_JUL_2021_ind); 
A09a_co2_flux_JUL_2021_corn =A09a_2021_corn_hrly.A09a_co2_flux_2021_corn(A09a_JUL_2021_ind);

figure()
plot(A09a_JUL_2021,A09a_co2_flux_JUL_2021_corn,'.','MarkerEdgeColor','black')
title('A09a or US-INi - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09a_co2_flux_JUL_2021_corn))

%%%%%%%%%%%%%  A09b or US-INj %%%%%%%%%%%%% 

% 2020 corn 
A09b_JUL_2020_ind = find(month(A09b_2020_corn_hrly.A09b_time_stamps_2020)==7); 
A09b_JUL_2020 = A09b_2020_corn_hrly.A09b_time_stamps_2020(A09b_JUL_2020_ind); 
A09b_co2_flux_JUL_2020_corn =A09b_2020_corn_hrly.A09b_co2_flux_2020_corn(A09b_JUL_2020_ind);

figure()
plot(A09b_JUL_2020,A09b_co2_flux_JUL_2020_corn,'.','MarkerEdgeColor','black')
title('A09b or US-INj - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09b_co2_flux_JUL_2020_corn))


%%%%%%%%%%%%%  A14a or US-INn - 2021 corn, 2019 corn %%%%%%%%%%%%% 

% 2021 corn 
A14a_JUL_2021_ind = find(month(A14a_2021_corn_hrly.A14a_time_stamps_2021)==7); 
A14a_JUL_2021 = A14a_2021_corn_hrly.A14a_time_stamps_2021(A14a_JUL_2021_ind); 
A14a_co2_flux_JUL_2021_corn =A14a_2021_corn_hrly.A14a_co2_flux_2021_corn(A14a_JUL_2021_ind);

figure()
plot(A14a_JUL_2021,A14a_co2_flux_JUL_2021_corn,'.','MarkerEdgeColor','black')
title('A14a or US-INn - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A14a_co2_flux_JUL_2021_corn))

% 2019 corn 
A14a_JUL_2019_ind = find(month(A14a_2019_corn_hrly.A14a_time_stamps_2019)==7); 
A14a_JUL_2019 = A14a_2019_corn_hrly.A14a_time_stamps_2019(A14a_JUL_2019_ind); 
A14a_co2_flux_JUL_2019_corn =A14a_2019_corn_hrly.A14a_co2_flux_2019_corn(A14a_JUL_2019_ind);

figure()
plot(A14a_JUL_2019,A14a_co2_flux_JUL_2019_corn,'.','MarkerEdgeColor','black')
title('A14a or US-INn - 2019 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A14a_co2_flux_JUL_2019_corn))

%% select only the month of AUGUST for each year of observations %%

%%%%%%%%%%%%% ANWa 2017 %%%%%%%%%%%%%
ANWa_AUG_2017_ind = find(month(ANWa_2017_hrly.ANWa_time_stamps_2017)==8); 
ANWa_AUG_2017 = ANWa_2017_hrly.ANWa_time_stamps_2017(ANWa_AUG_2017_ind); 
ANWa_co2_flux_AUG_2017 = ANWa_2017_hrly.ANWa_co2_flux_2017(ANWa_AUG_2017_ind);

figure()
plot(ANWa_AUG_2017,ANWa_co2_flux_AUG_2017,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2017 - Aug - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWa_co2_flux_AUG_2017))

%%%%%%%%%%%%% ANWa 2018 - AUG selection (both corn and soy) %%%%%%%%%%%%%
ANWa_AUG_2018_soy_ind = find(month(ANWa_2018_soy_hrly.ANWa_time_stamps_2018)==8); 
ANWa_AUG_2018_soy = ANWa_2018_soy_hrly.ANWa_time_stamps_2018(ANWa_AUG_2018_soy_ind); 
ANWa_co2_flux_AUG_2018_soy =ANWa_2018_soy_hrly.ANWa_co2_flux_2018_soy(ANWa_AUG_2018_soy_ind);

ANWa_AUG_2018_corn_ind = find(month(ANWa_2018_corn_hrly.ANWa_time_stamps_2018)==8); 
ANWa_AUG_2018_corn = ANWa_2018_corn_hrly.ANWa_time_stamps_2018(ANWa_AUG_2018_corn_ind); 
ANWa_co2_flux_AUG_2018_corn =ANWa_2018_corn_hrly.ANWa_co2_flux_2018_corn(ANWa_AUG_2018_corn_ind);

figure()
plot(ANWa_AUG_2018_soy,ANWa_co2_flux_AUG_2018_soy,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2018 - Soy Side')
ylabel('co2 flux (umol/m2s)')
grid on 

figure()
plot(ANWa_AUG_2018_corn,ANWa_co2_flux_AUG_2018_corn,'.','MarkerEdgeColor','black')
title('ANWa or US-INd - 2018 - Corn Side')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWa_co2_flux_AUG_2018_corn))
sum(~isnan(ANWa_co2_flux_AUG_2018_soy))

%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn, 2017 soy 

% 2017 
ANWb_AUG_2017_ind = find(month(ANWb_2017_soy_hrly.ANWb_time_stamps_2017)==8); 
ANWb_AUG_2017 = ANWb_2017_soy_hrly.ANWb_time_stamps_2017(ANWb_AUG_2017_ind); 
ANWb_co2_flux_AUG_2017_soy =ANWb_2017_soy_hrly.ANWb_co2_flux_2017_soy(ANWb_AUG_2017_ind);

figure()
plot(ANWb_AUG_2017,ANWb_co2_flux_AUG_2017_soy,'.')
title('ANWb or US-INe - 2017 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2018
ANWb_AUG_2018_ind = find(month(ANWb_2018_corn_hrly.ANWb_time_stamps_2018)==8); 
ANWb_AUG_2018 = ANWb_2018_corn_hrly.ANWb_time_stamps_2018(ANWb_AUG_2018_ind); 
ANWb_co2_flux_AUG_2018_corn =ANWb_2018_corn_hrly.ANWb_co2_flux_2018_corn(ANWb_AUG_2018_ind);

figure()
plot(ANWb_AUG_2018,ANWb_co2_flux_AUG_2018_corn,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2018 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_AUG_2018_corn))

% 2019
ANWb_AUG_2019_ind = find(month(ANWb_2019_soy_hrly.ANWb_time_stamps_2019)==8); 
ANWb_AUG_2019 = ANWb_2019_soy_hrly.ANWb_time_stamps_2019(ANWb_AUG_2019_ind); 
ANWb_co2_flux_AUG_2019_soy = ANWb_2019_soy_hrly.ANWb_co2_flux_2019_soy(ANWb_AUG_2019_ind);

figure()
plot(ANWb_AUG_2019,ANWb_co2_flux_AUG_2019_soy,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_AUG_2019_soy))

% 2020
ANWb_AUG_2020_ind = find(month(ANWb_2020_corn_hrly.ANWb_time_stamps_2020)==8); 
ANWb_AUG_2020 = ANWb_2020_corn_hrly.ANWb_time_stamps_2020(ANWb_AUG_2020_ind); 
ANWb_co2_flux_AUG_2020_corn = ANWb_2020_corn_hrly.ANWb_co2_flux_2020_corn(ANWb_AUG_2020_ind);

figure()
plot(ANWb_AUG_2020,ANWb_co2_flux_AUG_2020_corn,'.','MarkerEdgeColor','black')
title('ANWb or US-INe - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(ANWb_co2_flux_AUG_2020_corn))

%%%%%%%%%%%%% A09a or US-INi %%%%%%%%%%%%%
% 2019 soy, 2021 corn 

% 2019
A09a_AUG_2019_ind = find(month(A09a_2019_soy_hrly.A09a_time_stamps_2019)==8); 
A09a_AUG_2019 = A09a_2019_soy_hrly.A09a_time_stamps_2019(A09a_AUG_2019_ind); 
A09a_co2_flux_AUG_2019_soy = A09a_2019_soy_hrly.A09a_co2_flux_2019_soy(A09a_AUG_2019_ind);

figure()
plot(A09a_AUG_2019,A09a_co2_flux_AUG_2019_soy,'.','MarkerEdgeColor','black')
title('A09a or US-INi - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09a_co2_flux_AUG_2019_soy))

% 2021 
A09a_AUG_2021_ind = find(month(A09a_2021_corn_hrly.A09a_time_stamps_2021)==8); 
A09a_AUG_2021 = A09a_2021_corn_hrly.A09a_time_stamps_2021(A09a_AUG_2021_ind); 
A09a_co2_flux_AUG_2021_corn =A09a_2021_corn_hrly.A09a_co2_flux_2021_corn(A09a_AUG_2021_ind);

figure()
plot(A09a_AUG_2021,A09a_co2_flux_AUG_2021_corn,'.','MarkerEdgeColor','black')
title('A09a or US-INi - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09a_co2_flux_AUG_2021_corn))

%%%%%%%%%%%%% A09b or US-INj %%%%%%%%%%%%%
% 2020 
A09b_AUG_2020_ind = find(month(A09b_2020_corn_hrly.A09b_time_stamps_2020)==8); 
A09b_AUG_2020 = A09b_2020_corn_hrly.A09b_time_stamps_2020(A09b_AUG_2020_ind); 
A09b_co2_flux_AUG_2020_corn =A09b_2020_corn_hrly.A09b_co2_flux_2020_corn(A09b_AUG_2020_ind);

figure()
plot(A09b_AUG_2020,A09b_co2_flux_AUG_2020_corn,'.','MarkerEdgeColor','black')
title('A09b or US-INj - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A09b_co2_flux_AUG_2020_corn))

%%%%%%%%%%%%% A14a or US-INn %%%%%%%%%%%%%
% 2021 corn, 2019 corn 

% 2021
A14a_AUG_2021_ind = find(month(A14a_2021_corn_hrly.A14a_time_stamps_2021)==8); 
A14a_AUG_2021 = A14a_2021_corn_hrly.A14a_time_stamps_2021(A14a_AUG_2021_ind); 
A14a_co2_flux_AUG_2021_corn = A14a_2021_corn_hrly.A14a_co2_flux_2021_corn(A14a_AUG_2021_ind);

figure()
plot(A14a_AUG_2021,A14a_co2_flux_AUG_2021_corn,'.','MarkerEdgeColor','black')
title('A14a or US-INn - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A14a_co2_flux_AUG_2021_corn))

% 2019
A14a_AUG_2019_ind = find(month(A14a_2019_corn_hrly.A14a_time_stamps_2019)==8); 
A14a_AUG_2019 = A14a_2019_corn_hrly.A14a_time_stamps_2019(A14a_AUG_2019_ind); 
A14a_co2_flux_AUG_2019_corn =A14a_2019_corn_hrly.A14a_co2_flux_2019_corn(A14a_AUG_2019_ind);

figure()
plot(A14a_AUG_2019,A14a_co2_flux_AUG_2019_corn,'.','MarkerEdgeColor','black')
title('A14a or US-INn - 2019 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

sum(~isnan(A14a_co2_flux_AUG_2019_corn))


%% load in VPRM ouputs for each site %% 
% set input path 
vprm_outputs_path = '/storage/group/zrb5027/default/INFLUX/smm8236/flux_comparisons_vprm/vprm_out_point/';

% ANWa - 2018 just soy 
ANWa_soy_2018_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_NWa_Other_Crops_2018.csv'));
ANWa_soy_2018_vprm_dates = datetime(ANWa_soy_2018_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
ANWa_soy_2018_vprm_nee = ANWa_soy_2018_vprm.nee_vec_at_point; 

% ANWb - 2018 corn, 2019 soy, 2020 corn
ANWb_corn_2018_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_NWb_Corn_2018.csv'));
ANWb_corn_2018_vprm_dates = datetime(ANWb_corn_2018_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
ANWb_corn_2018_vprm_nee = ANWb_corn_2018_vprm.nee_vec_at_point; 

ANWb_soy_2019_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_NWb_Other_Crops_2019.csv'));
ANWb_soy_2019_vprm_dates = datetime(ANWb_soy_2019_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
ANWb_soy_2019_vprm_nee = ANWb_soy_2019_vprm.nee_vec_at_point; 

ANWb_corn_2020_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_NWb_Corn_2020.csv'));
ANWb_corn_2020_vprm_dates = datetime(ANWb_corn_2020_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
ANWb_corn_2020_vprm_nee = ANWb_corn_2020_vprm.nee_vec_at_point; 

% A09a - 2019 soy, 2021 corn
A09a_soy_2019_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_09a_Other_Crops_2019.csv'));
A09a_soy_2019_vprm_dates = datetime(A09a_soy_2019_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
A09a_soy_2019_vprm_nee = A09a_soy_2019_vprm.nee_vec_at_point; 

A09a_corn_2021_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_09a_Corn_2021.csv'));
A09a_corn_2021_vprm_dates = datetime(A09a_corn_2021_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
A09a_corn_2021_vprm_nee = A09a_corn_2021_vprm.nee_vec_at_point; 

% A09b - 2020 corn 
A09b_corn_2020_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_09b_Corn_2020.csv'));
A09b_corn_2020_vprm_dates = datetime(A09b_corn_2020_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
A09b_corn_2020_vprm_nee = A09b_corn_2020_vprm.nee_vec_at_point; 

% A14a - 2019 corn, 2021 corn
A14a_corn_2021_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_14a_Corn_2021.csv'));
A14a_corn_2021_vprm_dates = datetime(A14a_corn_2021_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
A14a_corn_2021_vprm_nee = A14a_corn_2021_vprm.nee_vec_at_point; 

A14a_corn_2019_vprm = readtable(strcat(vprm_outputs_path,'vprm_nee_Site_A_14a_Corn_2019.csv'));
A14a_corn_2019_vprm_dates = datetime(A14a_corn_2019_vprm.Var2,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
A14a_corn_2019_vprm_nee = A14a_corn_2019_vprm.nee_vec_at_point; 

%% select only july for each year of VPRM outputs %% 
% ANWa - 2018 just soy 
ind_JUL_ANWa_soy_2018_vprm_dates = find(month(ANWa_soy_2018_vprm_dates)==7); 
ANWa_soy_2018_vprm_nee_JUL = ANWa_soy_2018_vprm_nee(ind_JUL_ANWa_soy_2018_vprm_dates); 

% ANWb - 2018 corn, 2019 soy, 2020 corn
ind_JUL_ANWb_corn_2018_vprm_dates = find(month(ANWb_corn_2018_vprm_dates)==7); 
ANWb_corn_2018_vprm_nee_JUL = ANWb_corn_2018_vprm_nee(ind_JUL_ANWb_corn_2018_vprm_dates); 

ind_JUL_ANWb_soy_2019_vprm_dates = find(month(ANWb_soy_2019_vprm_dates)==7); 
ANWb_soy_2019_vprm_nee_JUL = ANWb_soy_2019_vprm_nee(ind_JUL_ANWb_soy_2019_vprm_dates); 

ind_JUL_ANWb_corn_2020_vprm_dates = find(month(ANWb_corn_2020_vprm_dates)==7); 
ANWb_corn_2020_vprm_nee_JUL = ANWb_corn_2020_vprm_nee(ind_JUL_ANWb_corn_2020_vprm_dates); 

% A09a - 2019 soy, 2021 corn
ind_JUL_A09a_soy_2019_vprm_dates = find(month(A09a_soy_2019_vprm_dates)==7); 
A09a_soy_2019_vprm_nee_JUL = A09a_soy_2019_vprm_nee(ind_JUL_A09a_soy_2019_vprm_dates); 

ind_JUL_A09a_corn_2021_vprm_dates = find(month(A09a_corn_2021_vprm_dates)==7); 
A09a_corn_2021_vprm_nee_JUL = A09a_corn_2021_vprm_nee(ind_JUL_A09a_corn_2021_vprm_dates); 

% A09b - 2020 corn 
ind_JUL_A09b_corn_2020_vprm_dates = find(month(A09b_corn_2020_vprm_dates)==7); 
A09b_corn_2020_vprm_nee_JUL = A09b_corn_2020_vprm_nee(ind_JUL_A09b_corn_2020_vprm_dates); 

% A14a - 2019 corn, 2021 corn
ind_JUL_A14a_corn_2021_vprm_dates = find(month(A14a_corn_2021_vprm_dates)==7); 
A14a_corn_2021_vprm_nee_JUL = A14a_corn_2021_vprm_nee(ind_JUL_A14a_corn_2021_vprm_dates); 

ind_JUL_A14a_corn_2019_vprm_dates = find(month(A14a_corn_2019_vprm_dates)==7); 
A14a_corn_2019_vprm_nee_JUL = A14a_corn_2019_vprm_nee(ind_JUL_A14a_corn_2019_vprm_dates); 


%% grab AUG for each year of vprm outputs 
% ANWa - 2018 just soy 
ind_AUG_ANWa_soy_2018_vprm_dates = find(month(ANWa_soy_2018_vprm_dates)==8); 
ANWa_soy_2018_vprm_nee_AUG = ANWa_soy_2018_vprm_nee(ind_AUG_ANWa_soy_2018_vprm_dates); 

% ANWb - 2018 corn, 2019 soy, 2020 corn
ind_AUG_ANWb_corn_2018_vprm_dates = find(month(ANWb_corn_2018_vprm_dates)==8); 
ANWb_corn_2018_vprm_nee_AUG = ANWb_corn_2018_vprm_nee(ind_AUG_ANWb_corn_2018_vprm_dates); 

ind_AUG_ANWb_soy_2019_vprm_dates = find(month(ANWb_soy_2019_vprm_dates)==8); 
ANWb_soy_2019_vprm_nee_AUG = ANWb_soy_2019_vprm_nee(ind_AUG_ANWb_soy_2019_vprm_dates); 

ind_AUG_ANWb_corn_2020_vprm_dates = find(month(ANWb_corn_2020_vprm_dates)==8); 
ANWb_corn_2020_vprm_nee_AUG = ANWb_corn_2020_vprm_nee(ind_AUG_ANWb_corn_2020_vprm_dates); 

% A09a - 2019 soy, 2021 corn
ind_AUG_A09a_soy_2019_vprm_dates = find(month(A09a_soy_2019_vprm_dates)==8); 
A09a_soy_2019_vprm_nee_AUG = A09a_soy_2019_vprm_nee(ind_AUG_A09a_soy_2019_vprm_dates); 

ind_AUG_A09a_corn_2021_vprm_dates = find(month(A09a_corn_2021_vprm_dates)==8); 
A09a_corn_2021_vprm_nee_AUG = A09a_corn_2021_vprm_nee(ind_AUG_A09a_corn_2021_vprm_dates); 

% A09b - 2020 corn 
ind_AUG_A09b_corn_2020_vprm_dates = find(month(A09b_corn_2020_vprm_dates)==8); 
A09b_corn_2020_vprm_nee_AUG = A09b_corn_2020_vprm_nee(ind_AUG_A09b_corn_2020_vprm_dates); 

% A14a - 2019 corn, 2021 corn
ind_AUG_A14a_corn_2021_vprm_dates = find(month(A14a_corn_2021_vprm_dates)==8); 
A14a_corn_2021_vprm_nee_AUG = A14a_corn_2021_vprm_nee(ind_AUG_A14a_corn_2021_vprm_dates); 

ind_AUG_A14a_corn_2019_vprm_dates = find(month(A14a_corn_2019_vprm_dates)==8); 
A14a_corn_2019_vprm_nee_AUG = A14a_corn_2019_vprm_nee(ind_AUG_A14a_corn_2019_vprm_dates); 


%% quick plots of VPRM NEE in july %% 
figure()
plot(A14a_corn_2019_vprm_dates(ind_JUL_A14a_corn_2019_vprm_dates),A14a_corn_2019_vprm_nee_JUL,'.','MarkerEdgeColor','black')
grid on 
ylabel('VPRM NEE')
xlabel('Time')
title('Site 14a - July 2019')

figure()
plot(A14a_corn_2021_vprm_dates(ind_JUL_A14a_corn_2021_vprm_dates),A14a_corn_2021_vprm_nee_JUL,'.','MarkerEdgeColor','black')
grid on 
ylabel('VPRM NEE')
xlabel('Time')
title('Site 14a - July 2021')

figure()
plot(A09a_soy_2019_vprm_dates(ind_JUL_A09a_soy_2019_vprm_dates),A09a_soy_2019_vprm_nee_JUL,'.','MarkerEdgeColor','black')
grid on 
ylabel('VPRM NEE')
xlabel('Time')
title('Site 09a - July 2019')

figure()
plot(A09a_soy_2019_vprm_dates(ind_JUL_A09a_soy_2019_vprm_dates),A09a_soy_2019_vprm_nee_JUL,'.','MarkerEdgeColor','black')
grid on 
ylabel('VPRM NEE')
xlabel('Time')
title('Site 09a - July 2019')


%% remove vprm data based on flux data availability - July %% 
% remove points from the model data when we are missing observations 

%%%%%%%%%%%%% ANWa %%%%%%%%%%%%%
% 2018 soy 
% check the lengths are equal 
length(ANWa_co2_flux_JUL_2018_soy)
length(ANWa_soy_2018_vprm_nee_JUL)
% remove model points when obs missing 
ANWa_soy_2018_vprm_nee_JUL(isnan(ANWa_co2_flux_JUL_2018_soy)==true) = NaN; 
% check that the number of nans is equal 
sum(~isnan(ANWa_soy_2018_vprm_nee_JUL))
sum(~isnan(ANWa_co2_flux_JUL_2018_soy))

%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%% 
% 2020 corn, 2019 soy, 2018 corn
% check the lengths are equal 
length(ANWb_co2_flux_JUL_2018_corn)
length(ANWb_corn_2018_vprm_nee_JUL)
% remove model points when obs missing 
ANWb_corn_2018_vprm_nee_JUL(isnan(ANWb_co2_flux_JUL_2018_corn)==true) = NaN; 
% check that the number of nans is equal 
sum(~isnan(ANWb_corn_2018_vprm_nee_JUL))
sum(~isnan(ANWb_co2_flux_JUL_2018_corn))

% check the lengths are equal 
length(ANWb_co2_flux_JUL_2019_soy)
length(ANWb_soy_2019_vprm_nee_JUL)
% remove model points when obs missing 
ANWb_soy_2019_vprm_nee_JUL(isnan(ANWb_co2_flux_JUL_2019_soy)==true) = NaN; 
% check that the number of nans is equal 
sum(~isnan(ANWb_soy_2019_vprm_nee_JUL))
sum(~isnan(ANWb_co2_flux_JUL_2019_soy))

% check the lengths are equal 
length(ANWb_co2_flux_JUL_2020_corn)
length(ANWb_corn_2020_vprm_nee_JUL)
% remove model points when obs missing 
ANWb_corn_2020_vprm_nee_JUL(isnan(ANWb_co2_flux_JUL_2020_corn)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(ANWb_corn_2020_vprm_nee_JUL))
sum(~isnan(ANWb_co2_flux_JUL_2020_corn))

%%%%%%%%%%%%% site A09a or US INi %%%%%%%%%%%%%
% check the lengths are equal 
length(A09a_co2_flux_JUL_2019_soy)
length(A09a_soy_2019_vprm_nee_JUL)
% remove model points when obs missing 
A09a_soy_2019_vprm_nee_JUL(isnan(A09a_co2_flux_JUL_2019_soy)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(A09a_soy_2019_vprm_nee_JUL))
sum(~isnan(A09a_co2_flux_JUL_2019_soy))

% check the lengths are equal 
length(A09a_co2_flux_JUL_2021_corn)
length(A09a_corn_2021_vprm_nee_JUL)
% remove model points when obs missing 
A09a_corn_2021_vprm_nee_JUL(isnan(A09a_co2_flux_JUL_2021_corn)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(A09a_corn_2021_vprm_nee_JUL))
sum(~isnan(A09a_co2_flux_JUL_2021_corn))

%%%%%%%%%%%%% A09b 2020 or US-INj %%%%%%%%%%%%%
% check the lengths are equal 
length(A09b_co2_flux_JUL_2020_corn)
length(A09b_corn_2020_vprm_nee_JUL)
% remove model points when obs missing 
A09b_corn_2020_vprm_nee_JUL(isnan(A09b_co2_flux_JUL_2020_corn)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(A09b_corn_2020_vprm_nee_JUL))
sum(~isnan(A09b_co2_flux_JUL_2020_corn))

%%%%%%%%%%%%% site A14a or US INn %%%%%%%%%%%%%
% check the lengths are equal 
length(A14a_co2_flux_JUL_2021_corn)
length(A14a_corn_2021_vprm_nee_JUL)
% remove model points when obs missing 
A14a_corn_2021_vprm_nee_JUL(isnan(A14a_co2_flux_JUL_2021_corn)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(A14a_corn_2021_vprm_nee_JUL))
sum(~isnan(A14a_co2_flux_JUL_2021_corn))

% check the lengths are equal 
length(A14a_co2_flux_JUL_2019_corn)
length(A14a_corn_2019_vprm_nee_JUL)
% remove model points when obs missing 
A14a_corn_2019_vprm_nee_JUL(isnan(A14a_co2_flux_JUL_2019_corn)==true) = NaN; 
% check the number of nans is equal 
sum(~isnan(A14a_corn_2019_vprm_nee_JUL))
sum(~isnan(A14a_co2_flux_JUL_2019_corn))

%% remove vprm data based on flux data availability - Aug %% 
% remove points from the model data when we are missing observations 

%%%%%%%%%%%%% ANWa %%%%%%%%%%%%%
% 2018 soy 
% check the lengths are equal 
length(ANWa_co2_flux_AUG_2018_soy)
length(ANWa_soy_2018_vprm_nee_AUG)
% remove model points when obs missing 
ANWa_soy_2018_vprm_nee_AUG(isnan(ANWa_co2_flux_AUG_2018_soy)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(ANWa_soy_2018_vprm_nee_AUG))
sum(~isnan(ANWa_co2_flux_AUG_2018_soy))

%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn
% check the lengths are equal 
length(ANWb_co2_flux_AUG_2018_corn)
length(ANWb_corn_2018_vprm_nee_AUG)
% remove model points when obs missing 
ANWb_corn_2018_vprm_nee_AUG(isnan(ANWb_co2_flux_AUG_2018_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(ANWb_corn_2018_vprm_nee_AUG))
sum(~isnan(ANWb_co2_flux_AUG_2018_corn))

% check the lengths are equal 
length(ANWb_co2_flux_AUG_2019_soy)
length(ANWb_soy_2019_vprm_nee_AUG)
% remove model points when obs missing 
ANWb_soy_2019_vprm_nee_AUG(isnan(ANWb_co2_flux_AUG_2019_soy)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(ANWb_soy_2019_vprm_nee_AUG))
sum(~isnan(ANWb_co2_flux_AUG_2019_soy))

% check the lengths are equal 
length(ANWb_co2_flux_AUG_2020_corn)
length(ANWb_corn_2020_vprm_nee_AUG)
% remove model points when obs missing 
ANWb_corn_2020_vprm_nee_AUG(isnan(ANWb_co2_flux_AUG_2020_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(ANWb_corn_2020_vprm_nee_AUG))
sum(~isnan(ANWb_co2_flux_AUG_2020_corn))

%%%%%%%%%%%%% site A09a or US INi %%%%%%%%%%%%%
% check the lengths are equal 
length(A09a_co2_flux_AUG_2019_soy)
length(A09a_soy_2019_vprm_nee_AUG)
% remove model points when obs missing 
A09a_soy_2019_vprm_nee_AUG(isnan(A09a_co2_flux_AUG_2019_soy)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(A09a_soy_2019_vprm_nee_AUG))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))

% check the lengths are equal 
length(A09a_co2_flux_AUG_2021_corn)
length(A09a_corn_2021_vprm_nee_AUG)
% remove model points when obs missing 
A09a_corn_2021_vprm_nee_AUG(isnan(A09a_co2_flux_AUG_2021_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(A09a_corn_2021_vprm_nee_AUG))
sum(~isnan(A09a_co2_flux_AUG_2021_corn))

%%%%%%%%%%%%% A09b 2020 or US-INj %%%%%%%%%%%%%
% check the lengths are equal 
length(A09b_co2_flux_AUG_2020_corn)
length(A09b_corn_2020_vprm_nee_AUG)
% remove model points when obs missing 
A09b_corn_2020_vprm_nee_AUG(isnan(A09b_co2_flux_AUG_2020_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(A09b_corn_2020_vprm_nee_AUG))
sum(~isnan(A09b_co2_flux_AUG_2020_corn))

%%%%%%%%%%%%% site A14a or US INn %%%%%%%%%%%%%
% check the lengths are equal 
length(A14a_co2_flux_AUG_2021_corn)
length(A14a_corn_2021_vprm_nee_AUG)
% remove model points when obs missing 
A14a_corn_2021_vprm_nee_AUG(isnan(A14a_co2_flux_AUG_2021_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(A14a_corn_2021_vprm_nee_AUG))
sum(~isnan(A14a_co2_flux_AUG_2021_corn))

% check the lengths are equal 
length(A14a_co2_flux_AUG_2019_corn)
length(A14a_corn_2019_vprm_nee_AUG)
% remove model points when obs missing 
A14a_corn_2019_vprm_nee_AUG(isnan(A14a_co2_flux_AUG_2019_corn)==true) = NaN; 
% check the number of nans is equal
sum(~isnan(A14a_corn_2019_vprm_nee_AUG))
sum(~isnan(A14a_co2_flux_AUG_2019_corn))

%% quick plots of the model and obs for july (for each hourly point) %%

% check the number of nans 
sum(~isnan(ANWa_co2_flux_JUL_2018_soy))
sum(~isnan(ANWb_co2_flux_JUL_2018_corn))
sum(~isnan(ANWb_co2_flux_JUL_2019_soy))
sum(~isnan(ANWb_co2_flux_JUL_2020_corn))
sum(~isnan(A09a_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_JUL_2021_corn))
sum(~isnan(A09b_co2_flux_JUL_2020_corn))
sum(~isnan(A14a_co2_flux_JUL_2021_corn))
sum(~isnan(A14a_co2_flux_JUL_2019_corn))

%%%%%%%%%%%%% ANWa 2018 %%%%%%%%%%%%%
% JUL selection (both corn and soy) 
figure()
plot(ANWa_JUL_2018_soy,ANWa_co2_flux_JUL_2018_soy,'.','MarkerEdgeColor','black')
hold on 
plot(ANWa_JUL_2018_soy,ANWa_soy_2018_vprm_nee_JUL,'.')
title('ANWa or US-INd - 2018 - soy')
ylabel('co2 flux (umol/m2s)')
grid on 
legend('obs','model')

%%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn, 2017 soy 
% 2017 - NO DATA IN JULY 
figure()
plot(ANWb_JUL_2017,ANWb_co2_flux_JUL_2017_soy,'.')

% 2018
figure()
plot(ANWb_JUL_2018,ANWb_co2_flux_JUL_2018_corn,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_JUL_2018,ANWb_corn_2018_vprm_nee_JUL,'.')
legend('obs','model')
title('ANWb or US-INe - 2018 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2019
figure()
plot(ANWb_JUL_2019,ANWb_co2_flux_JUL_2019_soy,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_JUL_2019,ANWb_soy_2019_vprm_nee_JUL,'.')
legend('obs','model')
title('ANWb or US-INe - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2020
figure()
plot(ANWb_JUL_2020,ANWb_co2_flux_JUL_2020_corn,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_JUL_2020,ANWb_corn_2020_vprm_nee_JUL,'.')
legend('obs','model')
title('ANWb or US-INe - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A09a or US-INi %%%%%%%%%%%%%
% 2019 soy, 2021 corn 
% 2019
figure()
plot(A09a_JUL_2019,A09a_co2_flux_JUL_2019_soy,'.','MarkerEdgeColor','black')
hold on 
plot(A09a_JUL_2019,A09a_soy_2019_vprm_nee_JUL,'.')
legend('obs','model')
title('A09a or US-INi - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2021 
figure()
plot(A09a_JUL_2021,A09a_co2_flux_JUL_2021_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A09a_JUL_2021,A09a_corn_2021_vprm_nee_JUL,'.')
title('A09a or US-INi - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A09b or US-INj %%%%%%%%%%%%%
% 2020 
figure()
plot(A09b_JUL_2020,A09b_co2_flux_JUL_2020_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A09b_JUL_2020,A09b_corn_2020_vprm_nee_JUL,'.')
legend('obs','model')
title('A09b or US-INj - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A14a or US-INn %%%%%%%%%%%%%
% 2021 corn, 2019 corn 
% 2021
figure()
plot(A14a_JUL_2021,A14a_co2_flux_JUL_2021_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A14a_JUL_2021,A14a_corn_2021_vprm_nee_JUL,'.')
legend('obs','model')
title('A14a or US-INn - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2019
figure()
plot(A14a_JUL_2019,A14a_co2_flux_JUL_2019_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A14a_JUL_2019,A14a_corn_2019_vprm_nee_JUL,'.')
legend('obs','model')
title('A14a or US-INn - 2019 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%% quick plots of the model and obs for august (for each hourly point) %%
% check the number of nans 
sum(~isnan(ANWa_co2_flux_AUG_2018_soy))
sum(~isnan(ANWb_co2_flux_AUG_2018_corn))
sum(~isnan(ANWb_co2_flux_AUG_2019_soy))
sum(~isnan(ANWb_co2_flux_AUG_2020_corn))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))
sum(~isnan(A09a_co2_flux_AUG_2021_corn))
sum(~isnan(A09b_co2_flux_AUG_2020_corn))
sum(~isnan(A14a_co2_flux_AUG_2021_corn))
sum(~isnan(A14a_co2_flux_AUG_2019_corn))

%%%%%%%%%%%%% ANWa 2018 %%%%%%%%%%%%%
% AUG selection (both corn and soy) 
figure()
plot(ANWa_AUG_2018_soy,ANWa_co2_flux_AUG_2018_soy,'.','MarkerEdgeColor','black')
hold on 
plot(ANWa_AUG_2018_soy,ANWa_soy_2018_vprm_nee_AUG,'.')
title('ANWa or US-INd - 2018 - soy')
ylabel('co2 flux (umol/m2s)')
grid on 
legend('obs','model')

%%%%%%%%%%%%% ANWb or US-INe %%%%%%%%%%%%%
% 2020 corn, 2019 soy, 2018 corn, 2017 soy 
% 2017 - no data in aug 
figure()
plot(ANWb_AUG_2017,ANWb_co2_flux_AUG_2017_soy,'.')

% 2018
figure()
plot(ANWb_AUG_2018,ANWb_co2_flux_AUG_2018_corn,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_AUG_2018,ANWb_corn_2018_vprm_nee_AUG,'.')
legend('obs','model')
title('ANWb or US-INe - 2018 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2019
figure()
plot(ANWb_AUG_2019,ANWb_co2_flux_AUG_2019_soy,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_AUG_2019,ANWb_soy_2019_vprm_nee_AUG,'.')
legend('obs','model')
title('ANWb or US-INe - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2020
figure()
plot(ANWb_AUG_2020,ANWb_co2_flux_AUG_2020_corn,'.','MarkerEdgeColor','black')
hold on 
plot(ANWb_AUG_2020,ANWb_corn_2020_vprm_nee_AUG,'.')
legend('obs','model')
title('ANWb or US-INe - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A09a or US-INi %%%%%%%%%%%%%
% 2019 soy, 2021 corn 
% 2019
figure()
plot(A09a_AUG_2019,A09a_co2_flux_AUG_2019_soy,'.','MarkerEdgeColor','black')
hold on 
plot(A09a_AUG_2019,A09a_soy_2019_vprm_nee_AUG,'.')
legend('obs','model')
title('A09a or US-INi - 2019 - Soy')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2021 
figure()
plot(A09a_AUG_2021,A09a_co2_flux_AUG_2021_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A09a_AUG_2021,A09a_corn_2021_vprm_nee_AUG,'.')
title('A09a or US-INi - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A09b or US-INj %%%%%%%%%%%%%
% 2020 
figure()
plot(A09b_AUG_2020,A09b_co2_flux_AUG_2020_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A09b_AUG_2020,A09b_corn_2020_vprm_nee_AUG,'.')
legend('obs','model')
title('A09b or US-INj - 2020 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%%%%%%%%%%%%% A14a or US-INn %%%%%%%%%%%%%
% 2021 corn, 2019 corn 
% 2021
figure()
plot(A14a_AUG_2021,A14a_co2_flux_AUG_2021_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A14a_AUG_2021,A14a_corn_2021_vprm_nee_AUG,'.')
legend('obs','model')
title('A14a or US-INn - 2021 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

% 2019
figure()
plot(A14a_AUG_2019,A14a_co2_flux_AUG_2019_corn,'.','MarkerEdgeColor','black')
hold on 
plot(A14a_AUG_2019,A14a_corn_2019_vprm_nee_AUG,'.')
legend('obs','model')
title('A14a or US-INn - 2019 - Corn')
ylabel('co2 flux (umol/m2s)')
grid on 

%% calculate diurnal cycles for july %% 

% test a reshape of the days 
%test_reshape_dates = reshape(ANWb_JUL_2019, 24, []); 

% use group summary - check with the code below to see if i get the
% same daily cycle just to be safe - okay i get the same cycles 
%ANWa_vprm_2018_jul_TT = timetable(ANWb_JUL_2018,ANWa_co2_flux_JUL_2018_soy,ANWa_soy_2018_vprm_nee_JUL); 
%test = groupsummary(ANWa_vprm_2018_jul_TT,"ANWb_JUL_2018","hourofday",{"mean","std","nummissing"})
%test_stderr = test.std_ANWa_co2_flux_JUL_2018_soy ./ sqrt(test.GroupCount - test.nummissing_ANWa_co2_flux_JUL_2018_soy)

%test.mean_ANWa_co2_flux_JUL_2018_soy - ANWa_flux_2018_soy_cycle

% ANWa - 2018 soy 
%ANWa_flux_2018_soy_cycle = mean(reshape(ANWa_co2_flux_JUL_2018_soy, 24, []),2,'omitnan');
%ANWa_vprm_2018_soy_cycle = mean(reshape(ANWa_soy_2018_vprm_nee_JUL, 24, []),2,'omitnan');
%ANWa_flux_2018_soy_cycle_stderr = std(reshape(ANWa_co2_flux_JUL_2018_soy, 24, []),0,2,'omitnan')./...
%    sqrt(sum(~isnan(reshape(ANWa_co2_flux_JUL_2018_soy, 24, [])),2));
%ANWa_vprm_2018_soy_cycle_stderr = std(reshape(ANWa_soy_2018_vprm_nee_JUL, 24, []),0,2,'omitnan')./...
%    sqrt(sum(~isnan(reshape(ANWa_soy_2018_vprm_nee_JUL, 24, [])),2));


%%%%%%%%%%%% ANWa - 2018 soy %%%%%%%%%%%%
% make timetable of obs and model 
ANWa_2018_jul_TT = timetable(ANWa_JUL_2018_soy,ANWa_co2_flux_JUL_2018_soy,ANWa_soy_2018_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
ANWa_jul_2018_soy_cycle = groupsummary(ANWa_2018_jul_TT,"ANWa_JUL_2018_soy","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWa_flux_jul_2018_soy_cycle_stderr = ANWa_jul_2018_soy_cycle.std_ANWa_co2_flux_JUL_2018_soy ./ ...
    sqrt(ANWa_jul_2018_soy_cycle.GroupCount - ANWa_jul_2018_soy_cycle.nummissing_ANWa_co2_flux_JUL_2018_soy); 
% calculate the standard error of the model 
ANWa_vprm_jul_2018_soy_cycle_stderr = ANWa_jul_2018_soy_cycle.std_ANWa_soy_2018_vprm_nee_JUL ./ ...
    sqrt(ANWa_jul_2018_soy_cycle.GroupCount - ANWa_jul_2018_soy_cycle.nummissing_ANWa_soy_2018_vprm_nee_JUL); 

%%%%%%%%%%%% ANWb %%%%%%%%%%%%
% ANWb - 2018 corn, 2019 soy, 2020 corn 

% make timetable of obs and model - 2018 corn 
ANWb_2018_jul_TT = timetable(ANWb_JUL_2018,ANWb_co2_flux_JUL_2018_corn,ANWb_corn_2018_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
ANWb_jul_2018_corn_cycle = groupsummary(ANWb_2018_jul_TT,"ANWb_JUL_2018","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_jul_2018_corn_cycle_stderr = ANWb_jul_2018_corn_cycle.std_ANWb_co2_flux_JUL_2018_corn ./ ...
    sqrt(ANWb_jul_2018_corn_cycle.GroupCount - ANWb_jul_2018_corn_cycle.nummissing_ANWb_co2_flux_JUL_2018_corn); 
% calculate the standard error of the model 
ANWb_vprm_jul_2018_corn_cycle_stderr = ANWb_jul_2018_corn_cycle.std_ANWb_corn_2018_vprm_nee_JUL ./ ...
    sqrt(ANWb_jul_2018_corn_cycle.GroupCount - ANWb_jul_2018_corn_cycle.nummissing_ANWb_corn_2018_vprm_nee_JUL); 

% make timetable of obs and model - 2019 soy
ANWb_2019_jul_TT = timetable(ANWb_JUL_2019,ANWb_co2_flux_JUL_2019_soy,ANWb_soy_2019_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
ANWb_jul_2019_soy_cycle = groupsummary(ANWb_2019_jul_TT,"ANWb_JUL_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_jul_2019_soy_cycle_stderr = ANWb_jul_2019_soy_cycle.std_ANWb_co2_flux_JUL_2019_soy ./ ...
    sqrt(ANWb_jul_2019_soy_cycle.GroupCount - ANWb_jul_2019_soy_cycle.nummissing_ANWb_co2_flux_JUL_2019_soy); 
% calculate the standard error of the model 
ANWb_vprm_jul_2019_soy_cycle_stderr = ANWb_jul_2019_soy_cycle.std_ANWb_soy_2019_vprm_nee_JUL ./ ...
    sqrt(ANWb_jul_2019_soy_cycle.GroupCount - ANWb_jul_2019_soy_cycle.nummissing_ANWb_soy_2019_vprm_nee_JUL); 

% make timetable of obs and model - 2020 corn
ANWb_2020_jul_TT = timetable(ANWb_JUL_2020,ANWb_co2_flux_JUL_2020_corn,ANWb_corn_2020_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
ANWb_jul_2020_corn_cycle = groupsummary(ANWb_2020_jul_TT,"ANWb_JUL_2020","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_jul_2020_corn_cycle_stderr = ANWb_jul_2020_corn_cycle.std_ANWb_co2_flux_JUL_2020_corn ./ ...
    sqrt(ANWb_jul_2020_corn_cycle.GroupCount - ANWb_jul_2020_corn_cycle.nummissing_ANWb_co2_flux_JUL_2020_corn); 
% calculate the standard error of the model 
ANWb_vprm_jul_2020_corn_cycle_stderr = ANWb_jul_2020_corn_cycle.std_ANWb_corn_2020_vprm_nee_JUL ./ ...
    sqrt(ANWb_jul_2020_corn_cycle.GroupCount - ANWb_jul_2020_corn_cycle.nummissing_ANWb_corn_2020_vprm_nee_JUL); 

%%%%%%%%%%%% A09a %%%%%%%%%%%%
% 2019 soy, 2021 corn 

% make timetable of obs and model - 2019 soy
A09a_2019_jul_TT = timetable(A09a_JUL_2019,A09a_co2_flux_JUL_2019_soy,A09a_soy_2019_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
A09a_jul_2019_soy_cycle = groupsummary(A09a_2019_jul_TT,"A09a_JUL_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09a_flux_jul_2019_soy_cycle_stderr = A09a_jul_2019_soy_cycle.std_A09a_co2_flux_JUL_2019_soy ./ ...
    sqrt(A09a_jul_2019_soy_cycle.GroupCount - A09a_jul_2019_soy_cycle.nummissing_A09a_co2_flux_JUL_2019_soy); 
% calculate the standard error of the model 
A09a_vprm_jul_2019_soy_cycle_stderr = A09a_jul_2019_soy_cycle.std_A09a_soy_2019_vprm_nee_JUL ./ ...
    sqrt(A09a_jul_2019_soy_cycle.GroupCount - A09a_jul_2019_soy_cycle.nummissing_A09a_soy_2019_vprm_nee_JUL); 

% make timetable of obs and model - 2021 corn
A09a_2021_jul_TT = timetable(A09a_JUL_2021,A09a_co2_flux_JUL_2021_corn,A09a_corn_2021_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
A09a_jul_2021_corn_cycle = groupsummary(A09a_2021_jul_TT,"A09a_JUL_2021","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09a_flux_jul_2021_corn_cycle_stderr = A09a_jul_2021_corn_cycle.std_A09a_co2_flux_JUL_2021_corn ./ ...
    sqrt(A09a_jul_2021_corn_cycle.GroupCount - A09a_jul_2021_corn_cycle.nummissing_A09a_co2_flux_JUL_2021_corn); 
% calculate the standard error of the model 
A09a_vprm_jul_2021_corn_cycle_stderr = A09a_jul_2021_corn_cycle.std_A09a_corn_2021_vprm_nee_JUL ./ ...
    sqrt(A09a_jul_2021_corn_cycle.GroupCount - A09a_jul_2021_corn_cycle.nummissing_A09a_corn_2021_vprm_nee_JUL); 

%%%%%%%%%%%% A09b - 2020 corn %%%%%%%%%%%%
% make timetable of obs and model - 2020 corn
A09b_2020_jul_TT = timetable(A09b_JUL_2020,A09b_co2_flux_JUL_2020_corn,A09b_corn_2020_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
A09b_jul_2020_corn_cycle = groupsummary(A09b_2020_jul_TT,"A09b_JUL_2020","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09b_flux_jul_2020_corn_cycle_stderr = A09b_jul_2020_corn_cycle.std_A09b_co2_flux_JUL_2020_corn ./ ...
    sqrt(A09b_jul_2020_corn_cycle.GroupCount - A09b_jul_2020_corn_cycle.nummissing_A09b_co2_flux_JUL_2020_corn); 
% calculate the standard error of the model 
A09b_vprm_jul_2020_corn_cycle_stderr = A09b_jul_2020_corn_cycle.std_A09b_corn_2020_vprm_nee_JUL ./ ...
    sqrt(A09b_jul_2020_corn_cycle.GroupCount - A09b_jul_2020_corn_cycle.nummissing_A09b_corn_2020_vprm_nee_JUL); 

%%%%%%%%%%%% A14a %%%%%%%%%%%%
% 2019 corn, 2021 corn 

% make timetable of obs and model - 2019 corn
A14a_2019_jul_TT = timetable(A14a_JUL_2019,A14a_co2_flux_JUL_2019_corn,A14a_corn_2019_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
A14a_jul_2019_corn_cycle = groupsummary(A14a_2019_jul_TT,"A14a_JUL_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A14a_flux_jul_2019_corn_cycle_stderr = A14a_jul_2019_corn_cycle.std_A14a_co2_flux_JUL_2019_corn ./ ...
    sqrt(A14a_jul_2019_corn_cycle.GroupCount - A14a_jul_2019_corn_cycle.nummissing_A14a_co2_flux_JUL_2019_corn); 
% calculate the standard error of the model 
A14a_vprm_jul_2019_corn_cycle_stderr = A14a_jul_2019_corn_cycle.std_A14a_corn_2019_vprm_nee_JUL ./ ...
    sqrt(A14a_jul_2019_corn_cycle.GroupCount - A14a_jul_2019_corn_cycle.nummissing_A14a_corn_2019_vprm_nee_JUL); 

% make timetable of obs and model - 2021 corn
A14a_2021_jul_TT = timetable(A14a_JUL_2021,A14a_co2_flux_JUL_2021_corn,A14a_corn_2021_vprm_nee_JUL); 
% calculate the mean daily cycle (obs and model)
A14a_jul_2021_corn_cycle = groupsummary(A14a_2021_jul_TT,"A14a_JUL_2021","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A14a_flux_jul_2021_corn_cycle_stderr = A14a_jul_2021_corn_cycle.std_A14a_co2_flux_JUL_2021_corn ./ ...
    sqrt(A14a_jul_2021_corn_cycle.GroupCount - A14a_jul_2021_corn_cycle.nummissing_A14a_co2_flux_JUL_2021_corn); 
% calculate the standard error of the model 
A14a_vprm_jul_2021_corn_cycle_stderr = A14a_jul_2021_corn_cycle.std_A14a_corn_2021_vprm_nee_JUL ./ ...
    sqrt(A14a_jul_2021_corn_cycle.GroupCount - A14a_jul_2021_corn_cycle.nummissing_A14a_corn_2021_vprm_nee_JUL); 



%% plot July diurnal cycles of co2 flux %% 

figure()
t = tiledlayout(3,3,"TileSpacing","compact"); 

nexttile
plot(0:1:23,ANWa_flux_2018_soy_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,ANWa_vprm_2018_soy_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INd 2018')
legend({'observations','model'},'AutoUpdate','off')
text(19,15,'(a)')
text(2,-55,'N=282')
ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,ANWb_flux_2019_soy_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,ANWb_vprm_2019_soy_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2019')
text(19,15,'(b)')
text(2,-55,'N=931')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,A09a_flux_2019_soy_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,A09a_vprm_2019_soy_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INi 2019')
text(19,15,'(c)')
text(2,-55,'N=396')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,ANWb_flux_2018_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,ANWb_vprm_2018_corn_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2018')
text(19,15,'(d)')
text(2,-55,'N=250')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,ANWb_flux_2020_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,ANWb_vprm_2020_corn_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2020')
text(19,15,'(e)')
text(2,-55,'N=1407')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,A09a_flux_2021_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,A09a_vprm_2021_corn_cycle,'.--','Color','r');
set(gca,'Xticklabel',[])
text(2,-40,'US-INi 2021')
text(19,15,'(f)')
text(2,-55,'N=87')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,A09b_flux_2020_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,A09b_vprm_2020_corn_cycle,'.--','Color','r');
%set(gca,'Xticklabel',[])
text(2,-40,'US-INj 2020')
text(19,15,'(g)')
text(2,-55,'N=355')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,A14a_flux_2019_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,A14a_vprm_2019_corn_cycle,'.--','Color','r');
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn 2019')
text(19,15,'(h)')
text(2,-55,'N=468')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

nexttile
plot(0:1:23,A14a_flux_2021_corn_cycle,'.-','Color','k'); 
hold on 
plot(0:1:23,A14a_vprm_2021_corn_cycle,'.--','Color','r');
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn 2021')
text(19,15,'(i)')
text(2,-55,'N=758')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-60 25])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})')
xlabel(t,'Hour of Day (UTC)')


%% plot july diurnal cycles of co2 flux with standard error %% 

sum(~isnan(ANWa_co2_flux_JUL_2018_soy))
sum(~isnan(ANWb_co2_flux_JUL_2018_corn))
sum(~isnan(ANWb_co2_flux_JUL_2019_soy))
sum(~isnan(ANWb_co2_flux_JUL_2020_corn))
sum(~isnan(A09a_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_JUL_2021_corn))
sum(~isnan(A09b_co2_flux_JUL_2020_corn))
sum(~isnan(A14a_co2_flux_JUL_2021_corn))
sum(~isnan(A14a_co2_flux_JUL_2019_corn))

figure()
t = tiledlayout(3,3,"TileSpacing","compact"); 

nexttile
errorbar(0:1:23,ANWa_flux_2018_soy_cycle,ANWa_flux_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_vprm_2018_soy_cycle,ANWa_vprm_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INd 2018')
legend({'observations','model'},'AutoUpdate','off')
text(19,25,'(a)')
text(2,-55,'N=53')
ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2019_soy_cycle,ANWb_flux_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2019_soy_cycle,ANWb_vprm_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2019')
text(19,25,'(b)')
text(2,-55,'N=370')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_flux_2019_soy_cycle,A09a_flux_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_vprm_2019_soy_cycle,A09a_vprm_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INi 2019')
text(19,25,'(c)')
text(2,-55,'N=77')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2018_corn_cycle,ANWb_flux_2018_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2018_corn_cycle,ANWb_vprm_2018_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2018')
text(19,25,'(d)')
text(2,-55,'N=98')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2020_corn_cycle,ANWb_flux_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2020_corn_cycle,ANWb_vprm_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2020')
text(19,25,'(e)')
text(2,-55,'N=467')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_flux_2021_corn_cycle,A09a_flux_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_vprm_2021_corn_cycle,A09a_vprm_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INi 2021')
text(19,25,'(f)')
text(2,-55,'N=12')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_flux_2020_corn_cycle,A09b_flux_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_vprm_2020_corn_cycle,A09b_vprm_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INj 2020')
text(19,25,'(g)')
text(2,-55,'N=172')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_flux_2019_corn_cycle,A14a_flux_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2019_corn_cycle,A14a_vprm_2019_corn_cycle_stderr,'.-','Color','r');
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn 2019')
text(19,25,'(h)')
text(2,-55,'N=190')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_flux_2021_corn_cycle,A14a_flux_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2021_corn_cycle,A14a_vprm_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn 2021')
text(19,25,'(i)')
text(2,-55,'N=188')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})')
xlabel(t,'Hour of Day (UTC)')




figure()
t = tiledlayout(3,3,"TileSpacing","compact"); 

nexttile
errorbar(0:1:23,ANWa_flux_2018_soy_cycle,ANWa_flux_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_vprm_2018_soy_cycle,ANWa_vprm_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INd 2018')
text(19,45,'(a)')
text(2,-55,'N=53')
ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2019_soy_cycle,ANWb_flux_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2019_soy_cycle,ANWb_vprm_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe 2019')
text(19,45,'(b)')
text(2,-55,'N=370')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_flux_2019_soy_cycle,A09a_flux_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_vprm_2019_soy_cycle,A09a_vprm_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi 2019')
text(19,45,'(c)')
text(2,-55,'N=77')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2018_corn_cycle,ANWb_flux_2018_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2018_corn_cycle,ANWb_vprm_2018_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe 2018')
text(19,45,'(d)')
text(2,-55,'N=98')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_flux_2020_corn_cycle,ANWb_flux_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_vprm_2020_corn_cycle,ANWb_vprm_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe 2020')
text(19,45,'(e)')
text(2,-55,'N=467')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_flux_2020_corn_cycle,A09b_flux_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_vprm_2020_corn_cycle,A09b_vprm_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INj 2020')
text(19,45,'(f)')
text(2,-55,'N=172')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_flux_2019_corn_cycle,A14a_flux_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2019_corn_cycle,A14a_vprm_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn 2019')
text(19,45,'(g)')
text(2,-55,'N=190')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_flux_2021_corn_cycle,A14a_flux_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2021_corn_cycle,A14a_vprm_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INn 2021')
text(19,45,'(h)')
text(2,-55,'N=188')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})')
xlabel(t,'Hour of Day (UTC)')

%% calculate AUG diurnal cycles 

%%%%%%%%%%%% ANWa %%%%%%%%%%%%
% 2018 soy 

% make timetable of obs and model 
ANWa_2018_aug_TT = timetable(ANWa_AUG_2018_soy,ANWa_co2_flux_AUG_2018_soy,ANWa_soy_2018_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
ANWa_aug_2018_soy_cycle = groupsummary(ANWa_2018_aug_TT,"ANWa_AUG_2018_soy","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWa_flux_aug_2018_soy_cycle_stderr = ANWa_aug_2018_soy_cycle.std_ANWa_co2_flux_AUG_2018_soy ./ ...
    sqrt(ANWa_aug_2018_soy_cycle.GroupCount - ANWa_aug_2018_soy_cycle.nummissing_ANWa_co2_flux_AUG_2018_soy); 
% calculate the standard error of the model 
ANWa_vprm_aug_2018_soy_cycle_stderr = ANWa_aug_2018_soy_cycle.std_ANWa_soy_2018_vprm_nee_AUG ./ ...
    sqrt(ANWa_aug_2018_soy_cycle.GroupCount - ANWa_aug_2018_soy_cycle.nummissing_ANWa_soy_2018_vprm_nee_AUG); 

%%%%%%%%%%%% ANWb %%%%%%%%%%%%
% 2018 corn, 2019 soy, 2020 corn 

% make timetable of obs and model - 2018 corn 
ANWb_2018_aug_TT = timetable(ANWb_AUG_2018,ANWb_co2_flux_AUG_2018_corn,ANWb_corn_2018_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
ANWb_aug_2018_corn_cycle = groupsummary(ANWb_2018_aug_TT,"ANWb_AUG_2018","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_aug_2018_corn_cycle_stderr = ANWb_aug_2018_corn_cycle.std_ANWb_co2_flux_AUG_2018_corn ./ ...
    sqrt(ANWb_aug_2018_corn_cycle.GroupCount - ANWb_aug_2018_corn_cycle.nummissing_ANWb_co2_flux_AUG_2018_corn); 
% calculate the standard error of the model 
ANWb_vprm_aug_2018_corn_cycle_stderr = ANWb_aug_2018_corn_cycle.std_ANWb_corn_2018_vprm_nee_AUG ./ ...
    sqrt(ANWb_aug_2018_corn_cycle.GroupCount - ANWb_aug_2018_corn_cycle.nummissing_ANWb_corn_2018_vprm_nee_AUG); 

% make timetable of obs and model - 2019 soy
ANWb_2019_aug_TT = timetable(ANWb_AUG_2019,ANWb_co2_flux_AUG_2019_soy,ANWb_soy_2019_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
ANWb_aug_2019_soy_cycle = groupsummary(ANWb_2019_aug_TT,"ANWb_AUG_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_aug_2019_soy_cycle_stderr = ANWb_aug_2019_soy_cycle.std_ANWb_co2_flux_AUG_2019_soy ./ ...
    sqrt(ANWb_aug_2019_soy_cycle.GroupCount - ANWb_aug_2019_soy_cycle.nummissing_ANWb_co2_flux_AUG_2019_soy); 
% calculate the standard error of the model 
ANWb_vprm_aug_2019_soy_cycle_stderr = ANWb_aug_2019_soy_cycle.std_ANWb_soy_2019_vprm_nee_AUG ./ ...
    sqrt(ANWb_aug_2019_soy_cycle.GroupCount - ANWb_aug_2019_soy_cycle.nummissing_ANWb_soy_2019_vprm_nee_AUG); 

% make timetable of obs and model - 2020 corn
ANWb_2020_aug_TT = timetable(ANWb_AUG_2020,ANWb_co2_flux_AUG_2020_corn,ANWb_corn_2020_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
ANWb_aug_2020_corn_cycle = groupsummary(ANWb_2020_aug_TT,"ANWb_AUG_2020","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
ANWb_flux_aug_2020_corn_cycle_stderr = ANWb_aug_2020_corn_cycle.std_ANWb_co2_flux_AUG_2020_corn ./ ...
    sqrt(ANWb_aug_2020_corn_cycle.GroupCount - ANWb_aug_2020_corn_cycle.nummissing_ANWb_co2_flux_AUG_2020_corn); 
% calculate the standard error of the model 
ANWb_vprm_aug_2020_corn_cycle_stderr = ANWb_aug_2020_corn_cycle.std_ANWb_corn_2020_vprm_nee_AUG ./ ...
    sqrt(ANWb_aug_2020_corn_cycle.GroupCount - ANWb_aug_2020_corn_cycle.nummissing_ANWb_corn_2020_vprm_nee_AUG); 

%%%%%%%%%%%% A09a %%%%%%%%%%%%
% 2019 soy, 2021 corn 

% make timetable of obs and model - 2019 soy
A09a_2019_aug_TT = timetable(A09a_AUG_2019,A09a_co2_flux_AUG_2019_soy,A09a_soy_2019_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
A09a_aug_2019_soy_cycle = groupsummary(A09a_2019_aug_TT,"A09a_AUG_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09a_flux_aug_2019_soy_cycle_stderr = A09a_aug_2019_soy_cycle.std_A09a_co2_flux_AUG_2019_soy ./ ...
    sqrt(A09a_aug_2019_soy_cycle.GroupCount - A09a_aug_2019_soy_cycle.nummissing_A09a_co2_flux_AUG_2019_soy); 
% calculate the standard error of the model 
A09a_vprm_aug_2019_soy_cycle_stderr = A09a_aug_2019_soy_cycle.std_A09a_soy_2019_vprm_nee_AUG ./ ...
    sqrt(A09a_aug_2019_soy_cycle.GroupCount - A09a_aug_2019_soy_cycle.nummissing_A09a_soy_2019_vprm_nee_AUG); 

% make timetable of obs and model - 2021 corn
A09a_2021_aug_TT = timetable(A09a_AUG_2021,A09a_co2_flux_AUG_2021_corn,A09a_corn_2021_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
A09a_aug_2021_corn_cycle = groupsummary(A09a_2021_aug_TT,"A09a_AUG_2021","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09a_flux_aug_2021_corn_cycle_stderr = A09a_aug_2021_corn_cycle.std_A09a_co2_flux_AUG_2021_corn ./ ...
    sqrt(A09a_aug_2021_corn_cycle.GroupCount - A09a_aug_2021_corn_cycle.nummissing_A09a_co2_flux_AUG_2021_corn); 
% calculate the standard error of the model 
A09a_vprm_aug_2021_corn_cycle_stderr = A09a_aug_2021_corn_cycle.std_A09a_corn_2021_vprm_nee_AUG ./ ...
    sqrt(A09a_aug_2021_corn_cycle.GroupCount - A09a_aug_2021_corn_cycle.nummissing_A09a_corn_2021_vprm_nee_AUG); 

%%%%%%%%%%%% A09b %%%%%%%%%%%%
% 2020 corn 

% make timetable of obs and model - 2020 corn
A09b_2020_aug_TT = timetable(A09b_AUG_2020,A09b_co2_flux_AUG_2020_corn,A09b_corn_2020_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
A09b_aug_2020_corn_cycle = groupsummary(A09b_2020_aug_TT,"A09b_AUG_2020","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A09b_flux_aug_2020_corn_cycle_stderr = A09b_aug_2020_corn_cycle.std_A09b_co2_flux_AUG_2020_corn ./ ...
    sqrt(A09b_aug_2020_corn_cycle.GroupCount - A09b_aug_2020_corn_cycle.nummissing_A09b_co2_flux_AUG_2020_corn); 
% calculate the standard error of the model 
A09b_vprm_aug_2020_corn_cycle_stderr = A09b_aug_2020_corn_cycle.std_A09b_corn_2020_vprm_nee_AUG ./ ...
    sqrt(A09b_aug_2020_corn_cycle.GroupCount - A09b_aug_2020_corn_cycle.nummissing_A09b_corn_2020_vprm_nee_AUG); 

%%%%%%%%%%%% A14a %%%%%%%%%%%% 
% 2019 corn, 2021 corn 

% make timetable of obs and model - 2019 corn
A14a_2019_aug_TT = timetable(A14a_AUG_2019,A14a_co2_flux_AUG_2019_corn,A14a_corn_2019_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
A14a_aug_2019_corn_cycle = groupsummary(A14a_2019_aug_TT,"A14a_AUG_2019","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A14a_flux_aug_2019_corn_cycle_stderr = A14a_aug_2019_corn_cycle.std_A14a_co2_flux_AUG_2019_corn ./ ...
    sqrt(A14a_aug_2019_corn_cycle.GroupCount - A14a_aug_2019_corn_cycle.nummissing_A14a_co2_flux_AUG_2019_corn); 
% calculate the standard error of the model 
A14a_vprm_aug_2019_corn_cycle_stderr = A14a_aug_2019_corn_cycle.std_A14a_corn_2019_vprm_nee_AUG ./ ...
    sqrt(A14a_aug_2019_corn_cycle.GroupCount - A14a_aug_2019_corn_cycle.nummissing_A14a_corn_2019_vprm_nee_AUG); 

% make timetable of obs and model - 2021 corn
A14a_2021_aug_TT = timetable(A14a_AUG_2021,A14a_co2_flux_AUG_2021_corn,A14a_corn_2021_vprm_nee_AUG); 
% calculate the mean daily cycle (obs and model)
A14a_aug_2021_corn_cycle = groupsummary(A14a_2021_aug_TT,"A14a_AUG_2021","hourofday",{"mean","std","nummissing"}); 
% calculate the standard error of obs 
A14a_flux_aug_2021_corn_cycle_stderr = A14a_aug_2021_corn_cycle.std_A14a_co2_flux_AUG_2021_corn ./ ...
    sqrt(A14a_aug_2021_corn_cycle.GroupCount - A14a_aug_2021_corn_cycle.nummissing_A14a_co2_flux_AUG_2021_corn); 
% calculate the standard error of the model 
A14a_vprm_aug_2021_corn_cycle_stderr = A14a_aug_2021_corn_cycle.std_A14a_corn_2021_vprm_nee_AUG ./ ...
    sqrt(A14a_aug_2021_corn_cycle.GroupCount - A14a_aug_2021_corn_cycle.nummissing_A14a_corn_2021_vprm_nee_AUG); 


%% plot july and aug together for a mega figure 

sum(~isnan(ANWa_co2_flux_AUG_2018_soy))
sum(~isnan(ANWb_co2_flux_AUG_2018_corn))
sum(~isnan(ANWb_co2_flux_AUG_2019_soy))
sum(~isnan(ANWb_co2_flux_AUG_2020_corn))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))
sum(~isnan(A09a_co2_flux_AUG_2021_corn))
sum(~isnan(A09b_co2_flux_AUG_2020_corn))
sum(~isnan(A14a_co2_flux_AUG_2021_corn))
sum(~isnan(A14a_co2_flux_AUG_2019_corn))

sum(~isnan(ANWa_co2_flux_JUL_2018_soy))
sum(~isnan(ANWa_co2_flux_AUG_2018_soy))
sum(~isnan(ANWb_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))
sum(~isnan(ANWb_co2_flux_JUL_2018_corn))
sum(~isnan(ANWb_co2_flux_JUL_2020_corn))

figure()
t = tiledlayout(5,3,"TileSpacing","compact"); 

nexttile
errorbar(0:1:23,ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy,ANWa_flux_jul_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL,ANWa_vprm_jul_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INd Jul 2018')
text(21,45,'(a)')
text(2,-55,'N=46')
ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy,ANWa_flux_aug_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG,ANWa_vprm_aug_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INd Aug 2018')
text(21,45,'(b)')
text(2,-55,'N=145')
%ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy,ANWb_flux_jul_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL,ANWb_vprm_jul_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2019')
text(21,45,'(c)')
text(2,-55,'N=373')
legend({'observations','model'},'AutoUpdate','off')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy,ANWb_flux_aug_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG,ANWb_vprm_aug_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Aug 2019')
text(21,45,'(d)')
text(2,-55,'N=276')
ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy,A09a_flux_jul_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL,A09a_vprm_jul_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Jul 2019')
text(21,45,'(e)')
text(2,-55,'N=70')
%ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy,A09a_flux_aug_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG,A09a_vprm_aug_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Aug 2019')
text(21,45,'(f)')
%ylabel('Soybean')
text(2,-55,'N=151')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn,ANWb_flux_jul_2018_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL,ANWb_vprm_jul_2018_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Jul 2018')
text(21,45,'(g)')
text(2,-55,'N=106')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
errorbar(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn,ANWb_flux_jul_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL,ANWb_vprm_jul_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2020')
text(21,45,'(h)')
text(2,-55,'N=432')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn,ANWb_flux_aug_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG,ANWb_vprm_aug_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Aug 2020')
text(21,45,'(i)')
text(2,-55,'N=422')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn,A09b_flux_jul_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL,A09b_vprm_jul_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INj Jul 2020')
text(21,45,'(j)')
text(2,-55,'N=159')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn,A09b_flux_aug_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG,A09b_vprm_aug_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INj Aug 2020')
text(21,45,'(k)')
text(2,-55,'N=47')
%ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
errorbar(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn,A14a_flux_jul_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL,A14a_vprm_jul_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn Jul 2019')
text(21,45,'(l)')
set(gca,'Yticklabel',[])
text(2,-55,'N=190')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

%{
nexttile
errorbar(0:1:23,A14a_flux_2019_corn_cycle,A14a_flux_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2019_corn_cycle,A14a_vprm_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn Aug 2019')
text(19,45,'(g)')
text(2,-55,'N=44')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)
%} 
nexttile
errorbar(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn,A14a_flux_jul_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL,A14a_vprm_jul_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INn Jul 2021')
text(21,45,'(m)')
text(2,-55,'N=183')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
ylabel('Corn')
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn,A14a_flux_aug_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG,A14a_vprm_aug_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INn Aug 2021')
text(21,45,'(n)')
text(2,-55,'N=197')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})')
xlabel(t,'Hour of Day (UTC)')

%% add y-axis to show number of points in hour %%

figure()
t = tiledlayout(5,3,"TileSpacing","compact"); 

nexttile
errorbar(0:1:23,ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy,ANWa_flux_jul_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL,ANWa_vprm_jul_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INd Jul 2018')
text(21,45,'(a)')
text(2,-55,'N=46')
ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy,ANWa_flux_aug_2018_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG,ANWa_vprm_aug_2018_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INd Aug 2018')
text(21,45,'(b)')
text(2,-55,'N=145')
%ylabel('Soybean')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy,ANWb_flux_jul_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL,ANWb_vprm_jul_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2019')
text(21,45,'(c)')
text(2,-55,'N=373')
legend({'observations','model'},'AutoUpdate','off')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy,ANWb_flux_aug_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG,ANWb_vprm_aug_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Aug 2019')
text(21,45,'(d)')
text(2,-55,'N=276')
ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy,A09a_flux_jul_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL,A09a_vprm_jul_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Jul 2019')
text(21,45,'(e)')
text(2,-55,'N=70')
%ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy,A09a_flux_aug_2019_soy_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG,A09a_vprm_aug_2019_soy_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Aug 2019')
text(21,45,'(f)')
%ylabel('Soybean')
text(2,-55,'N=151')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn,ANWb_flux_jul_2018_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL,ANWb_vprm_jul_2018_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Jul 2018')
text(21,45,'(g)')
text(2,-55,'N=106')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
errorbar(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn,ANWb_flux_jul_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL,ANWb_vprm_jul_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2020')
text(21,45,'(h)')
text(2,-55,'N=432')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn,ANWb_flux_aug_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG,ANWb_vprm_aug_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Aug 2020')
text(21,45,'(i)')
text(2,-55,'N=422')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn,A09b_flux_jul_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL,A09b_vprm_jul_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INj Jul 2020')
text(21,45,'(j)')
text(2,-55,'N=159')
ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn,A09b_flux_aug_2020_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG,A09b_vprm_aug_2020_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INj Aug 2020')
text(21,45,'(k)')
text(2,-55,'N=47')
%ylabel('Corn')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
errorbar(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn,A14a_flux_jul_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL,A14a_vprm_jul_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn Jul 2019')
text(21,45,'(l)')
set(gca,'Yticklabel',[])
text(2,-55,'N=190')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

%{
nexttile
errorbar(0:1:23,A14a_flux_2019_corn_cycle,A14a_flux_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2019_corn_cycle,A14a_vprm_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn Aug 2019')
text(19,45,'(g)')
text(2,-55,'N=44')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)
%} 
nexttile
errorbar(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn,A14a_flux_jul_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL,A14a_vprm_jul_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INn Jul 2021')
text(21,45,'(m)')
text(2,-55,'N=183')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
ylabel('Corn')
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
errorbar(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn,A14a_flux_aug_2021_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG,A14a_vprm_aug_2021_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INn Aug 2021')
text(21,45,'(n)')
text(2,-55,'N=197')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})')
xlabel(t,'Hour of Day (UTC)')


%% make july and august figure with shading %% 

sum(~isnan(ANWa_co2_flux_AUG_2018_soy))
sum(~isnan(ANWb_co2_flux_AUG_2018_corn))
sum(~isnan(ANWb_co2_flux_AUG_2019_soy))
sum(~isnan(ANWb_co2_flux_AUG_2020_corn))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))
sum(~isnan(A09a_co2_flux_AUG_2021_corn))
sum(~isnan(A09b_co2_flux_AUG_2020_corn))
sum(~isnan(A14a_co2_flux_AUG_2021_corn))
sum(~isnan(A14a_co2_flux_AUG_2019_corn))

sum(~isnan(ANWa_co2_flux_JUL_2018_soy))
sum(~isnan(ANWa_co2_flux_AUG_2018_soy))
sum(~isnan(ANWb_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_JUL_2019_soy))
sum(~isnan(A09a_co2_flux_AUG_2019_soy))
sum(~isnan(ANWb_co2_flux_JUL_2018_corn))
sum(~isnan(ANWb_co2_flux_JUL_2020_corn))


figure()
t = tiledlayout(5,3,"TileSpacing","compact"); 

nexttile
%yyaxis left
boundedline((0:1:23),ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy,ANWa_flux_jul_2018_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL,ANWa_vprm_jul_2018_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INd Jul 2018')
text(20,40,'(a)')
text(2,-55,'N=42')
ylabel('Soybean')
xlim([0 23])
ylim([-80 60])
yline(0)
grid on 
box on 
%yyaxis right 
%plot((0:1:23),ANWa_jul_2018_soy_cycle.Groupcount - ANWa_jul_2018_soy_cycle.nummissing_ANWa_soy_2018_vprm_nee_JUL)
%ylim([0 10])

nexttile
boundedline(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy,ANWa_flux_aug_2018_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG,ANWa_vprm_aug_2018_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INd Aug 2018')
text(20,40,'(b)')
text(2,-55,'N=124')
%ylabel('Soybean')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy,ANWb_flux_jul_2019_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL,ANWb_vprm_jul_2019_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2019')
text(20,40,'(c)')
text(2,-55,'N=303')
legend({'','observations','','model'},'AutoUpdate','off')
legend box off
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy,ANWb_flux_aug_2019_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG,ANWb_vprm_aug_2019_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Aug 2019')
text(20,40,'(d)')
text(2,-55,'N=218')
ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy,A09a_flux_jul_2019_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL,A09a_vprm_jul_2019_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Jul 2019')
text(20,40,'(e)')
text(2,-55,'N=63')
%ylabel('Soybean')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy,A09a_flux_aug_2019_soy_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG,A09a_vprm_aug_2019_soy_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
%legend({'observations','model'},'AutoUpdate','off')
text(2,-40,'US-INi Aug 2019')
text(20,40,'(f)')
%ylabel('Soybean')
text(2,-55,'N=141')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn,ANWb_flux_jul_2018_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL,ANWb_vprm_jul_2018_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
text(2,-40,'US-INe Jul 2018')
text(20,40,'(g)')
text(2,-55,'N=83')
ylabel('Corn')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
boundedline(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn,ANWb_flux_jul_2020_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL,ANWb_vprm_jul_2020_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Jul 2020')
text(20,40,'(h)')
text(2,-55,'N=432')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn,ANWb_flux_aug_2020_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG,ANWb_vprm_aug_2020_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INe Aug 2020')
text(20,40,'(i)')
text(2,-55,'N=422')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn,A09b_flux_jul_2020_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL,A09b_vprm_jul_2020_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INj Jul 2020')
text(20,40,'(j)')
text(2,-55,'N=152')
ylabel('Corn')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn,A09b_flux_aug_2020_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG,A09b_vprm_aug_2020_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INj Aug 2020')
text(20,40,'(k)')
text(2,-55,'N=46')
%ylabel('Corn')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)


nexttile
boundedline(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn,A14a_flux_jul_2019_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL,A14a_vprm_jul_2019_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
%set(gca,'Xticklabel',[]
text(2,-40,'US-INn Jul 2019')
text(20,40,'(l)')
set(gca,'Yticklabel',[])
text(2,-55,'N=165')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

%{
nexttile
errorbar(0:1:23,A14a_flux_2019_corn_cycle,A14a_flux_2019_corn_cycle_stderr,'.-','Color','k','LineWidth',1); 
hold on 
errorbar(0:1:23,A14a_vprm_2019_corn_cycle,A14a_vprm_2019_corn_cycle_stderr,'.-','Color','r','LineWidth',1);
%set(gca,'Xticklabel',[])
text(2,-40,'US-INn Aug 2019')
text(19,45,'(g)')
text(2,-55,'N=44')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
xlim([0 23])
ylim([-80 60])
yline(0)
%} 
nexttile
boundedline(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn,A14a_flux_jul_2021_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL,A14a_vprm_jul_2021_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
%set(gca,'Xticklabel',[])
%set(gca,'Yticklabel',[])
text(2,-40,'US-INn Jul 2021')
text(20,40,'(m)')
text(2,-55,'N=154')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
ylabel('Corn')
xlim([0 23])
ylim([-80 60])
yline(0)

nexttile
boundedline(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn,A14a_flux_aug_2021_corn_cycle_stderr,'nan','gap','alpha','-k.','Linewidth',1); 
hold on 
boundedline(0:1:23,A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG,A14a_vprm_aug_2021_corn_cycle_stderr,'nan','gap','alpha','-r.','Linewidth',1);
%set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
text(2,-40,'US-INn Aug 2021')
text(20,40,'(n)')
text(2,-55,'N=175')
%ylabel('CO_2 flux (/mu mol m^{-2} s^{-1}')
grid on 
box on 
xlim([0 23])
ylim([-80 60])
yline(0)

ylabel(t,'Mean CO_2 Flux (\mumolm^{-2}s^{-1})','FontSize',12)
xlabel(t,'Hour of Day (UTC)','FontSize',12) 


%% calculate the percent differences between model and obs peak drawdown and resp %% 

drawdown_peaks_obs = [min(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy),...
    min(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy), ...
    min(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy), ...
    min(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy), ...
    min(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy), ...
    min(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy), ...
    min(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn), ...
    min(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn), ...
    min(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn), ...
    min(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn), ...
    min(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn), ...
    min(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn), ...
    min(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn), ...
    min(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn)]; 

min(drawdown_peaks_obs)
max(drawdown_peaks_obs)

resp_peaks_obs = [max(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy),...
    max(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy), ...
    max(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy), ...
    max(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy), ...
    max(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy), ...
    max(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy), ...
    max(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn), ...
    max(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn), ...
    max(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn), ...
    max(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn), ...
    max(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn), ...
    max(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn), ...
    max(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn), ...
    max(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn)]; 

min(resp_peaks_obs)
max(resp_peaks_obs)


percent_diff_model_obs_drawdown_peak = [ abs((min(ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL)-min(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))/min(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))*100, ...
    abs((min(ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG)-min(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))/min(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))*100, ...
    abs((min(ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL)-min(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))/min(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))*100, ...
    abs((min(ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG)-min(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))/min(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))*100, ...
    abs((min(A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL)-min(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))/min(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))*100, ...
    abs((min(A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG)-min(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))/min(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))*100, ...
    abs((min(ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL)-min(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))/min(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))*100, ...
    abs((min(ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL)-min(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))/min(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))*100, ...
    abs((min(ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG)-min(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))/min(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))*100, ...
    abs((min(A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL)-min(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))/min(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))*100, ...
    abs((min(A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG)-min(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))/min(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))*100, ...
    abs((min(A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL)-min(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))/min(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))*100, ...
    abs((min(A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL)-min(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))/min(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))*100, ...
    abs((min(A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG)-min(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))/min(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))*100]; 

min(percent_diff_model_obs_drawdown_peak)
max(percent_diff_model_obs_drawdown_peak)

% remove us-inn / 14a 2019 

percent_diff_model_obs_drawdown_peak_2 = [ abs((min(ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL)-min(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))/min(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))*100, ...
    abs((min(ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG)-min(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))/min(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))*100, ...
    abs((min(ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL)-min(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))/min(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))*100, ...
    abs((min(ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG)-min(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))/min(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))*100, ...
    abs((min(A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL)-min(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))/min(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))*100, ...
    abs((min(A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG)-min(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))/min(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))*100, ...
    abs((min(ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL)-min(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))/min(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))*100, ...
    abs((min(ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL)-min(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))/min(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))*100, ...
    abs((min(ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG)-min(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))/min(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))*100, ...
    abs((min(A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL)-min(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))/min(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))*100, ...
    abs((min(A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG)-min(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))/min(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))*100, ...
    abs((min(A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL)-min(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))/min(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))*100, ...
    abs((min(A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG)-min(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))/min(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))*100]; 

min(percent_diff_model_obs_drawdown_peak_2)
max(percent_diff_model_obs_drawdown_peak_2)



percent_diff_model_obs_resp_peak = [ abs((max(ANWa_jul_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_JUL)-max(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))/max(ANWa_jul_2018_soy_cycle.mean_ANWa_co2_flux_JUL_2018_soy))*100, ...
    abs((max(ANWa_aug_2018_soy_cycle.mean_ANWa_soy_2018_vprm_nee_AUG)-max(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))/max(ANWa_aug_2018_soy_cycle.mean_ANWa_co2_flux_AUG_2018_soy))*100, ...
    abs((max(ANWb_jul_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_JUL)-max(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))/max(ANWb_jul_2019_soy_cycle.mean_ANWb_co2_flux_JUL_2019_soy))*100, ...
    abs((max(ANWb_aug_2019_soy_cycle.mean_ANWb_soy_2019_vprm_nee_AUG)-max(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))/max(ANWb_aug_2019_soy_cycle.mean_ANWb_co2_flux_AUG_2019_soy))*100, ...
    abs((max(A09a_jul_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_JUL)-max(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))/max(A09a_jul_2019_soy_cycle.mean_A09a_co2_flux_JUL_2019_soy))*100, ...
    abs((max(A09a_aug_2019_soy_cycle.mean_A09a_soy_2019_vprm_nee_AUG)-max(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))/max(A09a_aug_2019_soy_cycle.mean_A09a_co2_flux_AUG_2019_soy))*100, ...
    abs((max(ANWb_jul_2018_corn_cycle.mean_ANWb_corn_2018_vprm_nee_JUL)-max(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))/max(ANWb_jul_2018_corn_cycle.mean_ANWb_co2_flux_JUL_2018_corn))*100, ...
    abs((max(ANWb_jul_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_JUL)-max(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))/max(ANWb_jul_2020_corn_cycle.mean_ANWb_co2_flux_JUL_2020_corn))*100, ...
    abs((max(ANWb_aug_2020_corn_cycle.mean_ANWb_corn_2020_vprm_nee_AUG)-max(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))/max(ANWb_aug_2020_corn_cycle.mean_ANWb_co2_flux_AUG_2020_corn))*100, ...
    abs((max(A09b_jul_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_JUL)-max(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))/max(A09b_jul_2020_corn_cycle.mean_A09b_co2_flux_JUL_2020_corn))*100, ...
    abs((max(A09b_aug_2020_corn_cycle.mean_A09b_corn_2020_vprm_nee_AUG)-max(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))/max(A09b_aug_2020_corn_cycle.mean_A09b_co2_flux_AUG_2020_corn))*100, ...
    abs((max(A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL)-max(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))/max(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))*100, ...
    abs((max(A14a_jul_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_JUL)-max(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))/max(A14a_jul_2021_corn_cycle.mean_A14a_co2_flux_JUL_2021_corn))*100, ...
    abs((max(A14a_aug_2021_corn_cycle.mean_A14a_corn_2021_vprm_nee_AUG)-max(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))/max(A14a_aug_2021_corn_cycle.mean_A14a_co2_flux_AUG_2021_corn))*100]; 

max(percent_diff_model_obs_resp_peak)
min(percent_diff_model_obs_resp_peak)

abs((max(A14a_jul_2019_corn_cycle.mean_A14a_corn_2019_vprm_nee_JUL)-max(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))/max(A14a_jul_2019_corn_cycle.mean_A14a_co2_flux_JUL_2019_corn))*100, ...


