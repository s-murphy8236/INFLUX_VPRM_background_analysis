Summary of scripts and files in this directory 

— Scripts —
domain_landcover.R - calculate the percentage of the domain land cover that each model Plant Functional Type (PFT) makes up 

domain_PFT_fraction_plots.m - plot the percentage of the domain that each model Plant Functional Type (PFT) makes up 

extract_vals_from_VPRM_outputs.R - extract PFT-specific VPRM NEE results at flux tower sites and save to a csv file for further analysis 

inf_func_VPRM_convolution_script_UPDATED20220809.R - convolve tower influence functions and VPRM outputs (regridded to match influence function spatial resolution) and save out the modeled biogenic CO2 mole fractions to CSV files, this script runs one tower for one month 

inf_func_raster_stackings.R - read in the tower influence functions from the netCDF files, turn them into rasters in WGS1984 coordinates (to match the VPRM NEE outputs), and stack them to be used more easily in the VPRM-influence function convolutions, this script runs one tower for one month 

regrid_VPRM_outputs_for_convolution_UPDATED20220803.R - read in VPRM outputs, regridd them to match the wrf outputs (i.e., convert the 266x354 VPRM NEE layers to 99x99 to match the influence function grids), this script runs one tower for one month 

stack_VPRM_outputs.R - purpose: read in VPRM NEE in the 99x99 grid, stack VPRM files to make the convolution easier, this script runs one month at a time 

VPRM_gridded_script_UPDATED.R - run VPRM (updated structure from Gourdji et al., 2022)

model_data_comparison_CLEANEDUP.m - analysis for model observation comparisons of CO2 enhancements, as well as analysis on model mechanisms impacting its ability to capture observations 

obs_vprm_comparison_JULAUG_cycles.m - analysis code to compare observed daily cycles of co2 flux at agricultural sites to modeled fluxes from VPRM 

flux_footprint_filtering_fraction_calculation_UPDATED.m - calculate fraction of the footprint attributable to the crop of interest (corn/soy) using footprint model from Kljun et al., 2015. 

flux_footprint_filtering_fraction_calculation_UPDATED_ANWa_2018.m - calculate fraction of the footprint attributable to the crop of interest (corn/soy) using footprint model from Kljun et al., 2015. This code is specifically for one site during 2018, when the site was split by corn on one side and soy on another. 


— Other files —
site_map_csvs - csv files that indicate the location of the crop of interest vs other land cover at each agricultural flux sites used in this study (NOTE: the names of the sites on the files are different from the Ameriflux site names used in this paper, see site name key included below)

vprm_params_20210727.csv - csv file with parameters for the VPRM runs (parameters from Gourdji et al., 2022)

PFT_yearly_fractions.csv - fraction of the domain land cover that each model Plant Functional Type (PFT) makes up 

temporary_flux_files_*.zip - zip files containing CSVs of CO2 flux data used for analysis (this is temporarily available here until the data is published)

— Flux Site Name Key —

A09a = Ameriflux US-INi
A09b = Ameriflux US-INj
A14a = Ameriflux US-INn
Anwa = Ameriflux US-INd
Anwb = Ameriflux US-INe



References noted here:
Gourdji, S. M., Karion, A., Lopez-coto, I., Ghosh, S., Mueller, K. L., Zhou, Y., et al. (2022). A Modified Vegetation Photosynthesis and Respiration Model (VPRM) for the Eastern USA and Canada, Evaluated With Comparison to Atmospheric Observations and Other Biospheric Models. Journal of Geophysical Research: Biogeosciences. https://doi.org/10.1029/2021JG006290

Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP). Geoscientific Model Development, 8(11), 3695–3713. https://doi.org/10.5194/gmd-8-3695-2015
