#### PolyRx main script ####
##  This main script sources Polypharmacy (PolyRx) scripts to load data, prep data, and create figures and tables: end of project (eop)
## Sys.Date() #"2025-03-25"
## author: christian brieghel

## READ ME: 
## Scripts should be run in listed order, but IPI_eop.R, DALYCARE_entities_eop.R, PolyRx_CCI_eop.R
## saves csv files to work directory, so they may be loaded (and above scripts by-passed) in 
## PolyRx_data_prep.R

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
# SAVE = TRUE #Saves all Figure tables if TRUE
# setwd('insert_your_wd_for_where_to_save_output_if_SAVE = T') #
setwd('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/') # E.g.

source('/ngc/projects2/dalyca_r/chribr_r/PolyRx/scripts/eop/upload/IPI_eop.R') # calculates IPI scores: 5 min. OK! Saves as:
# By-passes the above script that saves csv files as
# IPI = read_csv2('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/IPI_2025.csv') 

source('/ngc/projects2/dalyca_r/chribr_r/PolyRx/scripts/eop/upload/DALYCARE_entities_eop.R') # creates dalycare entities #ETA 0.1 min # Needs update
source('/ngc/projects2/dalyca_r/chribr_r/PolyRx/scripts/eop/upload/PolyRx_CCI_eop.R') #Calculates CCI scores and polypharmacy # ETA 6-8 min 

# By-passes the above scripts that saves csv files loaded in data_prep.R
source('/ngc/projects2/dalyca_r/chribr_r/PolyRx/scripts/eop/PolyRx_data_prep2.R') # preps all data # ETA

# Saves all figures and tables to work directory here: 
getwd()
source('/ngc/projects2/dalyca_r/chribr_r/PolyRx/scripts/eop/PolyRx_Figures_eop.R') # creates and saves all figures and tables

## End of script