
#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## install necessary packages
bash analysis/packages.sh

## first create a folder for saving the occupancy model results and tables
mkdir -p ../⁨speciesRoles_saved⁩/⁨occupancy⁩/
mkdir -p ../⁨speciesRoles_saved⁩/⁨occupancy⁩/tables

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github
Rscript ../dataPrep/dataPrep.R

##***************************************************************
## interaction flexibility calculating and analysis
##***************************************************************
## pollinator partner variability
Rscript analysis/variability/1nulls.R "pols" 999
Rscript analysis/variability/2partner.R "pols" "abund"
## pollinator role variability
Rscript analysis/variability/3role.R
## pollinator structural variability
Rscript analysis/variability/4nestedContr.R
## path analysis
Rscript analysis/variability/5piecewiseSEM.R

##***************************************************************
## Plotting supplementary figures
##***************************************************************
Rscript analysis/variability/plotting/allValues.R

##***************************************************************
## Occupancy models
##***************************************************************
## Running the models
## the next script takes quite a few hours to run
Rscript analysis/occupancy/1runModels.R
## calculating posteriors
Rscript analysis/occupancy/2posteriors.R
