# Discharge_area_scaling

#data
The MATLAB file named KentuckyData.mat consists of climatic and streamflow data along with the upstream and downstream lengths of the channelized pixels.
Column 1- Rainfall (mm/day)
Column 2- Avg. temperature (deg. Celcius)
Column 3- PET (mm/day)
Column 4- Streamflow (mm/day)
Column 5- upstream length (km)
Column 6- downstream length (km)
Column 7- flow direction
The file named Ken_parameters.mat consists of calibrated model parameters for the Kentucky basin. 

%% For the events-based analysis, the synthetic rainfall of known characteristics were generated. 
First column of -Data_kensynt4ev1.mat consists of 4-hr synthetic uniform effective rainfall events of 5 mm magnitude for 23 nested basins in Kentucky,  Cloumn 4-7 are the same as above.
First column of -Data_kensynioev1.mat consists of 4-hr synthetic effective rainfall that gradually decreases towards downstream
Similarly, the first column of -Data_kensyndoev1.mat consists of 4-hr synthetic effective rainfall that gradually decreases towards upstream.
Excel file named BasinAreaAll.xlsx consists of basin areas of all the Kentucky basins.

# Model scripts
MATLAB script under the name DB_GHRM.m consist of Dynamic Budyko-GHRM model framework to simulate the discharge time series and calculates the scaling exponent for observed and modelled discharge.
DBM is a MATLAB function that is called in DB_GHRM script to estimate the effective rainfall (ER) to be routed using GHRM.
GHRM_event.m calculates the scaling exponent for event-based simulations with different synthetic effective rainfall as input.
