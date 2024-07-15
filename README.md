
Sebastian Sippel
24.10.2023

This is a code repository for dynamical adjustment, implementing the method presented in: 
Sippel, S., Meinshausen, N., Merrifield, A., Lehner, F., Pendergrass, A. G., Fischer, E. M., and Knutti, R. (2019) Uncovering the forced climate response from a single ensemble member using statistical learning. Journal of Climate 32, 5677-5699, doi:10.1175/JCLI-D-18-0882.1.

The code is used in the following paper that is accepted at _Weather and Climate Dynamics_: Sippel, S., Barnes, C., Cadiou, C., Fischer, E., Kew, S., Kretschmer, M., Philip, S., Shepherd, T. G., Singh, J., Vautard, R., and Yiou, P.: Could an extremely cold central European winter such as 1963
happen again despite climate change? _Weather and Climate Dynamics_ (accepted), 2024. Preprint: https://doi.org/10.5194/egusphere-2023-2523.

Short summary of the functions in each script:

"_00_preprocessing_4DYNADJ.R"   Function repository for preprocessing based on time series (xts) package in R for clean data handling.
"_00_preprocessing_SPACETIME_4DYNADJ.R" 
"_01_GRIDCELL_DYNADJ.R"
"_01_SPACE_DYNADJ.R"
"_02_ANALYSISFUN_DYNADJ.R"
"_02_ANALYSISFUN_DYNADJ_RB.R"



The script "_03_simple_test_script.R" implements a simple version (file paths need to be adjusted).
The data files are available from: 

FTP Server: data.iac.ethz.ch/sippels/dynamical_adjustment_elasticnet
Username: climphys
Password: friendly

# Internal use:
#cp /net/h2o/climphys1/sippels/_DATA/cesm122/monthly_control/PSL_Europe_2000y_NDJF_anom.nc /net/radon/climphys/FTP/sippels/dynamical_adjustment_elasticnet/
#cp /net/h2o/climphys1/sippels/_DATA/cesm122/monthly_control/TREFHT_Europe_2000y_NDJF_anom.nc /net/radon/climphys/FTP/sippels/dynamical_adjustment_elasticnet/
