# Stat Learning Dynamical Adjustment - 

setwd("/home/meranna/R/SL_Dynamical_Adjustment")

# packages
require(ncdf4)
require(xts)
require(zoo)
require(raster)

install.packages("ncdf.tools")
require("ncdf.tools")

install.packages("glmnet")
require("glmnet")

#install.packages("R.matlab")
#require("R.matlab")

source("_00_preprocessing_4DYNADJ.R")
source("_01_GRIDCELL_DYNADJ.R")
source("_02_ANALYSISFUN_DYNADJ.R")
source("_02_ANALYSISFUN_DYNADJ_RB.R")

## read in netcdfs

# Control Run series of Januarys
ncin <- nc_open("b.e104.B_1850_CN.f19_g16.control-im64.340.cam2.h0.0340-1321.ac.nc")
# b.e104.B_1850_CN.f19_g16.control-im64.340.cam2.h0.0340-1321.ac.nc"   # Atlantic centered map

print(ncin)

## method a to read in variables 
  # TREFHT <- ncvar_get(ncin,"TREFHT")
  # str(TREFHT)
  # image(TREFHT[,,1])

## method b to read in variables 
TREFHT=brick("b.e104.B_1850_CN.f19_g16.control-im64.340.cam2.h0.0340-1321.ac.nc", varname="TREFHT")-0
PSL=brick("b.e104.B_1850_CN.f19_g16.control-im64.340.cam2.h0.0340-1321.ac.nc", varname="PSL")-0

Tflip=brick("b.e104.B_1850_CN.f19_g16.2xCO2.420.cam2.h0.1000-2923.ac.nc", varname="Tflip")-0
Pflip=brick("b.e104.B_1850_CN.f19_g16.2xCO2.420.cam2.h0.1000-2923.ac.nc", varname="Pflip")-0

#Tflip=brick("b.e104.B_1850_CN.f19_g16.4xCO2.420.cam2.h0.5000-6329.ac.nc", varname="Tflip")-0
#Pflip=brick("b.e104.B_1850_CN.f19_g16.4xCO2.420.cam2.h0.5000-6329.ac.nc", varname="Pflip")-0

#Tflip=brick("b.e104.B_1850_CN.f19_g16.8xCO2.420.cam2.h0.4000-5589.ac.nc", varname="Tflip")-0
#Pflip=brick("b.e104.B_1850_CN.f19_g16.8xCO2.420.cam2.h0.4000-5589.ac.nc", varname="Pflip")-0

## Construct Inputs 
  # Y.train = TREFHT time series for i.e. London (51N,0)
  # X = PSL, subsetted over the DA domain

# fix ncar date lag issue for Y.train
time.info.T = as.Date(substring(names(TREFHT), first =  2), "%Y.%m.%d")-15  
Y.train.ctl = xts(c(extract(TREFHT, data.frame(x=37.5,y=61.5789))), order.by = time.info.T) # CTL

time.info.T = as.Date(substring(names(Tflip), first =  2), "%Y.%m.%d")-15  
Y.train.2x = xts(c(extract(Tflip, data.frame(x=37.5,y=61.5789))), order.by = time.info.T) # 2xCO2

X.train = crop(PSL,extent(c(-60,45,25,90))) # CTL slp
X = crop(Pflip,extent(c(-60,45,25,90))) # 2xCO2 SLP

# fix ncar date lag issue in X
names(X) = paste("X", format(as.Date(substring(names(X), 2), "%Y.%m.%d")-15, "%Y.%m.%d"), sep="") # replace time
names(X.train) = paste("X.train", format(as.Date(substring(names(X.train), 2), "%Y.%m.%d")-15, "%Y.%m.%d"), sep="") # replace time

# cut to analogue selection periods
X.ctl.cut = extract.period.RB(RB = X.train, years = 340:1321, months = 1:12)
X.2x.cut = extract.period.RB(RB = X, years = 1942:2923, months = 1:12)

# For Apples to Apples Comparison
# CTL: 340-01 to 1321-12
# 2x: 1942-01 to 2923-12
# 4x: 5348-01 to 6329-12
# 8x: 4608-01 to 5589-12

# Initialize and Call Functions in GRIDCELL file
#extract.years.from.ts(Y.train,4608:5589)
#extract.period.RB(RB = X, years = 4608:5589, months = 1:12) 

Y.hat = train.dyn.adj.elastnet.annual(Y.train = extract.years.from.ts(Y.train.ctl,340:1321), X = X.2x.cut , X.train = X.ctl.cut, X.date.seq = NULL, 
                                     train.years = NULL, train.months = 1:12, add.mon = 0,
                                     alpha = 0, nfolds = 10, ret.cv.model=F, s = "lambda.min",
                                     lags.to.aggregate = NULL, n.pc = 10)
                                     

## Output to Matlab

#writeMat("MLDA_4x_v3.mat",T_zrh = Y.train ,T_dyn_zrh = Y.hat)