
# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment simple script:
# ----------------------------------------------------------------------------------------------

# Sebastian Sippel
# 18.02.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(fields)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_00_preprocessing_4DYNADJ.R")
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_00_preprocessing_SPACETIME_4DYNADJ.R")
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_01_GRIDCELL_DYNADJ.R")
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_01_SPACE_DYNADJ.R")
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ.R")
source("/net/h2o/climphys1/sippels/_code/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ_RB.R")

# other functions to read:
source("/net/h2o/climphys1/sippels/_code/tools/convert.to.eurocentric.R")
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_code/tools/project_raster.R")

source("/net/h2o/climphys1/sippels/_code/tools/gridcorts.R")



# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_DATA/cesm122/monthly_control/")

PSL_EUROPE = brick("PSL_Europe_2000y_NDJF_anom.nc") + 0
TREFHT_EUROPE = brick("TREFHT_Europe_2000y_NDJF_anom.nc") + 0

swiss.extent = extent(c(-10, 30, 25, 65))
PSL_CH = crop(PSL_EUROPE, swiss.extent)


# ------------------------------------------------------------------------------------------ 
# 2. TRAIN MODEL FOR ONE TIME SERIES:
# ------------------------------------------------------------------------------------------ 

# select Zurich grid point:
TREFHT_ZH = xts(c(raster::extract(x = TREFHT_EUROPE, data.frame(cbind(x = 8.5, y = 47.4)))), order.by = extract.date.RB(RB = TREFHT_EUROPE))


# run dynamical adjustment training:
TREFHT_ZH_hat = train.dyn.adj.elastnet.annual(Y.train = TREFHT_ZH, X = PSL_CH, train.years = 1000:2000, train.months = 1, add.mon = 1, alpha = 0.2, nfolds = 10, s = "lambda.1se", 
                                              lags.to.aggregate = list(1), n.pc = 10, cv.parallel = T)
TREFHT_ZH_hat1 = extract.years.from.ts(ts = extract.months.from.ts(TREFHT_ZH_hat, months=1), years = 2001:3000)
TREFHT_ZH1 = extract.years.from.ts(ts = extract.months.from.ts(TREFHT_ZH, months=1), years = 2001:3000)

## Assess how well this works:
plot(x = as.numeric(TREFHT_ZH_hat1), y = as.numeric(TREFHT_ZH1), xlab = "Predicted", ylab = "Observed")
cor(as.numeric(TREFHT_ZH_hat1), as.numeric(TREFHT_ZH1)) ^ 2




# ------------------------------------------------------------------------------------------ 
# 2. TRAIN MODEL FOR AN ENTIRE TEMPERATURE FIELD:
# ------------------------------------------------------------------------------------------ 

library(foreach)
library(doParallel)
registerDoParallel(cores = 36)

CEU.TREFHT.extent = extent(c(-10, 30, 40, 60))
CEU.PSL.extent = extent(c(-30, 50, 20, 80))
PSL_CEU = crop(PSL_EUROPE, CEU.PSL.extent)
TREFHT_CEU = crop(TREFHT_EUROPE, CEU.TREFHT.extent)


# run dynamical adjustment on the whole field:
TREFHT_CEU_hat = train.dyn.adj.elastnet.annual.RB(Y.RB = TREFHT_CEU, X.RB = PSL_CEU, train.years = 1000:2000, train.months = 1, add.mon = 1, alpha = 0.2, nfolds = 10, s = "lambda.1se", 
                                              lags.to.aggregate = list(1), n.pc = 10, nr.cores = 36, x.domain = 20, y.domain = 20)

TREFHT_CEU_hat1 = extract.period.RB(RB = TREFHT_CEU_hat, years = 2001:3000, months = 1)
TREFHT_CEU1 = extract.period.RB(RB = TREFHT_CEU, years = 2001:3000, months = 1)


test.cor = gridcorts(rasterstack = brick(list(TREFHT_CEU_hat1, TREFHT_CEU1)), method = "pearson", type = "corel")


## Plot correlation between prediction and true values:
library(rworldmap)

plot(test.cor)
lines(coastsCoarse)














