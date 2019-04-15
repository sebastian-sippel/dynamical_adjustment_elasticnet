# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for Prediction step in dynamical adjustment
# Based on the Elastic Net
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018

require(xts)
require(zoo)
require(raster)
require(ncdf4)
require(foreach)
require(doParallel)
registerDoParallel(cores = 44)










#' @title Train a parallel, shared-memory (MEMORY-FRIENDLY!) model for a RASTERBRICK climatic target variable Y based on circulation pattern X
#' @description 
#' @param Y.RB Rasterbrick with consistent time format
#' @param X.RB 
#' @param X.RB.train
#' @param train.years subset of years to train on for dynamical adjustment (if NULL: all years are used)
#' @details
#' Some details about the function
#' @return A deseasonalized (daily or monthly) time series in xts format.
#' @references None
#' @author Sebastian Sippel
#' @examples
train.dyn.adj.elastnet.annual.RB <- function(Y.RB, X.RB, X.RB.train = NULL, 
                                                    train.years = NULL, train.months = 1:12, add.mon = 2,
                                                    alpha = 0.2, nfolds = 10, s = "lambda.min", 
                                                    lags.to.aggregate = list(1:2, 3:7, c(8:30)), n.pc = 10, nr.cores = 46) {
  
  # (1) Subset relevant years for Y in training sample and KICK out NA grid cells:
  # -----------------------------------------------
  if(is.null(train.years)) train.years = unique(as.numeric(format(as.Date(substring(names(Y.RB), first = 2, last = 11), format = "%Y.%m.%d"), "%Y")))
  Y.RB.values = t(values(subset.years.from.RB(cr = Y.RB, years = train.years)))                 # Get Y values (only!) for training time 
  Y.RB.date = as.Date(substring(names(Y.RB), first = 2, last = 11), format = "%Y.%m.%d")        # Get dates for training entie y time series
  Y.RB.date = Y.RB.date[which(as.numeric(format(Y.RB.date, "%Y")) %in% train.years)]            # Get only dates fpr training time series
  
  NONA.idx = which(apply(X = Y.RB.values, MARGIN=2, FUN=function(x) !all(is.na(x))))             # Get NA values in Y raster
  Y.RB.values.nona = Y.RB.values[,NONA.idx]                                                     # New Y dataset with NA's removed
  Y.RB.coord.nona = coordinates(Y.RB)[NONA.idx,]                                                # Coordinates of Y dataset.
  
  # (2) Sanity CHECKS and (possible) early return: 
  # ---------------------------------------------------------------------
  if (dim(Y.RB.values.nona)[2] == 0) {    # are all values NA in Y.RB.values.nona?
    Yhat.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, nlayers(X.RB)))
    values(Yhat.RB) = NA; names(Yhat.RB) = names(X.RB);  # ALL VALUES SHOULD BE NA !
    return(Yhat.RB)
  }
  # remove unnecessary variables to avoid exceeding memory:
  rm(Y.RB.values)

  ## (3) RUN IN PARALLEL Dynamical adjustment:
  ## ---------------------------------------------------------------------

  # 3.1 Definition of bigmemory matrices:
    Y.values.nona.bm <- as.big.matrix(x = Y.RB.values.nona, type = "double", separated = FALSE, shared = TRUE,
                            backingfile = "Y.bin", descriptorfile = "Y.desc")
  
    X.bm = as.big.matrix(x = t(values(X.RB)), type = "double", separated = FALSE, shared = TRUE,
                                   backingfile = "X.bin", descriptorfile = "X.desc")
    
    X.train.bm.bin = F
    if (!is.null(X.RB.train)) {
      X.train.bm.bin = T
      X.train.bm = as.big.matrix(x = t(values(X.RB.train)), type = "double", separated = FALSE, shared = TRUE,
                                 backingfile = "Xtrain.bin", descriptorfile = "Xtrain.desc")
    }
    
  # 3.2 Define small raster object to use in parallel loops:
  X.coords = coordinates(X.RB)
  X.template = raster(subset(X.RB, 1))
  values(X.template) <- 1:(dim(X.coords)[1])
  
  # 3.3 Set up cluster based on SNOW/SOCKET architecture (because of smaller memory load):
  ## ---------------------------------------------------------------------------------------
  cl <- parallel::makeCluster(nr.cores, outfile="")
  doParallel::registerDoParallel(cl)
  clusterExport(cl=cl, varlist = list("train.dyn.adj.elastnet.annual", "add.lag.predictors", "rep.row"))
  
  # 3.4. RUN DYNAMICAL ADJUSTMENT:
  ## ---------------------------------------------------------------------------------------
  
  # Start the clock!
  # ptm <- proc.time()
  RB.values.notrend = foreach(i=1:(dim(Y.RB.values.nona)[2]), 
                              .packages = c("bigmemory", "glmnet", "xts", "raster", "zoo")) %dopar% {
    print(i)
    
    # (a) Attach bigmemory matrices:
    Y <- bigmemory::attach.big.matrix("Y.desc")
    X <- bigmemory::attach.big.matrix("X.desc")
    if (X.train.bm.bin == T) {
      Xtrain <- bigmemory::attach.big.matrix("Xtrain.desc")
    } else if (X.train.bm.bin == F) {
      Xtrain = NULL
    }
    
    # (b) Subset/Define new raster based on present data point:
    cur.extent = extent(c(Y.RB.coord.nona[i,1] - 15, Y.RB.coord.nona[i,1] + 15, Y.RB.coord.nona[i,2] - 15, Y.RB.coord.nona[i,2] + 15))
    cur.X.idx = values(crop(X.template, cur.extent))
    cur.X.RB = X[,cur.X.idx]
    cur.Xtrain = NULL
    if (X.train.bm.bin == T) cur.Xtrain = Xtrain[,cur.X.idx]
    # X.date.seq = as.Date(substring(dimnames(X[,cur.X.idx])[[1]], first=2, last=11), format="%Y.%m.%d")
    # cur.X.RB = brick(crop(SLP.raster.template, cur.extent), values = F, nl= 24837)
    # values(cur.X.RB) = t(X[,cur.X.idx])
    
    # (c) Run Dynamical Adjustment:
    Yhat = tryCatch({
      # print(A) # produce error to see if this can work...
    Yhat = train.dyn.adj.elastnet.annual(Y.train = xts(Y[,i], order.by = Y.RB.date), X = cur.X.RB, X.train = cur.Xtrain, X.date.seq = NULL, 
                                         train.years = train.years, train.months = train.months, add.mon = add.mon, 
                                         alpha = alpha, nfolds = nfolds, ret.cv.model=F, s = s,
                                         lags.to.aggregate = lags.to.aggregate, n.pc = n.pc)
    }, error = function(err) {
      print("ERROR HAS BEEN CAUGHT")
      Yhat = rep(NA, dim(cur.X.RB)[1])
    })
    return(Yhat)
  }
  # Stop the clock
  # proc.time() - ptm
  parallel::stopCluster(cl = cl)
  
  RB.values.notrend_mat = sapply(X = RB.values.notrend, FUN=function(x) as.numeric(x))
  
  ## (3) Derive NEW Yhat raster and return:
  ## ---------------------------------------------------------
  Yhat.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, nlayers(X.RB)))
  values(Yhat.RB)[NONA.idx,] = t(RB.values.notrend_mat)     # Put predictions into Yhat raster raster (!!)
  names(Yhat.RB) = names(X.RB)
  
  # CONTINUE HERE WITH GOOD TRAINING+PREDICTION SOLUTION:
  # -> "train.dyn.adj.elastnet.annual" selects prediction years based on Y.ts but should be for entire "X" array
  # test.obs = t(values(Y.RB))[,NONA.idx[1]]
  # test.obs = xts(test.obs, as.Date(substring(names(Y.RB), first = 2, last = 11), format = "%Y.%m.%d"))
  #    cor(as.numeric(period.mean.ts(ts = test, years = 1970:2017, months = "DJF")),
  #         as.numeric(period.mean.ts(ts = test.obs, years = 1970:2017, months = "DJF")), use = "complete.obs") ^ 2
  
  return(Yhat.RB)
}



#### ancillary functions:

write.dailyRB.to.ncdf <- function(RB, out.file = "test.nc",
                                  var.name = "temperature", var.units = "degC",
                                  origin = "1950-01-01", add.abs.file = NULL, to.monthly = F, to.seasonal = F, to.annual = F, window.width = 5) {
  require(ncdf4)
  
  RB.values = values(RB)
  RB.coord = coordinates(RB)
  date.format = paste("days since ", origin, " 12:00", sep="")
  time_step_daily = as.numeric(as.Date(substring(names(RB), 2, 11), format="%Y.%m.%d")) + as.numeric(as.Date(0) - as.Date(0, origin= as.Date(origin)))
  
  # define dimensions
  londim <- ncdim_def("lon", "degrees_east", as.double(unique(RB.coord[,1])))
  latdim <- ncdim_def("lat", "degrees_north", as.double(unique(RB.coord[,2])))
  dim_time <- ncdim_def('time', date.format, as.double(time_step_daily), unlim=T) 
  
  # Create a new Variable:
  var_out <- ncvar_def(var.name, var.units, list(londim, latdim,dim_time), 9.e20) 
  ncid_out <- nc_create(out.file, var_out) 
  ncvar_put(ncid_out, var_out, RB.values, start=c(1, 1, 1), count=c(length(unique(RB.coord[,1])), length(unique(RB.coord[,2])), length(time_step_daily)))
  nc_close(ncid_out)
  
  system("module load cdo")
  if (!is.null(add.abs.file)) {
    out.file.abs = paste(substring(out.file, 1, nchar(out.file)-3), "_abs.nc", sep="")
    system(paste("cdo -ymonadd ", out.file, " ", add.abs.file, " ", out.file.abs, sep=""))
  }
  
  if (to.monthly == T) {
    # get monthly mean:
    system(paste("cdo -monmean ", out.file, " ", substring(out.file, 1, nchar(out.file)-3), "_monmean.nc", sep=""))
    # get monthly max:
    system(paste("cdo -monmax -runmean,",window.width, " " , out.file, " ", substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,".nc", sep=""))
    # get monthly min:
    system(paste("cdo -monmin -runmean,",window.width, " " , out.file, " ", substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,".nc", sep=""))
  }
  
  if (to.monthly == T & !is.null(add.abs.file)) {
    # get monthly mean:
    system(paste("cdo -monmean ", out.file.abs, " ", substring(out.file, 1, nchar(out.file)-3), "_monmean_abs.nc", sep=""))
    # get monthly max:
    system(paste("cdo -monmax -runmean,",window.width, " " , out.file.abs, " ", substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,"_abs.nc", sep=""))
    # get monthly min:
    system(paste("cdo -monmin -runmean,",window.width, " " , out.file.abs, " ", substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,"_abs.nc", sep=""))
  }
  
  if (to.seasonal == T) {
    # get monthly mean:
    system(paste("cdo -seasmean ", substring(out.file, 1, nchar(out.file)-3), "_monmean.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmean.nc", sep=""))
    # get monthly max:
    system(paste("cdo -seasmax ",  substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,".nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmax", window.width,".nc", sep=""))
    # get monthly min:
    system(paste("cdo -seasmin ",  substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,".nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmin", window.width,".nc", sep=""))
  }
  
  if (to.seasonal == T & !is.null(add.abs.file)) {
    # get monthly mean:
    system(paste("cdo -seasmean ", substring(out.file, 1, nchar(out.file)-3), "_monmean_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmean_abs.nc", sep=""))
    # get monthly max:
    system(paste("cdo -seasmax ",  substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,"_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmax", window.width,"_abs.nc", sep=""))
    # get monthly min:
    system(paste("cdo -seasmin ",  substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,"_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_seasmin", window.width,"_abs.nc", sep=""))
  }
  
  if (to.annual == T) {
    # get monthly mean:
    system(paste("cdo -yearmean ", substring(out.file, 1, nchar(out.file)-3), "_monmean.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmean.nc", sep=""))
    # get monthly max:
    system(paste("cdo -yearmax ",  substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,".nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmax", window.width,".nc", sep=""))
    # get monthly min:
    system(paste("cdo -yearmin ",  substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,".nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmin", window.width,".nc", sep=""))
  }
  
  if (to.annual == T & !is.null(add.abs.file)) {
    # get monthly mean:
    system(paste("cdo -yearmean ", substring(out.file, 1, nchar(out.file)-3), "_monmean_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmean_abs.nc", sep=""))
    # get monthly max:
    system(paste("cdo -yearmax ",  substring(out.file, 1, nchar(out.file)-3), "_monmax", window.width,"_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmax", window.width,"_abs.nc", sep=""))
    # get monthly min:
    system(paste("cdo -yearmin ",  substring(out.file, 1, nchar(out.file)-3), "_monmin", window.width,"_abs.nc", " ", substring(out.file, 1, nchar(out.file)-3), "_yearmin", window.width,"_abs.nc", sep=""))
  }
  
  return("NetCDF file succesfully written to Disk.")
}





shift.extent.360deg <- function(raster.object, add.360 = T) {
  cur.extent = as.matrix(extent(raster.object))
  
  # add360 -> shift from left to right!
  if (add.360 == T) {
    cur.extent[1,] = cur.extent[1,] + 360
    extent(raster.object) = cur.extent
  } else if (add.360 == F) {
    cur.extent[1,] = cur.extent[1,] - 360
    extent(raster.object) = cur.extent
  }
  return(raster.object)
}



