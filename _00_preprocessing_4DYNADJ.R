# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for preprocessing 
# based on time series (xts) package in R for clean data handling
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018

require(xts)
require(zoo)



## DESEASONALIZE TIME SERIES
#' @title Deseasonalize a time series
#' @description This is a generic function to obtain a deseasonalized time series based on monthly mean values
#' @param ts Time series to be deseasonalized
#' @param ref.year Years to use as a reference period
#' @param filt Should the monthly seasonal cycle be smoothed?
#' @param ret.seas.comp Should the seasonal component be 
#' @details
#' Some details about the function
#' @return A deseasonalized (daily or monthly) time series in xts format.
#' @references None
#' @author Sebastian Sippel
#' @examples
#' x = xts(x = rnorm(10^4), order.by = as.Date(1:10000, origin="1970-01-01"))
#' x.deseas = deseasonalize.ts(ts = x, ref.year = 1971:2000, filt = 1, ret.seas.comp = F, na.rm = T)
deseasonalize.ts <- function(ts, ref.year = 1971:2000, filt = 1, ret.seas.comp = F, na.rm = T) {
  
  # Convert time series to xts object:
  if (!is.xts(ts)) {
    xts.ts = as.xts(ts)
  } else {
    xts.ts = ts
  }
  
  # Return iff all values are NA
  if (all(is.na(xts.ts))) return(xts.ts)
  
  # Get monthly means via zoo::aggregate:
  
  # aggregate ref time series to monthly:
  ref.ts.filt = apply.monthly(x = window(xts.ts, start = (paste(head(ref.year, 1), "-01-01", sep="")), 
                                         end = (paste(tail(ref.year, 1), "-12-31", sep=""))), FUN=mean, na.rm = na.rm)
  ref.ts.filt = rollapply(data = ref.ts.filt, width = filt, FUN = mean, na.rm=T)
  ref.ts.monmean = aggregate(x = ref.ts.filt, by = format(as.Date(time(as.xts(ref.ts.filt))), "%m"), FUN=mean, na.rm = na.rm)
  
  # make new time series that contains monthly values:
  xts.ts.seasonal = xts(x = as.numeric(ref.ts.monmean)[as.numeric(format(as.Date(time(xts.ts)), "%m"))], order.by = index(xts.ts))
  
  if (ret.seas.comp == T) {
    return(xts.ts.seasonal)
  } else {
    return(xts.ts - xts.ts.seasonal)
  }
}


## SUBTRACT TREND VIA LOESS SMOOTHER
#' @title Subtract smooth monthly trend from a time series using LOESS/LOWESS smoother
#' @description Generic function to subtract smooth LOESS/LOWESS trend from a time series
#' @param ts Time series to be deseasonalized
#' @param freq Is the time series to be detrended in daily or monthly resolution?
#' @param filt Should the trend to be subtracted from each month be smoothed?
#' @param alpha Environment alpha from LOESS smoother 
#' @param degree Degree of polynomial for loess smoother
#' @param ret.trend.comp Should the trend component be returned?
#' @details
#' Some details about the function
#' @return A deseasonalized (daily or monthly) time series in xts format.
#' @references None
#' @author Sebastian Sippel
#' @examples
#' x = xts(x = rnorm(10^4), order.by = as.Date(1:10000, origin="1970-01-01"))
#' x.trend.comp = subtract.loess.ts(ts = x, freq = "daily", filt = 3, ret.trend.comp = T)
#' x.notrend = subtract.loess.ts(ts = x, freq = "daily", filt = 3, ret.trend.comp = F)
subtract.loess.ts <- function(ts, freq = "daily", filt = 3, alpha = 0.75, degree = 1, fill.NA = T, ret.trend.comp = F) {
  
  if (!is.xts(ts)) {
    xts.ts = as.xts(ts)
  } else {
    xts.ts = ts
  } 
  
  if (all(is.na(xts.ts))) return(xts.ts)
  
  # i. aggregate time series to monthly data (if it is not monthly yet):
  if (freq == "daily") {
    xts.ts.monmean = apply.monthly(x = xts.ts, FUN=mean, na.rm=fill.NA)
  } else if (freq == "monthly") {
    xts.ts.monmean = xts.ts
  }
  
  # ii. filter time series as desired for estimation of loess raster:
  # ref.ts.filt = rollmean(x = xts.ts.monmean, k = filt, fill = NA)
  ref.ts.filt = rollapply(data = xts.ts.monmean, width = filt, FUN = mean, na.rm=fill.NA)
  
  # iii. split into monthly value run through every month and calculate loess trend:
  f.split = as.factor(format(as.Date(time(as.xts(ref.ts.filt))), "%m"))
  mon.split = split(x = as.numeric(ref.ts.filt), f = f.split)
  loess.fitted = lapply(X = mon.split, FUN=function(x, alpha, degree, fill.NA) {
    if (fill.NA == F) {
      loess.fit = loess(c(x) ~ c(1:length(x)), span = alpha, degree = degree)
    } else if (fill.NA == T) {
      loess.fit = loess(c(x) ~ c(1:length(x)), span = alpha, degree = degree, control=loess.control(surface="direct"))
    }
    return(predict(object = loess.fit, newdata = c(1:length(x))))
  }, alpha = alpha, degree = degree, fill.NA = fill.NA)
  
  xts.ts.loess = xts(unsplit(value = loess.fitted, f = f.split), order.by = time(xts.ts.monmean))
  
  if (freq == "daily") {  # disaggregate to daily... 
    xts.ts.loess = disaggregate.monthly.to.daily.ts(monthly.ts = xts.ts.loess)
  } 
  
  # Return respective time series:
  if (ret.trend.comp == T) {
    return(xts.ts.loess)
  } else {
    return(xts.ts - xts.ts.loess)
  }
}





## ---------------------------------------------------------------------------------------- ##
## Other convenience functions:
## ---------------------------------------------------------------------------------------- ##


disaggregate.monthly.to.daily.ts <- function(monthly.ts) {
  ix.monthly=format(as.Date(time(monthly.ts)), "%Y-%m")
  date.seq=seq.Date(as.Date(paste(min(format(as.Date(time(monthly.ts)), "%Y-%m")), "-01", sep="")),max(as.Date(time(monthly.ts))), 1)
  ix.daily=format(date.seq, "%Y-%m")
  disagg.idx= match(x = ix.daily, table = ix.monthly)
  
  ret.ts = xts(x = monthly.ts[disagg.idx], order.by = date.seq)
  return(ret.ts)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}








