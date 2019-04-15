# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for preprocessing 
# based on time series (xts) package in R for clean data handling
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018

require(xts)
require(zoo)



# Convert raster name to time series:
extract.ts.from.RB <- function(RB, x, y) {
  
  test=c(extract(RB, y = data.frame(x = x, y = y)))
  names(test) = names(RB)
  
  cur.dates = as.Date(substring(text = names(test), first = 2, last = 11), format = "%Y.%m.%d")
  return(xts(test, order.by = cur.dates))
}


# CONVERT RB TO XTS for use in calc function (!!)
convert.RBts.to.ts <- function(x) {
  cur.dates = as.Date(substring(text = names(x), first = 2, last = 11), format = "%Y.%m.%d")
  return(xts(x, order.by = cur.dates))
}


## SPACE-TIME COMBINATION convenience function:
subset.dates.from.RB <- function(cr, date.seq) {
  cr.date  = as.Date(substring(names(cr), first = 2), format = "%Y.%m.%d")
  match.idx = match(x = date.seq, table = cr.date)
  return(subset(cr, na.omit(match.idx)))
}

## 
subset.years.from.RB <- function(cr, years) {
  cr.years  = as.numeric(substring(names(cr), first = 2, last = 5))
  # match.idx = match(x = years, table = cr.years)
  match.idx = which(cr.years %in% years)
  return(subset(cr, na.omit(match.idx)))
}



# set all grid cells NA that contain more than X% NA's:
rm.NA.from.RB <- function(RB, NA.frac = 0.1) {
  
  test=calc(RB, fun=function(x) (!all(is.na(x)) & ((length(which(is.na(x))) / length(x)) > NA.frac) ))
  values(RB)[which(values(test) == 1),]=NA   # CHECK if there is any not-NA value: any(!is.na(values(RB)[which(values(test) == 1),]))
  return(RB)  
}
  
  


## EOF detrend RB object:
detrend.EOF.RB <- function(cr, freq = "daily", filt = 3, alpha = 1, degree = 1, n.EOF.detrend = 20) {
  
  # 1. estimate area and define indices:
  area.vec = values(raster::area(cr))
  cr.val=values(cr) * sqrt(area.vec)
  time.idx=as.Date(substring(names(cr), 2), format="%Y.%m.%d")
  
  # 2. SVD transform:
  PSL.SVD.transform = svd(x = cr.val)
  print("SVD transform done!")
  
  # 3. Define how many EOFs are to be detrended:
  # n.EOF.detrend = which(cumsum( PSL.SVD.transform$d^2 / sum(PSL.SVD.transform$d^2)) > prop.var)[1]
  # if (prop.var == 1) n.EOF.detrend = 20
  # n.EOF.detrend = 20
  
  # 4. Project ALL VALUES on EOFs:
  Xproj.svd_ALL = t(PSL.SVD.transform$u) %*% (cr.val)  # projection of FULL time period onto EOF's 
  print("1st matrix multiplication done!")
  # works because -> t(U) = U^-1 for orthonormal matrix; and A = U x D x V; thus t(U) x A = D x V
  
  # 5. Detrend each EOF with loess smoother:
  Xproj.svd_ALL_detrended = t(sapply(X = 1:(n.EOF.detrend), FUN=function(pc.idx) {
    cur.daily.ts = xts(Xproj.svd_ALL[pc.idx,], order.by = time.idx)   # for example: pc.idx = 11
    cur.loess.fit.daily = subtract.loess.ts(ts = cur.daily.ts, freq = freq, filt = filt, alpha = alpha, degree = degree, ret.trend.comp = T)
    # plot(cur.loess.fit.daily)
    res.ts.daily = cur.daily.ts - cur.loess.fit.daily
    return(res.ts.daily)
  }))
  
  # 6. fill with low-variance EOFs:     
  # if (prop.var == 1) {
  #  Xproj.svd_ALL_detrended_filled = Xproj.svd_ALL_detrended
  #} else {
  Xproj.svd_ALL_detrended_filled = rbind(Xproj.svd_ALL_detrended, Xproj.svd_ALL[(n.EOF.detrend+1):(dim(Xproj.svd_ALL)[1]),])
  #}
  
  # 7. project back onto full spatial pattern (and transform back to raster):
  PSL_reconstr = (PSL.SVD.transform$u %*% Xproj.svd_ALL_detrended_filled) / sqrt(area.vec)
  print("2nd matrix multiplication done!")
  psl.raster_reconstr = cr
  values(psl.raster_reconstr) = PSL_reconstr
  return(psl.raster_reconstr)
}




extract.date.RB <- function(RB) {
  return(as.Date(substring(text = names(RB), first = 2, last = 11), format = "%Y.%m.%d"))
}


# -----------------------------------------------------------
## DETREND RB in parallel:
# -----------------------------------------------------------
require(foreach)
require(doParallel)
registerDoParallel(cores = 35)

subtract.loess.ts.RB <- function(RB, nr.cores=1,
                                 freq = "daily", filt = 3, alpha = 1, degree = 1) {

  registerDoParallel(cores = nr.cores)
  
  RB.values=t(values(RB))
  RB.date  = as.Date(substring(names(RB), first = 2, last = 11), format = "%Y.%m.%d")
  
  RB.values.notrend=foreach(i=1:(dim(RB.values)[2])) %dopar% {
  # RB.values.notrend=foreach(i=1:(100)) %dopar% {
    print(i)
    ts=xts(RB.values[,i], order.by = RB.date)
    ts.notrend=subtract.loess.ts(ts=ts, freq = freq, filt = filt, alpha = alpha, degree = degree)
    return(as.numeric(ts.notrend))
  }
  # convert list to matrix:
  RB.values.new=sapply(RB.values.notrend, FUN=function(x) x)

  values(RB) = t(RB.values.new)
  return(RB)
}


# DESEAONALIZE TIME SERIES IN PARALLEL:




