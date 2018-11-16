# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for ANALYSIS 
# based on time series (xts) package in R for clean data handling
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018


### get period mean for raster field...:
# EOBS_TGhat_anom_LDT

## WRITE DIFFERENT FUNCTIONS TO OPERATE ON RB:

# RB.to.monthly <- function(RB) {
#  RB.values = values(RB)
#  RB.dates = as.Date(substring(names(RB), first = 2, last = 11), "%Y.%m.%d")
#  RB.values.monthly = apply(X = RB.values, MARGIN = 1, FUN=function(x) {
#    return(apply.monthly(x = xts(x, RB.dates), FUN = mean, na.rm=T))
#  })
# }





# Function to calculate area-weighted mean on Rasterbrick:
# --------------------------------------------------------
fldmean.RB <- function(RB, w = "area", mask = NULL, maskvalue = NA, ret.xts = T) {
  
  # check for potential mask:
  if (!is.null(mask)) RB = mask(RB, mask = mask, maskvalue = maskvalue)
  
  RB.values = values(RB)
  area.vec = values(area(RB))
  
  if (w == "area") {
    RB.mean.ts = apply(X = RB.values, MARGIN = 2, FUN=function(x) {
      weighted.mean(x, w = area.vec, na.rm = T)
    })  
  } else if (w == "none") {
    RB.mean.ts = cellStats(RB, stat="mean", na.rm = T)
  }
  
  if (ret.xts == T) {
    ret.ts = xts(x = RB.mean.ts, order.by = as.Date(substring(names(RB), first = 2, last = 11), "%Y.%m.%d"))
  } else if (ret.xts == F) {
    ret.ts = RB.mean.ts
  }

  return(ret.ts)
}



# Get period mean from Rasterbrick:
# --------------------------------------------------------
# RB = EOBS_TG_anom
# years = 1960:2017
# months = c(12, 1, 2)
# RB.values = values(RB)


extract.period.RB <- function(RB, years, months) {
  
  date.seq = as.Date(substring(text = names(RB), first = 2, last = 11), "%Y.%m.%d")
  time.idx = which(as.numeric(format(date.seq, "%Y")) %in% years & as.numeric(format(date.seq, "%m")) %in% months)

  return(subset(RB, time.idx))
}


period.mean.RB <- function(RB, months = "annual") {
  
  date.seq = as.Date(substring(text = names(RB), first = 2, last = 11), "%Y.%m.%d")
  years.seq = as.numeric(format(date.seq, "%Y"))
  months.seq = as.numeric(format(date.seq, "%m"))

  if (months == "annual") {
    ret.RB = brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-07-01", sep=""))
    return(ret.RB)
}
  
  # extract.ts.from.RB(TGhat_LDT_pred20CR.ERAI_ridge_pres.sfc)
}

  
  
  
  
  
  
  
  

