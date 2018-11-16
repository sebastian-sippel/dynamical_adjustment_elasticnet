# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for ANALYSIS 
# based on time series (xts) package in R for clean data handling
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018

get.year.from.ts <- function(ts) return(as.numeric(format(as.Date(time(ts)), "%Y")))
get.month.from.ts <- function(ts) return(as.numeric(format(as.Date(time(ts)), "%m")))


extract.months.from.ts <- function(ts, months=1) {
  return(ts[which(get.month.from.ts(ts) %in% months)])
}


extract.years.from.ts <- function(ts, years) {
  return(ts[which(get.year.from.ts(ts) %in% years)])
}


period.mean.ts <- function(ts, years, months = "annual", FUN=mean, na.rm = T) {
  # 1. extract relevant period from time series:
  period.ts = extract.years.from.ts(ts = ts, years = years)
  # 2. extract relevant months:
  if (months[1] == "annual") {
    return(apply.yearly(x = period.ts, FUN = FUN, na.rm=na.rm))
  } else {
    period.ts.monthly = apply.monthly(period.ts, FUN=mean)
    if (months[1] == "DJF" | months[1] == "NDJFM" | months[1] == "ONDJFMA") {
      if (months[1] == "DJF") cur.months = c(12, 1, 2)
      if (months[1] == "NDJFM") cur.months = c(11, 12, 1, 2, 3)
      if (months[1] == "ONDJFMA") cur.months = c(10, 11, 12, 1, 2, 3, 4)
      cur.ts = extract.months.from.ts(ts = period.ts.monthly, months = cur.months)
      ann.mean.ts = apply.yearly(xts(as.numeric(cur.ts), order.by = as.Date(time(cur.ts)) + 180), FUN=mean, na.rm=na.rm) # change time stamp to +6months and apply.yearly
      ret.ts = xts(as.numeric(ann.mean.ts), order.by = as.Date(time(ann.mean.ts))-180) # # change time stemp backwards to -6months
      return(ret.ts)
    } else {
      cur.ts = extract.months.from.ts(ts = period.ts.monthly, months = months)
      return(apply.yearly(cur.ts, FUN=mean, na.rm = na.rm))
    }
  }
}


R2.ts <- function(ts1, ts2, years, months = "annual") {
  return(cor(as.numeric(period.mean.ts(ts1, years=years, months=months)), as.numeric(period.mean.ts(ts2, years=years, months=months)), use = "complete.obs")^2)
}

require(hydroGOF)
mse.ts <- function(ts1, ts2, years, months = "annual") {
  return(hydroGOF::mse(as.numeric(period.mean.ts(ts1, years=years, months=months)), as.numeric(period.mean.ts(ts2, years=years, months=months))))
}
