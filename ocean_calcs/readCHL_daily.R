readCHL_daily <- function(date,time.resolution = "daily", xylim = NULL, ..., inputfiles = NULL, latest = TRUE,lon180 = FALSE,
                          nobsonly = FALSE,
                          returnfiles = FALSE) {
  if (is.null(inputfiles)) {
    ## memoize this call
    files <- ocfiles("daily", product = "MODISA", varnam = "CHL", type = "L3m")
  } else {
    files <- inputfiles
  }
  if (returnfiles)
    return(files)
  if (missing(date)) {
    if (latest)  date <- max(files$date) else date <- min(files$date)
  }
  date <- raadtools:::timedateFrom(date)
  files <- .processFiles(date, files, "daily")
  rl <- lapply(files$fullname, raster::raster, varname = "chlor_a")
  if (!is.null(xylim)) rl <- lapply(rl, raster::crop, raster::extent(xylim))
  raster::setZ(raster::brick(rl), date)
}

readCHL_monthly <- function(date,time.resolution = "monthly", xylim = NULL, ..., inputfiles = NULL, latest = TRUE,lon180 = FALSE,
                          nobsonly = FALSE,
                          returnfiles = FALSE) {
  if (is.null(inputfiles)) {
    ## memoize this call
    files <- ocfiles("monthly", product = "MODISA", varnam = "CHL", type = "L3m")
  } else {
    files <- inputfiles
  }
  if (returnfiles)
    return(files)
  if (missing(date)) {
    if (latest)  date <- max(files$date) else date <- min(files$date)
  }
  date <- raadtools:::timedateFrom(date)
  files <- .processFiles(date, files, "monthly")
  rl <- lapply(files$fullname, raster::raster, varname = "chlor_a")
  if (!is.null(xylim)) rl <- lapply(rl, raster::crop, raster::extent(xylim))
  raster::setZ(raster::brick(rl), date)
}

readCHL_weekly <- function(date,time.resolution = "weekly", xylim = NULL, ..., inputfiles = NULL, latest = TRUE,lon180 = FALSE,
                          nobsonly = FALSE,
                          returnfiles = FALSE) {
  if (is.null(inputfiles)) {
    ## memoize this call
    files <- ocfiles(time.resolution = "weekly", product = "MODISA", varnam = "CHL", type = "L3m")
  } else {
    files <- inputfiles
  }
  if (returnfiles)
    return(files)
  if (missing(date)) {
    if (latest)  date <- max(files$date) else date <- min(files$date)
  }
  date <- raadtools:::timedateFrom(date)
  files <- .processFiles(date, files, "weekly")
  rl <- lapply(files$fullname, raster::raster, varname = "chlor_a")
  if (!is.null(xylim)) rl <- lapply(rl, raster::crop, raster::extent(xylim))
  raster::setZ(raster::brick(rl), date)
}
