#' @title readPSD_NAS.R
#' Read NASA-AMES format particle size distribution files
#'
#' @author Christopher Rapp
#'
#' @description
#' Reads one or more files in the standard NASA-AMES (1001) format, as
#' commonly used for ambient/airborne aerosol measurements, and parses each
#' into a data.frame. Metadata lines preceding the variable header are used
#' to recover the file's reference date and to relabel columns (start/end
#' time, binned diameters, pressure, temperature). Only columns representing
#' the "mean" of a measurement (plus the time and flag columns) are
#' retained. Time columns, stored as fractional days since a reference date
#' per the NASA-AMES convention, are converted to \code{POSIXct}.
#'
#' @param filepath Character string giving either the path to a single
#'   NASA-AMES file, or the path to a directory containing one or more such
#'   files. If \code{filepath} is a directory, all files immediately inside
#'   it are read (via \code{list.files(filepath, full.names = TRUE)}); if
#'   no files are found there, \code{filepath} itself is treated as the file
#'   to read.
#'
#' @returns A list of \code{data.frame}s, one per file read, each containing
#'   a \code{starttime} and \code{endtime} column (converted to
#'   \code{POSIXct}), the retained "mean" measurement columns (binned
#'   diameter columns renamed to their numeric diameter value, pressure and
#'   temperature columns relabeled with detected units where present), and
#'   the trailing flag column.
#'
#' @export
readPSD_NAS <- function(filepath){

  # List files if multiple
  files <- list.files(filepath, full.names = T)

  if (length(files) == 0){
    files.nas = filepath
  }

  if (length(files) >= 1){
    files.nas <- files
  }

  data.nas <- lapply(files.nas, function(x){

    # Read file
    tmp.nas <- read.delim(x, header = F)

    # Identify where the variable headers are located
    header.ix = as.numeric(stringr::str_extract(tmp.nas[1, 1], "(\\d+)\\s\\d+", group = T))

    # Extract
    header.c <- tmp.nas[header.ix, ]

    # Retrieve metadata
    meta.data <- as.vector(tmp.nas[1:(as.numeric(header.ix) - 1), 1])

    # Create a dataframe
    data <- tmp.nas[as.numeric(header.ix):nrow(tmp.nas), ]
    data <- read.table(text = data, header = T)

    # Identify standard variable labels
    var1 <- stringr::str_which(meta.data, "end_time")
    var2 <- stringr::str_which(meta.data, "numflag")

    # Variable names
    variables.c <- tmp.nas[var1:var2, ]

    # Time variables
    # Based on standard NASA-AMES format
    {
      date.range <- meta.data[7]
      time.format <- meta.data[9]

      # Regular expressions to catch dates
      date1 <- stringr::str_extract(date.range, "\\d{4}\\s{1}\\d{2}\\s{1}\\d{2}")
      date2 <- stringr::str_extract(date.range, "\\s{1}(\\d{4}\\s{1}\\d{2}\\s{1}\\d{2})", group = T)

      # Format to UTC
      date1 <- as.POSIXct(date1, format = "%Y %m %d", tz = "UTC")
      date2 <- as.POSIXct(date2, format = "%Y %m %d", tz = "UTC")
    }

    if (length(variables.c) != ncol(data)){
      variables.c <- append(paste0("start_time of measurement, ", time.format), variables.c)
    }

    # For PSD's retrieve the mean
    {
      keep.c <- stringr::str_which(variables.c, "mean")

      # This keeps the two time columns
      keep.c <- append(1:2, keep.c)

      # This keeps the flag variable
      keep.c <- append(keep.c, length(variables.c))

      # Assumes data.table format
      data <- data[, keep.c]

      # Which variables are labeled with mean
      tmp.c <- variables.c[keep.c]
    }

    # Binned data
    {
      # Use regex
      bin.ix <- stringr::str_which(tmp.c, "D=\\d+")
      binned.c <- tmp.c[bin.ix]

      # Extract bins
      binned.c <- as.numeric(stringr::str_extract(binned.c, "(\\d*\\.\\d*)"))

      # Setnames as the bins
      tmp.c[bin.ix] <- binned.c
    }

    # Temperature and pressure if present
    {
      # Pressure
      if (any(stringr::str_which(tmp.c, "hPa|pressure|Pressure|pres."))){

        # Rename
        tmp.c[stringr::str_which(tmp.c, "hPa|pressure|Pressure|pres.")] <- "pressure (hPa)"
      }

      # Temperature
      if (any(stringr::str_which(tmp.c, "temperature|Temperature|Temp."))){

        # Rename
        tmp.c[stringr::str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature"

        # If temperature exceeds 70, it is in Kelvin
        if (mean(data[, stringr::str_which(tmp.c, "temperature|Temperature|Temp.")]) > 70){
          tmp.c[stringr::str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (K)"
        } else {
          tmp.c[stringr::str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (C)"
        }
      }
    }

    # Rename time variables
    tmp.c[1] <- "starttime"
    tmp.c[2] <- "endtime"

    # Rename dataframe
    data.table::setnames(data, new = tmp.c)

    # Format time to POSIX
    {
      data$starttime <- as.POSIXct(date1) + as.difftime(data$starttime, units = "days")
      data$endtime <- as.POSIXct(date1) + as.difftime(data$endtime, units = "days")
    }

    return(data)
  })

  return(data.nas)
}