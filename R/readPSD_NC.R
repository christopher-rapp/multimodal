#' @title readPSD_NC
#' Read NetCDF-format particle size distribution files
#'
#' @author Christopher Rapp
#'
#' @description
#' Reads one or more NetCDF (\code{.nc}) files containing particle size
#' distribution data and parses each into a data.frame. Timestamps are
#' decoded from the file's \code{time} variable and units attribute using
#' the CF time convention, binned mobility diameters are recovered from
#' \code{diameter_mobility_bounds} as bin midpoints, and quality-control
#' flags in \code{qc_dN_dlogDp} are applied to zero out flagged
#' concentration values in \code{dN_dlogDp}. Columns that are entirely
#' \code{NA} after masking are dropped.
#'
#' @param filepath Character string giving either the path to a single
#'   NetCDF file, or the path to a directory containing one or more such
#'   files. If \code{filepath} is a directory, all files immediately inside
#'   it are read (via \code{list.files(filepath, full.names = TRUE)}); if
#'   no files are found there, \code{filepath} itself is treated as the file
#'   to read.
#'
#' @returns A list of \code{data.frame}s, one per file read, each with a
#'   \code{Time} column (\code{POSIXct}, decoded via \code{CFtime}) followed
#'   by one column per binned mobility diameter (column names are the
#'   rounded bin midpoint diameters), containing QC-masked
#'   \code{dN_dlogDp} concentration values.
#'
#' @export
readPSD_NC <- function(filepath){

  # List files if multiple
  files <- list.files(filepath, full.names = T)

  if (length(files) == 0){
    files.nc = filepath
  }

  if (length(files) >= 1){
    files.nc <- files
  }

  data.nc <- lapply(files.nc, function(x){

    tmp.nc <- nc_open(x)
    on.exit(nc_close(tmp.nc), add = TRUE)

    # Time variable
    {
      time <- ncvar_get(tmp.nc, "time")
      time.unit <- ncatt_get(tmp.nc, "time", "units")

      # Use CFtime package to retrieve timestamp
      time <- CFtime::CFtime(time.unit$value, calendar = "standard", time) # convert time to CFtime class

      # Convert to POSIX format
      time <- lubridate::as_datetime(time$as_timestamp())
    }

    # Concentration values
    {
      conc <- t(ncvar_get(tmp.nc, 'dN_dlogDp'))
      conc.qc <- t(ncvar_get(tmp.nc, 'qc_dN_dlogDp'))
    }

    # Mobility diameters
    {
      bins <- ncvar_get(tmp.nc, 'diameter_mobility_bounds')
      bins = round(rowMeans(t(bins)), 2)
    }

    # Create dataframe
    {
      tmp.df <- as.data.frame(conc)

      # Setnames as diameter midpoints
      colnames(tmp.df) <- round(bins, 2)

      # Apply QC from conc.qc
      tmp.df <- tmp.df * as.data.frame(!conc.qc)

      # Remove empty columns
      tmp.df <- purrr::discard(tmp.df, ~all(is.na(.)))

      # Add time column to dataframe
      tmp.df <- cbind("Time" = time, tmp.df)
    }
  })

  return(data.nc)
}




