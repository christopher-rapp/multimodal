#' @title readPSD_BMI
#' Read BMI-format particle size distribution files
#'
#' @author Christopher Rapp
#'
#' @description
#' Recursively scans a directory for BMI instrument output files (MONO_DATA,
#' SEMS_CONC, SEMS_AUX, SEMS_RAW, SEMS_VOLTS naming patterns) and reads all
#' detected \code{SEMS_CONC} files into a list of data.frames. Column names
#' are stripped of their \code{Bin###_} prefix, a \code{UTC Time} (and, if
#' applicable, \code{Local Time}) column is constructed from the file's date
#' and time fields, and a \code{Size Range} column is added reflecting the
#' span of binned diameters present in that file.
#'
#' @param filepath Character string giving the path to a directory to
#'   search recursively (\code{list.files(..., recursive = TRUE)}) for BMI
#'   instrument files.
#' @param tz Character string giving the time zone of the timestamps in the
#'   files, e.g. \code{"UTC"} or an Olson name such as
#'   \code{"America/Los_Angeles"}. If \code{"UTC"}, timestamps are read
#'   directly into a \code{UTC Time} column. Otherwise, timestamps are first
#'   read into a \code{Local Time} column at the given time zone and then
#'   converted to a separate \code{UTC Time} column via
#'   \code{lubridate::with_tz()}.
#'
#' @returns A list of \code{data.frame}s, one per non-empty \code{SEMS_CONC}
#'   file found, each with: a \code{Units} column (\code{"dNdlogDp"}), a
#'   \code{UTC Time} column (and \code{Local Time} if \code{tz != "UTC"}), a
#'   \code{Size Range} column (max minus min binned diameter in that file),
#'   and columns reordered so the detected \code{POSIXt} time column comes
#'   first. Files that produce zero-row results are dropped from the
#'   returned list.
#'
#' @export
readPSD_BMI <- function(filepath, tz){

  # -------------------------------------------------------------------------- #
  ##### SECTION: File Wrangling #####
  #'

  {
    # List all level0 csv files
    files.bmi <- list.files(path = filepath,
                            recursive = T,
                            full.names = T)

    # Detect files of specific type
    files.SEMS_CONC <- files.bmi[grepl(paste("SEMS_CONC", collapse = '|'), ignore.case = T, files.bmi)]
  }

  # ------------------------------------------------------------------------ #
  ##### SECTION: Create Export Structure #####
  #'

  data.ls <- lapply(files.SEMS_CONC, function(x) {

    # Use data.table's fread to read in data
    # Fastest method in reading CSV's and least memory consuming
    # fill MUST equal FALSE with SPIN files
    tmp.df <- data.table::fread(paste0(x),
                                na.strings = c("", "NA", "NaN"),
                                skip = "#StartDate",
                                fill = TRUE,
                                strip.white = TRUE,
                                stringsAsFactors = FALSE)

    if (length(str_which(colnames(tmp.df), "Bin\\d{1,3}_")) != 0){
      setnames(tmp.df, new = str_replace(colnames(tmp.df), "Bin\\d{1,3}_", ""))
    } else {
      tmp.df <- NULL
    }

    if (nrow(tmp.df) >= 1){

      # Note the units for the binned data
      tmp.df <- tmp.df %>%
        mutate("Units" = "dNdlogDp")

      if (tz == "UTC") {

        # Create time objects
        tmp.df$`UTC Time` = as.POSIXct(paste0(tmp.df$`#StartDate`, ' ', tmp.df$`StartTime`),
                                       format = '%y%m%d %H:%M:%S', tz = tz)

      } else {

        # Create time objects
        tmp.df$`Local Time` = as.POSIXct(paste0(tmp.df$`#StartDate`, ' ', tmp.df$`StartTime`),
                                         format = '%y%m%d %H:%M:%S', tz = tz)

        # Create time objects
        tmp.df$`UTC Time` = with_tz(tmp.df$`Local Time`, tzone = "UTC")
      }

      # Retrieve midpoint diameters and associated columns
      bins.ix = str_which(colnames(tmp.df), "\\d+\\.\\d+")
      bins.nm = as.numeric(colnames(tmp.df)[bins.ix])

      # Calculate the size range within the file
      tmp.df <- tmp.df %>%
        mutate(`Size Range` = bins.nm[length(bins.nm)] - bins.nm[1])
    }

    # Identify column classes of dataframe for usage later
    col.classes <- lapply(tmp.df, class)

    # Find columns containing the POSIX class
    time.index = c(grep("POSIXt", col.classes))

    # Structure so time column is in the front
    tmp.df <- tmp.df %>%
      dplyr::select(all_of(time.index), everything())
  })

  # Remove empty data.frames
  data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)

  return(data.ls)
}
