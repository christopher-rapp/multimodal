#' @title readPSD_TSI
#' Read TSI-format particle size distribution files
#'
#' @author Christopher Rapp
#'
#' @description
#' Recursively scans a directory for TSI instrument output CSV files and
#' reads each into a data.frame. Handles a known TSI export defect where the
#' \code{"Laser, Flow,"} status-flag label embeds extra unescaped commas,
#' which otherwise misaligns column parsing; when a row's status flag
#' contains \code{"flow"}/\code{"Flow"} instead of the expected numeric
#' \code{td(s)} value, that row's trailing columns are shifted left to
#' correct the offset. Column headers containing non-UTF8 characters (from
#' degree/micro/greater-than-or-equal symbols in unit labels) are cleaned.
#' A \code{UTC Time} (and, if applicable, \code{Local Time}) column is
#' constructed from the file's date and time fields, and a \code{Size
#' Range} column is added reflecting the span of binned diameters present
#' in that file.
#'
#' @param filepath Character string giving the path to a directory to
#'   search recursively (\code{list.files(..., recursive = TRUE)}) for TSI
#'   instrument CSV files.
#' @param tz Character string giving the time zone of the timestamps in the
#'   files, e.g. \code{"UTC"} or an Olson name such as
#'   \code{"America/Los_Angeles"}. If \code{"UTC"}, timestamps are read
#'   directly into a \code{UTC Time} column. Otherwise, timestamps are first
#'   read into a \code{Local Time} column at the given time zone and then
#'   converted to a separate \code{UTC Time} column via
#'   \code{lubridate::with_tz()}.
#'
#' @returns A list of \code{data.frame}s, one per non-empty file found, each
#'   with: a \code{Units} column (\code{"dNdlogDp"}), a \code{UTC Time}
#'   column (and \code{Local Time} if \code{tz != "UTC"}), a \code{Size
#'   Range} column (max minus min binned diameter in that file), and
#'   columns reordered so the detected \code{POSIXt} time column comes
#'   first.
#'
#' @export
readPSD_TSI <- function(filepath, tz){

  # -------------------------------------------------------------------------- #
  ##### SECTION: File Wrangling #####
  #'

  # List all level0 csv files
  files <- list.files(path = filepath,
                          recursive = T,
                          full.names = T)

  if (length(files) == 0){
    files.tsi = filepath
  }

  if (length(files) >= 1){
    files.tsi <- files
  }

  data.ls <- lapply(files.tsi, function(x) {

    # This is a uncommon issue with TSI instruments
    if (any(stringr::str_detect(readLines(paste0(x)), "Laser, Flow,")) == T){

      print(paste0("WARNING: ", "Unallowed delimiter detected in status flag for ", x, ". Removing commas and retrying..."))

      tmp.catch <- paste0(stringr::str_replace(readLines(paste0(x)), "Laser, Flow,", "Laser Flow"), collapse = "\n")

      tmp.df = data.table::fread(tmp.catch, header = TRUE, stringsAsFactors = FALSE,
                     fill = TRUE, skip = 'Sample #', check.names = FALSE)

      print("Done")
    } else {

      tmp.df <- data.table::fread(paste0(x), header = TRUE, stringsAsFactors = FALSE,
                      fill = TRUE, skip = 'Sample #', check.names = FALSE)
    }

    issues = NULL
    issues = stringr::str_which(tmp.df$`td(s)`, "flow")
    issues = append(issues, str_which(tmp.df$`td(s)`, "Flow"))

    if (length(issues) > 0){

      # This is used to correct the embedded delimiter issue with
      # Status Flag "CPC Laser, Flow, or Temp out of range causes this problem
      for (i in issues) {

        # Index two columns ahead and shift values to the left
        start = which(c(colnames(tmp.df)) == "td(s)") + 2
        end = ncol(tmp.df)

        for (j in start:end){

          tmp.df[i, (j - 2)] <- tmp.df[i, j, with = FALSE]
        }
      }

      rm(i, j, start, end)
    }

    # Extract string of column names
    tmp.nm = colnames(tmp.df)

    # Change any non-valid locale column headers
    if (any(!validUTF8(tmp.nm)) == TRUE){

      # Replace non-locale characters
      tmp.nm = stringr::str_replace(tmp.nm, '<', '')
      tmp.nm = stringr::str_replace(tmp.nm, 'cm�', 'cc')
      tmp.nm = stringr::str_replace(tmp.nm, '≥', '')
      tmp.nm = stringr::str_replace(tmp.nm, '\\?', '')
      tmp.nm = stringr::str_replace(tmp.nm, '\\. ', '.')
      tmp.nm = stringr::str_replace(tmp.nm, ' \\(', '(')

      # Rename columns without erroneous column characters
      setnames(tmp.df, tmp.nm)
    }

    tmp.df <- tmp.df %>%
      mutate("Units" = "dNdlogDp")

    if (tz == "UTC") {

      # Create time objects
      tmp.df$`UTC Time` = as.POSIXct(paste0(tmp.df$Date, ' ', tmp.df$`Start Time`),
                                     format = '%m/%d/%y %H:%M:%S', tz = tz)

    } else {

      # Create time objects
      tmp.df$`Local Time` = as.POSIXct(paste0(tmp.df$Date, ' ', tmp.df$`Start Time`),
                                       format = '%m/%d/%y %H:%M:%S', tz = tz)

      # Create time objects
      tmp.df$`UTC Time` = with_tz(tmp.df$`Local Time`, tzone = "UTC")
    }

    # Retrieve median diameters and associated columns
    bins.ix = stringr::str_which(colnames(tmp.df), "\\d+\\.\\d+")
    bins.nm = as.numeric(colnames(tmp.df)[bins.ix])

    # Calculate the size range within the file
    tmp.df <- tmp.df %>%
      mutate(`Size Range` = bins.nm[length(bins.nm)] - bins.nm[1])

    # Identify column classes of dataframe for usage later
    col.classes <- lapply(tmp.df, class)

    # Find columns containing the POSIX class
    time.index = c(grep("POSIXt", col.classes))

    tmp.df <- tmp.df %>%
      select(all_of(time.index), everything())

    return(tmp.df)
  })

  # Remove empty data.frames
  data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)
}
