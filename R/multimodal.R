#' Title: multimodal.fitting.R
#'
#' @author Christopher N. Rapp
#'
#' @description
#' An R function that fits a lognormal multimodal particle size distribution or MPSD to
#' a measured dataset using an iterative subtractive Levenberg–Marquardt NLS or LM-NLS
#' algorithm. The curve fitting is based solely on the theoretical lognormal PSD
#' assumed for both ambient and laboratory aerosol.
#'
#' @references Atmospheric Measurement Techniques Preprint: https://doi.org/10.5194/egusphere-2025-4222
#'
#' @param data A \code{data.frame} (or object coercible to
#'   \code{data.table}) containing exactly one \code{POSIXct} timestamp
#'   column and one or more binned concentration columns whose names are
#'   the diameter midpoints of the corresponding bins (e.g. \code{"15.4"}).
#'   Concentration values are assumed to be in units of \code{dN/dlogDp}.
#'   Particle diameter units are arbitrary but must be consistent with
#'   \code{lower.limit} and \code{upper.limit}.
#' @param log.path Character string giving the directory in which the run log
#'   should be written. A timestamped \code{.log} file is created inside this
#'   directory. Defaults to \code{tempdir()} if not supplied.
#' @param frequency Numeric. Time-averaging interval, in minutes, used to group
#'   observations (via \code{lubridate::round_date()}) before fitting each group
#'   separately. If \code{NA}, \code{NULL}, or \code{NaN}, no temporal averaging is
#'   performed and the full dataset is treated as a single time group.
#' @param max.iterations Integer. Maximum number of candidate-peak iterations attempted
#'   per time group before the subtractive fitting loop is forcibly terminated, whether
#'   or not convergence (see \code{FVU.tolerance}) has been reached. Default \code{20}.
#' @param max.modes Integer. Maximum number of lognormal modes that may be fit to a
#'   single time group. Also used to size the plotting color palette. Default \code{5}.
#' @param lower.limit Numeric. Lower bound (in the same units as the bin \code{Dp}
#'   values) of the diameter range over which fitted mode curves are interpolated and
#'   predicted. Default is \code{10} nanometers.
#' @param upper.limit Numeric. Upper bound (in the same units as the bin \code{Dp}
#'   values) of the diameter range over which fitted mode curves are interpolated and
#'   predicted. Default is \code{1000} nanometers.
#' @param verbose Logical. If \code{TRUE}, progress and diagnostic messages are printed
#'   to the console in addition to being written to the log file; if \code{FALSE}, they
#'   are written to the log file only. Default \code{FALSE}.
#' @param plotting Logical. If \code{TRUE}, a three-panel diagnostic \code{ggplot}
#'   (fitted distribution with individual modes, residuals, and predicted-vs-actual
#'   concentration) is generated for each time group via \code{patchwork}. Default
#'   \code{TRUE}.
#' @param labeling Logical. If \code{TRUE} and \code{plotting = TRUE}, each fitted mode
#'   is annotated with its mode number at its peak in the top diagnostic panel. Ignored
#'   if \code{plotting = FALSE}. Default \code{TRUE}.
#' @param smoothing Logical. If \code{TRUE}, the residual \code{dN/dlogDp} vector is
#'   smoothed with \code{smooth.spline()} prior to peak detection with
#'   \code{pracma::findpeaks()}. This generally improves robustness to noisy ambient data. If
#'   smoothing fails or is set to \code{FALSE}, peaks are identified from the raw
#'   residuals. Default \code{TRUE}.
#' @param NMRSE.threshold Numeric. Maximum acceptable normalized RMSE (NRMSE) of the
#'   final combined fit; if the achieved NRMSE exceeds this value, the returned
#'   \code{pass} flag is set to \code{FALSE}. Default \code{0.05}.
#' @param FVU.threshold Numeric. Maximum acceptable Fraction of Variance Unexplained
#'   (FVU) of the final combined fit, expressed as a percentage (0-100); if exceeded,
#'   the returned \code{pass} flag is set to \code{FALSE}. Default \code{10}.
#' @param FVU.tolerance Numeric. Convergence tolerance for FVU, expressed as a
#'   percentage (0-100), used as the stopping criterion for the mode-fitting loop: once
#'   the running FVU drops below this value, no further modes are added to prevent over-fitting.
#'   Default \code{1}.
#' @param low.concentration Numeric. Minimum total particle number concentration
#'   (summed across all bins) required to attempt model fitting for a given time group.
#'   Groups below this threshold are skipped, returned with \code{pass = FALSE}, and a
#'   warning is written to the log. Default \code{100}.
#'
#' @returns A named list, with one element per time group (as defined by
#'   \code{frequency}), where each element is itself a 6-element list containing:
#'   \describe{
#'     \item{pass}{Logical flag indicating whether the fit for this time group met all
#'       quality thresholds (\code{NMRSE.threshold}, \code{FVU.threshold}, and a minimum
#'       concentration of \code{low.concentration}). Useful for filtering results.}
#'     \item{plot}{A \code{patchwork}/\code{ggplot} object with the diagnostic plot for
#'       this time group (\code{NULL} if \code{plotting = FALSE} or the fit failed).}
#'     \item{data}{A \code{data.frame} of predicted vs. actual \code{dN/dlogDp} and
#'       \code{dN} values, residuals, and ratios, by diameter bin.}
#'     \item{predict}{A \code{data.frame} of the interpolated model curve (each mode and
#'       their sum) evaluated between \code{lower.limit} and \code{upper.limit}.}
#'     \item{fits}{A \code{data.frame}/list summarizing each mode's fitted parameters
#'       (\code{N}, \code{GSD}, \code{Dpg}), fit statistics (BIC, RSS, TSS, R2, p-values),
#'       and peak characteristics, along with the underlying \code{nls} model objects.}
#'     \item{evaluation}{A one-row \code{data.frame} of overall goodness-of-fit
#'       statistics for the total MPSD (Pearson correlation, RMSE, NRMSE, dN RMSE,
#'       dN NRMSE, paired t-test p-value, chi-squared p-value) and elapsed fitting time in seconds.}
#'   }
#'
#' @importFrom stats smooth.spline predict coef deviance BIC var cor t.test
#'   chisq.test na.omit relevel
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import patchwork
#'
#' @export

multimodal.fitting <- function(data, log.path, plotting, labeling, frequency, max.iterations, max.modes, smoothing, lower.limit, upper.limit, NMRSE.threshold, FVU.threshold, FVU.tolerance, low.concentration, verbose){

  # -------------------------------------------------------------------------- #
  # SECTION: Default Arguments ----
  # -------------------------------------------------------------------------- #

  if(missing(labeling)) labeling <- T
  if(missing(plotting)) plotting <- T
  if(missing(max.iterations)) max.iterations <- 20
  if(missing(max.modes)) max.modes <- 5
  if(missing(smoothing)) smoothing <- T
  if(missing(lower.limit)) lower.limit <- 10
  if(missing(upper.limit)) upper.limit <- 1000
  if(missing(NMRSE.threshold)) NMRSE.threshold <- 0.05
  if(missing(FVU.threshold)) FVU.threshold <- 10
  if(missing(FVU.tolerance)) FVU.tolerance <- 1
  if(missing(verbose)) verbose <- FALSE
  if(missing(log.path)) log.path <- tempdir()
  if(missing(low.concentration)) low.concentration <- 100

  # -------------------------------------------------------------------------- #
  # SECTION: Helper Function #####
  #' Used in generating the output list

  model.output <- function(pass = FALSE, plot = NULL, data = NULL, predict = NULL, fits = NULL, evaluation = NULL, benchmark.start = NULL){

    # Benchmarking
    benchmark.end <- proc.time()["elapsed"]

    if (!is.null(benchmark.start)) {
      elapsed <- benchmark.end - benchmark.start
    } else {
      elapsed <- NA_real_
    }

    # Ensure evaluation exists
    if (is.null(evaluation)){

      # Dataframe for failed fittings
      evaluation <- data.frame(
        `Pearson Correlation` = NA_real_,
        `RMSE` = NA_real_,
        `NRMSE` = NA_real_,
        `dN RMSE` = NA_real_,
        `dN NRMSE` = NA_real_,
        `Students T Test` = NA_real_,
        `Chi-Squared` = NA_real_,
        `Elapsed Time` = elapsed,
        check.names = FALSE
      )

    } else {
      evaluation$`Elapsed Time` <- elapsed
    }

    # Force exact 6-element structure
    result <- list(
      pass = pass,
      plot = plot,
      data = data,
      predict = predict,
      fits = fits,
      evaluation = evaluation
    )

    return(result)
  }

  # -------------------------------------------------------------------------- #
  # SECTION: Logger #####
  #'

  {
    # Identify column classes of dataframe for usage later
    col.classes <- lapply(data, class)

    # Convert to data table for column indexing
    data <- data.table::as.data.table(data, check.names = FALSE)

    # Find columns containing the POSIX class
    time.index = c(grep("POSIXt", col.classes))[1]

    # Proceed if multiple time values exist
    if (length(time.index) >= 1){

      # Numeric timestamp of first time observation
      timestamp = unlist(data[1, time.index, with = F])

      # Generate string from this timestamp for the logfile nam
      filename <- format(as.POSIXct(timestamp, origin = "1970-01-01"), "%Y%m%d%H%M%S")
    } else {
      filename <- "NA"
    }

    log.path <- file.path(paste0(log.path, "/multimodal", filename, "_",  format(Sys.time(), "%Y%m%d%H%M%S"), ".log"))

    if (verbose){
      print(paste0("Log Path: ", log.path))
    }

    LOG <- logr::log_open(log.path, show_notes = F)
  }

  # -------------------------------------------------------------------------- #
  # SECTION: Initialize Storage Lists #####
  #'

  flag.control <- TRUE
  export.gg <- list()
  export.df <- list()
  export.ft <- list()
  export.lm <- list()
  export.pf <- list()

  # -------------------------------------------------------------------------- #
  # SECTION: Pre Processing #####
  #'

  {
    # Identify column classes of dataframe for usage later
    col.classes <- lapply(data, class)

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: POSIXct Identification #####
    #'

    {
      # Find columns containing the POSIX class
      time.index = c(grep("POSIXt", col.classes))

      if (length(time.index) >= 1){

        # See if the time zone is UTC
        if (any(which(stringr::str_detect(colnames(data), pattern = c("UTC|GMT"))))){

          timezone <- "UTC"

          # Select UTC instance
          time.index <- which(stringr::str_detect(colnames(data), pattern = c("UTC|GMT")))
        } else{

          timezone <- "Unknown"

          # Select first instance
          time.index <- time.index[1]
        }

        tmp.print <- paste0("Current Dataset Time: ",
                            lubridate::as_datetime(as.numeric(dplyr::first(data[, time.index, with = F]))),
                            " ", timezone)

        if (verbose){

          # Print to console
          print(tmp.print)

          # Print to log
          logr::log_print(tmp.print, console = FALSE)
        } else {
          # Print to log
          logr::log_print(tmp.print, console = FALSE)
        }

      } else {

        tmp.print <- "The data set does not contain any detectable POSIXct POSIXt time objects, please format your data accordingly"

        # Print to console
        print(tmp.print)

        # Print to log
        logr::log_print(tmp.print, console = FALSE)
      }
    }

    # ------------------------------------------------------------------------ #
    # SUBSECTION: Sampling Frequency #####
    #'

    {
      # Calculate the time difference between samples in minutes
      tmp.diff <- difftime(data[[time.index]][-1], data[[time.index]][-nrow(data)], units = "mins")
      time.interval <- round(mean(tmp.diff, na.rm = T), 1)

      tmp.print <- paste0("Dataset sampling frequency is ", time.interval, " min")

      if (verbose){
        # Print to console
        print(tmp.print)

        # Print to log
        logr::log_print(tmp.print, console = FALSE)
      } else {
        # Print to log
        logr::log_print(tmp.print, console = FALSE)
      }
    }

    # ------------------------------------------------------------------------ #
    # SUBSECTION: Temporal Averaging #####
    #'

    {
      if (any(c(is.na(frequency), is.null(frequency), is.nan(frequency)))){
        Group = lubridate::as_datetime(as.numeric(dplyr::first(data[, time.index, with = F])))
      } else {
        Group = lubridate::round_date(data[[time.index]], paste(frequency, " mins"))
      }

      # Averaging frequency for use in multimodal analysis
      tmp.df <- data %>%
        dplyr::mutate(`Group` = Group, .after = everything())

      # Split by group
      data.ls <- split(tmp.df, tmp.df$Group)

      # Remove empty lists to prevent errors downstream
      data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)
    }
  }

  # -------------------------------------------------------------------------- #
  # SECTION: Primary Loop #####
  #' Primary loop of the time averaged list
  #' Benchmarking

  export.list <- lapply(data.ls, function(z){

    # Benchmarking start
    benchmark.start <- proc.time()["elapsed"]

    # Used for logging
    date.time.c <- lubridate::as_datetime(as.numeric(dplyr::first(z[, time.index, with = F])))
    date.time.end.c <- lubridate::as_datetime(as.numeric(dplyr::last(z[, time.index, with = F])))

    # Number of observations per group
    n.obs <- nrow(z)

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Bin Identification and dlogDp #####
    #'

    {
      # Retrieve binned diameters and associated columns
      bins.ix = stringr::str_which(colnames(z), "\\d+(\\.\\d+)?$")
      bins.nm = as.numeric(colnames(z)[bins.ix])

      # Subset all non binned data from dataframe
      tmp.df <- z %>%
        dplyr::select(!all_of(bins.ix))

      # Calculate dlogDp
      # Uses forward difference where the gap between starting bin and next one is used
      {
        # Stagger the bins i.e trim end points for each to subtract
        tmp1 <- as.numeric(bins.nm[-c(length(bins.nm))])/1000
        tmp2 <- as.numeric(bins.nm[-1])/1000

        # Find log difference between bin diameters
        dlogDp = log10(tmp2) - log10(tmp1)
      }

      # Create a vector to select columns
      # Due to needing to calculate dlogDp using forward differences the last column is dropped
      # As nucleation/coagulation processes are more sensitive, fitting preferentially needs lower sizes
      select.c <- bins.ix[1:(length(bins.ix) - 1)]
    }

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Calculate dN #####
    #'

    {
      # Subset data and create a numeric matrix
      dNdlogDp.mat <- as.matrix(z[, select.c, with = FALSE])

      # Logical break if mismatch between width of dlogDp and dNdlogDp matrix
      stopifnot(length(dlogDp) == ncol(dNdlogDp.mat))

      # Calculate dN by multiplying dNdlogDp by dlogDp
      dN.mat <- sweep(dNdlogDp.mat, 2, dlogDp, "*")

      # This is to prevent single row data from throwing an error
      if (is.null(nrow(dNdlogDp.mat))){
        dNdlogDp <- dNdlogDp.mat
      } else {
        dNdlogDp <- colMeans(dNdlogDp.mat, na.rm = T)
      }

      # This is to prevent single row data from throwing an error
      if (is.null(nrow(dN.mat))){
        dN <- dN.mat
      } else {
        dN <- colMeans(dN.mat, na.rm = T)
      }

      # Test if sum of all particle density is lower than 100 n/cc by default
      if (sum(dN.mat, na.rm = T) < low.concentration){

        flag.control <- FALSE

        tmp.print <- paste0(date.time.c, ": Total concentration is lower than ", low.concentration, " n/cc, model fitting is not possible for this dataset")

        if (verbose){
          # Print to console
          print(tmp.print)

          # Print to log
          logr::log_print(tmp.print, console = FALSE)
        } else {
          # Print to log
          logr::log_print(tmp.print, console = FALSE)
        }

        # Empty  values for easier binding for total dataset later
        {
          stats.c <- c("Pearson Correlation", "RMSE", "NRMSE", "dN RMSE", "dN NRMSE", "Students T Test", "Chi-Squared")
          stats.nm <- rep(NA, 7)

          export.pf <- as.data.frame(t(stats.nm))
          colnames(export.pf) <- stats.c
          rownames(export.pf) <- NULL

          export.df <- data.frame(
            "Dp" = as.numeric(names(dN)),
            "Predicted dNdlogDp" = NA,
            "Predicted dN" = NA,
            "Actual dNdlogDp" = dNdlogDp,
            "Actual dN" = dN,
            "Residual dNdlogDp" = NA,
            "Residual dN" = NA,
            "Ratio" = NA,
            check.names = F
          )
        }

        return(
          model.output(
            pass = flag.control,
            plot = NULL,
            data = export.df,
            predict = NULL,
            fits = NULL,
            evaluation = export.pf,
            benchmark.start = benchmark.start
          )
        )
      }
    }

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Initialize Model Parameters #####
    #'

    {
      # Create temporary dataset used to initialize the fitting algorithm
      tmp.data = data.frame(Dp = as.numeric(names(dN)),
                            dlogDp = dlogDp,
                            dN = dN,
                            dNdlogDp = dNdlogDp)

      # Remove data containing NAs
      tmp.data <- na.omit(tmp.data)

      # Empty lists
      peak.fitting <- list()
      model.fitting <- list()
      export.lm <- list()
      slopes <- list()

      # Initial parameters to initialize while loop
      FVU = 100
      FVU.i = 100
      i = 1 # iteration counter
      m = 1 # used in peak indexing
    }

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Iterative Peak Identified Levenberg-Marquardt NLS Algorithm #####
    #'

    {
      # Start multi-conditional while loop
      while (FVU > FVU.tolerance & i < max.iterations & length(purrr::compact(model.fitting)) < max.modes){

        {
          # First iteration default
          if (i == 1){
            tmp.data$residual = tmp.data$dNdlogDp
          }

          # ------------------------------------------------------------------ #
          ### Peak Identification #####
          #'

          {
            # Create a smoothed residual vector for more robust peak identification
            smoothingSpline <- tryCatch(
              expr = smooth.spline(tmp.data$Dp, tmp.data$residual, cv = T),
              error = function(e) NULL
            )

            smooth.break <- is.null(smoothingSpline)

            if (smooth.break || isFALSE(smoothing)){

              # Peak identification
              peaks.df <- as.data.frame(pracma::findpeaks(tmp.data$residual,
                                                          minpeakdistance = 5,
                                                          sortstr = T
              ))
            } else if (isFALSE(smooth.break) && isTRUE(smoothing)){

              peaks.df <- as.data.frame(pracma::findpeaks(smoothingSpline$y,
                                                          minpeakdistance = 5,
                                                          sortstr = T
              ))
            }

            # Use indices to match corresponding diameters
            peaks.df$Max <- peaks.df$V1
            peaks.df$Lower <- tmp.data$Dp[peaks.df$V3]
            peaks.df$Mode <- tmp.data$Dp[peaks.df$V2]
            peaks.df$Upper <- tmp.data$Dp[peaks.df$V4]
            peaks.df$Width = peaks.df$V4 - peaks.df$V3
            peaks.df$Range = peaks.df$Upper - peaks.df$Lower

            # Stop code from continuing if the all peaks have been analyzed
            if (m > nrow(peaks.df)) {

              tmp.print <- paste0(date.time.c, ": Terminating loop, no remaining peaks")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              break
            }

            # Parameters used to select data
            peaks.df <- peaks.df[m, ]
            range.c <- peaks.df$V3:peaks.df$V4
          }

          # ------------------------------------------------------------------ #
          # Total Number Concentration (N) Calculation #####
          #'

          {
            # Find both dNdlogDp and dlogDp to convert to dN
            tmp.dNdlogDp = tmp.data$residual[range.c]
            tmp.dlogDp = tmp.data$dlogDp[range.c]

            # Find the original diameters corresponding to the peak range
            tmp.Dp = tmp.data$Dp[range.c]

            # Total concentration
            if (sum(tmp.dNdlogDp*tmp.dlogDp) > 0){

              # Add column for total number concentration (N)
              peaks.df$N = sum(tmp.dNdlogDp*tmp.dlogDp)
            } else {

              # Stop if negative concentration
              tmp.print <- paste0(date.time.c, ": Current Loop: ", i, ", Negative Concentration!")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              # Increase iteration counter (i) and peak index (m)
              m = m + 1
              i = i + 1

              next
            }
          }

          # ------------------------------------------------------------------ #
          # Geometric Standard Deviation (GSD) Calculation #####
          #'

          {
            # Prevent negative values from appearing in GSD calculation
            # This also prevents curves that are predicting negative concentrations
            tmp.dNdlogDp[which(tmp.dNdlogDp < 0)] <- 0
            tmp.dlogDp[which(tmp.dlogDp < 0)] <- 0

            # Calculate GSD
            # Warnings are suppressed for function output
            GSD = suppressWarnings({
              10^(sqrt(sum(tmp.dNdlogDp*((log10(tmp.Dp) - log10(peaks.df$Mode))^2))/(peaks.df$N - 1)))
            })

            # If GSD fails, increase iteration counter and try the next peak
            if (is.nan(GSD)){

              tmp.print <- paste0(date.time.c, ": Current Loop: ", i, ", GSD Error!")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              # Increase iteration counter (i) and peak index (m)
              m = m + 1
              i = i + 1

              next
            } else {

              peaks.df$GSD = GSD
            }
          }


          # ------------------------------------------------------------------ #
          # LM-NLS Model #####
          #' Try fitting an NLS model to the residuals (or original data if first iteration)
          #' Uses starting values from the calculations above

          {
            # Temporary data.frame for LM-NLS model
            # This modifies the available "view" NLS has when trying to curve fit
            tmp <- data.frame(x = tmp.data$Dp[range.c], y = tmp.data$residual[range.c])

            # Run LM-NLS with error handling for graceful termination
            NLS.MODEL <- tryCatch(
              expr = suppressWarnings(minpack.lm::nlsLM(y ~ dNdlogDp.PSD(x, N, GSD, Dpg), data = tmp,
                                                        start = list(N = peaks.df$N, GSD = peaks.df$GSD, Dpg = peaks.df$Mode),
                                                        trace = F)),
              error = function(e) NULL
            )

            # Set to logical
            NLS.break <- is.null(NLS.MODEL)

            # Stop code from continuing if NLS failed
            if (NLS.break) {

              tmp.print <- paste0(date.time.c, ": Terminating loop ", i, ", model fitting failure")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              break
            }

            # Use NLS output to predict data at a resolution of 0.01 units
            tmp.predict = data.frame(x = seq(lower.limit, upper.limit, by = 0.01))
            tmp.predict[[paste0("Mode ", i)]] <- predict(NLS.MODEL, newdata = tmp.predict)
          }

          # ------------------------------------------------------------------ #
          # Quality Control Pass #####
          #

          {
            # Diameter from predicted
            tmp.xL <- round(tmp.predict$x, 2)

            # Diameter from measured
            tmp.xR <- round(tmp.data$Dp, 2)

            # Find matches in predict
            index.match <- which(tmp.xL %in% tmp.xR)

            # Retrieve diameter and corresponding dNdlogDp value from predict
            predict.diameter = tmp.predict[index.match, 1]
            predict.dNdlogDp = tmp.predict[index.match, 2]

            # Throw error if mismatch in dimensions
            # This is due to prediction range not reflecting measured data range
            # Although this could be easily replaced, it is currently enforced to
            # force end-user to at least recognize the size range of their dataset
            if (length(predict.diameter) != nrow(tmp.data)){

              tmp.print <- paste0(date.time.c, ": Error, please modify lower and upper limits to accommadate data set")

              # Print to console
              print(tmp.print)

              # Print to log
              logr::log_print(tmp.print, console = FALSE)

              break
            }

            # Statistical significance of fitted model parameters for the identified
            # peak
            p.N <- summary(NLS.MODEL)$parameters[, 4][1]
            p.GSD <- summary(NLS.MODEL)$parameters[, 4][2]
            p.Dpg <- summary(NLS.MODEL)$parameters[, 4][3]

            # Prevent erroneous over prediction
            peak.height.ratio = max(tmp.predict[, 2])/peaks.df$Max

            # Threshold of 1.5 developed with datasets described in https://doi.org/10.5194/egusphere-2025-4222
            if (peak.height.ratio > 1.5){

              # Mask residuals within identified peak range
              tmp.data$residual[which(tmp.data$Dp < peaks.df$Upper & tmp.data$Dp > peaks.df$Lower)] <- 0

              # Print statement
              tmp.print <- paste0(date.time.c, ": Current Loop: ", i, ", Peak estimation is 150% higher than data, resetting to zero and retrying")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              # Increase iteration counter (i) and peak index (m)
              m = m + 1
              i = i + 1

              next
            } else if (isTRUE(any(p.N >= 0.05, p.GSD >= 0.05, p.Dpg >= 0.05))){

              # Print statement
              tmp.print <- paste0(date.time.c, ": Current Loop: ", i, ", LM-NLS significance greater than 0.05")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                logr::log_print(tmp.print, console = FALSE)
              }

              # Increase iteration counter (i) and peak index (m)
              m = m + 1
              i = i + 1

              next
            } else { # Successful fit

              # Calculate residual within loop
              # Note this includes a masked residual
              tmp.diff <- tmp.data$residual - predict.dNdlogDp

              # Replace
              tmp.data$residual <- tmp.diff

              # Calculate FVU
              FVU.i = round((var(tmp.data$residual)/var(tmp.data$dNdlogDp))*100, 2)

              # Add to peaks.df for exporting
              peaks.df$FVU <- FVU.i

              # Add model values to output list
              model.fitting[[i]] <- tmp.predict

              # Add peak values to output list
              peak.fitting[[i]] <- peaks.df

              # Add model values to output list
              export.lm[[i]] <- NLS.MODEL
            }

            # If newest fit's FVU is an improvement, replace existing
            if (FVU.i < FVU){

              FVU <- FVU.i
            }

            tmp.print <- paste0(date.time.c, ": Current Loop: ", i, ", Remaining Variance: ", FVU, "%",
                                ", # of Modes: ", nrow(data.table::rbindlist(purrr::compact(peak.fitting))))

            # Increase iteration counter (i)
            i <- i + 1

            if (verbose){

              # Print to console
              print(tmp.print)

              # Print to log
              logr::log_print(tmp.print, console = F)
            } else {
              # Print to log
              logr::log_print(tmp.print, console = FALSE)
            }
          }
        }
      }
    }

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Post Processing #####
    #'

    {
      # ---------------------------------------------------------------------- #
      ### Fits (export.lm) #####
      # Fit quality and peak information

      {
        # Bind list output
        tmp.ls <- purrr::compact(export.lm)

        if (length(tmp.ls) == 0){

          # If model fails to initialize due to early overestimation (noise in dataset typically)
          # Export actual data with filled NA's for anything the model would have provided
          export.df <- data.frame(
            "Dp" = tmp.data$Dp,
            "Predicted dNdlogDp" = NA,
            "Predicted dN" = NA,
            "Actual dNdlogDp" = tmp.data$dNdlogDp,
            "Actual dN" = tmp.data$dN,
            "Residual dNdlogDp" = NA,
            "Residual dN" = NA,
            "Ratio" = NA,
            check.names = F
          )

          return(
            model.output(
              pass = FALSE,
              plot = NULL,
              data = export.df,
              predict = NULL,
              fits = NULL,
              evaluation = export.pf,
              benchmark.start = benchmark.start
            )
          )
        }

        # Create blank dataframe
        tmp.df <- data.frame(matrix(nrow = length(tmp.ls), ncol = 11))

        # Fill dataframe
        for (s in 1:length(tmp.ls)){

          tmp.df[s, 1] <- paste0("Mode ", s)
          tmp.df[s, 2] <- c(coef(tmp.ls[[s]]))[1] # N
          tmp.df[s, 3] <- c(coef(tmp.ls[[s]]))[2] # GSD
          tmp.df[s, 4] <- c(coef(tmp.ls[[s]]))[3] # Dpg
          tmp.df[s, 5] <- BIC(tmp.ls[[s]]) # Bayesian information criterion
          tmp.df[s, 6] <- deviance(tmp.ls[[s]]) # Residual sum of squares
          tmp.df[s, 7] <- sum((tmp.ls[[s]]$m$getEnv()$y - mean(tmp.ls[[s]]$m$getEnv()$y))^2) # Total sum of squares
          tmp.df[s, 8] <- 1 - tmp.df[s, 6]/tmp.df[s, 7] # R2
          tmp.df[s, 9] <- signif(summary(tmp.ls[[s]])$parameters[, 4][1], 1) # N T pval
          tmp.df[s, 10] <- signif(summary(tmp.ls[[s]])$parameters[, 4][2], 2) # GSD T pval
          tmp.df[s, 11] <- signif(summary(tmp.ls[[s]])$parameters[, 4][3], 3) # Dpg T pval
        }

        # Setnames of data.frame
        data.table::setnames(tmp.df, new = c("Mode Label", "N", "GSD", "Dpg", "BIC", "RSS", "TSS", "R2", "N T pval", "GSD T pval", "Dpg T pval"))

        # Combine peak information into 1 dataframe
        peaks.df <- data.table::rbindlist(purrr::compact(peak.fitting))

        # Subset
        peaks.df <- peaks.df %>%
          dplyr::select(c(Max:Width), `FVU`)

        # Export summary of both fitting values and peak parameters
        export.lm <- cbind(tmp.df, peaks.df) %>%
          dplyr::relocate(Max:Width, .after = Dpg)
      }

      # ---------------------------------------------------------------------- #
      # Predict (export.ft) #####
      # Data.frame of predicted values (dNdlogDp) per mode

      {
        model.fitting <- purrr::compact(model.fitting)
        output <- do.call(cbind, model.fitting)
        output <- output[append(1, stringr::str_which(colnames(output), 'Mode'))]
        output <- round(output, 2)

        data.table::setnames(output, old = colnames(output)[stringr::str_which(colnames(output), pattern = "Mode \\d{1,2}")],
                 new = paste0("Mode ", seq_along((stringr::str_which(colnames(output), pattern = "Mode \\d{1,2}")))))

        data.table::setnames(output, old = "x", new = "Dp")

        if (ncol(output) > 2){
          output$dNdlogDp <- rowSums(output[, 2:ncol(output)])
        } else {
          output$dNdlogDp <- output[, 2]
        }

        # Predicted values for each mode
        export.ft <- output
      }

      # ---------------------------------------------------------------------- #
      # Data (export.df) #####
      # Calculate various statistical measures providing "goodness-of-fit"

      {
        # Matching predictions to measurements
        {
          # Find the matched diameters so measurements to predicted binned diameters is correct
          match.ix <- which(output$Dp %in% tmp.data$Dp)

          # Diameters
          Dp <- output$Dp[match.ix]
          dlogDp <- tmp.data$dlogDp

          # Predicted Lognormal Concentration
          predicted.dNdlogDp <- output$dNdlogDp[match.ix]

          # Predicted Concentration
          predicted.dN <- predicted.dNdlogDp*dlogDp

          # Actual lognormal concentration
          actual.dNdlogDp <- tmp.data$dNdlogDp

          # Actual concentration
          actual.dN <- tmp.data$dN
        }

        export.df <- data.frame(
          "Dp" = Dp,
          "Predicted dNdlogDp" = predicted.dNdlogDp,
          "Predicted dN" = predicted.dN,
          "Actual dNdlogDp" = actual.dNdlogDp,
          "Actual dN" = actual.dN,
          "Residual dNdlogDp" = actual.dNdlogDp - predicted.dNdlogDp,
          "Residual dN" = actual.dN - predicted.dN,
          "Ratio" = predicted.dNdlogDp/actual.dNdlogDp,
          check.names = F
        )
      }

      # ---------------------------------------------------------------------- #
      # Model Performance (export.pf) #####
      # Calculate various statistical measures evaluating "goodness-of-fit"

      {
        # Pearson Correlation
        stats.R2 = round(cor(predicted.dN, actual.dN)^2, 4)

        # FVU
        # As the residuals are sometimes set to 0 during the loop, the actual FVU needs to be calculated
        FVU = round((1 - stats.R2), 4)

        # RMSE
        stats.RMSE = round(Metrics::rmse(actual.dNdlogDp, predicted.dNdlogDp), 2)
        dN.RMSE = round(Metrics::rmse(actual.dN, predicted.dN), 2)

        # Max min normalized RMSE
        NRMSE = stats.RMSE/(max(actual.dNdlogDp) - min(actual.dNdlogDp))
        dN.NRMSE = dN.RMSE/(max(actual.dN) - min(actual.dN))

        # Significance testing
        # NOTE: This is not a recommended evaluation method as it has demonstrated to
        # pass or fail even the best fittings
        # USE AT YOUR OWN RISK
        stats.STTEST <- t.test(actual.dNdlogDp, predicted.dNdlogDp, paired = T)
        stats.STTEST <- round(stats.STTEST$p.value, 4)

        stats.CHI <- suppressWarnings(chisq.test(actual.dNdlogDp, predicted.dNdlogDp))
        stats.CHI <- round(stats.CHI$p.value, 4)

        stats.c <- c("Pearson Correlation", "RMSE", "NRMSE", "dN RMSE", "dN NRMSE", "Students T Test", "Chi-Squared")
        stats.nm <- c(stats.R2, stats.RMSE, NRMSE, dN.RMSE, dN.NRMSE, stats.STTEST, stats.CHI)

        # Use data from above to
        export.pf <- as.data.frame(t(stats.nm))
        colnames(export.pf) <- stats.c
        rownames(export.pf) <- NULL

        if (verbose){
          print(paste0("Concentration RMSE: ", stats.RMSE, " n/cc"))
        }

        # FVU check
        if (FVU > FVU.threshold/100){
          flag.control <- FALSE
        }

        # NRMSE check
        if (NRMSE > NMRSE.threshold){
          flag.control <- FALSE
        }

        rm(stats.nm, stats.c)
      }
    }

    # ------------------------------------------------------------------------ #
    ## SUBSECTION: Plotting #####
    #'

    if (plotting){

      # Pivot data into long format
      # For each Dp, there is a concentration for each mode
      plot.df <- output %>%
        tidyr::pivot_longer(
          cols = !c("Dp"),
          names_to = "Mode",
          values_to = "Concentration"
        )

      # Replace label with total for plotting
      plot.df[plot.df == "dNdlogDp"] <- "Total"

      # Set levels for modes so plotting in correct order
      plot.df <- plot.df %>%
        dplyr::mutate(`Mode` = relevel(factor(`Mode`, ordered = F), ref = "Total"))

      # X AXIS
      {
        # Use log distributed values rather than numerical to improve performance
        # Needs several orders of magnitude less data to produce smooth curves
        x = 10^seq(log10(lower.limit), log10(upper.limit), length.out = 1000)

        x.limits = 10^((log10(lower.limit)):(log10(upper.limit)))

        x.breaks <- NULL
        x.labels <- NULL
        for (i in round(log10(x.limits))){
          x.breaks <- append(x.breaks, i)
          x.labels <- append(x.labels, as.character(10^(i)))
        }
      }

      # Y AXIS
      {
        #
        y.limits = c(-1*max(c(export.df$`Predicted dNdlogDp`, export.df$`Actual dNdlogDp`))/10,
                     max(c(export.df$`Predicted dNdlogDp`, export.df$`Actual dNdlogDp`))/10)

        y.breaks = sort(unique(c(pretty(export.df$`Predicted dNdlogDp`), pretty(export.df$`Actual dNdlogDp`))))

        # When the breaks are uneven, force a reset
        if (any(diff(diff(y.breaks)) < 0)){
          y.breaks <- pretty(y.breaks, n = 10)
        } else {
          y.breaks
        }
      }

      # Colorblind friendly palette
      {
        cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        # Gradient to allow more colors up to max.modes
        cbPalette <- grDevices::colorRampPalette(cbPalette)(max.modes)

        # Total line is black by default and the first mode plotted
        cm.palette <- c("black", cbPalette)
      }

      # Trim excess data to prevent ggplot warnings
      plot.df <- plot.df %>%
        dplyr::filter(`Dp` <= x.limits[length(x.limits)] & `Dp` >= x.limits[1])

      # Use the RMSE as the shaded error region around the total fit line
      # Trim excess data to prevent ggplot warnings
      error <- output %>%
        dplyr::mutate(Min = dplyr::if_else(dNdlogDp - stats.RMSE < 0, 0, dNdlogDp - stats.RMSE)) %>%
        dplyr::mutate(Max = dNdlogDp + stats.RMSE) %>%
        dplyr::filter(`Dp` <= x.limits[length(x.limits)] & `Dp` >= x.limits[1])

      # Plot title
      if (n.obs == 1 & any(c(is.na(frequency), is.null(frequency), is.nan(frequency)))){
        plot.title = paste0(date.time.c)
      } else {
        plot.title = paste0(date.time.c, " - ", date.time.end.c)
      }

      # Top concentration values per mode to pin the annotation
      text.labels <- plot.df %>%
        dplyr::group_by(Mode) %>%
        dplyr::arrange(desc(Concentration)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

      # Don't annotate the total line and format to (#)
      text.labels <- text.labels %>%
        dplyr::filter(Mode != "Total") %>%
        dplyr::mutate(label = paste0("(", stringr::str_replace(Mode, "Mode ", ""), ")"))

      # Empty string if labeling is turned off
      if (labeling == FALSE){
        text.labels <- text.labels %>%
          dplyr::mutate(label = "")
      }

      # ---------------------------------------------------------------------- #
      ### 3-Panel Plot #####
      #'

      {
        # Shared base theme — panel/legend/grid styling common to all three panels
        base.theme <- theme(
          plot.subtitle         = element_text(hjust = 0.5, size = 12),
          panel.background      = element_rect(fill = "white"),
          panel.grid.major      = element_line(colour = "grey80", linewidth = 0.1),
          panel.grid.minor      = element_line(colour = "grey80", linewidth = 0.01),
          panel.border          = element_rect(colour = "black", fill = NA),
          axis.title.y          = element_text(angle = 90, vjust = 2),
          legend.background     = element_blank(),
          legend.box.background = element_blank(),
          legend.key            = element_blank(),
          legend.title          = element_blank(),
          legend.position       = "right"
        )

        top.gg <- ggplot(tmp.data, aes(x = Dp, y = dNdlogDp)) +
          geom_point(shape = 1) +
          geom_line(data = plot.df, aes(x = Dp, y = Concentration, color = Mode), linewidth = 0.5) +
          geom_ribbon(data = error, aes(x = Dp, ymin = Min, ymax = Max), inherit.aes = F, fill = "gray75", alpha = 0.5) +
          scale_x_log10(breaks = 10^x.breaks, labels = x.labels, limits = c(min(x.limits), max(x.limits))) +
          scale_y_continuous(breaks = y.breaks, limits = range(y.breaks)) +
          ggrepel::geom_text_repel(data = text.labels, aes(x = Dp, y = Concentration, label = label),
                                   nudge_x = 0.1, nudge_y = y.limits[2] * 0.1, box.padding = 1) +
          labs(title = plot.title,
               subtitle = paste0("Pass: ", flag.control, ", Concentration RMSE: ", stats.RMSE, " n/cc",
                                 ", NRMSE: ", round(NRMSE, 2), ", (n = ", n.obs, ")")) +
          ylab(expression("dN/dlog"[10] * "D"[p] * "  [" * "cm"^-3 * "]")) +
          xlab(expression("D"[p] * "  [nm]")) +
          scale_color_manual(values = cm.palette) +
          base.theme +
          theme(
            plot.title  = element_text(hjust = 0.5, size = 16),
            plot.margin = margin(1, 1, 0.1, 1, unit = "cm")
          ) +
          guides(x = guide_axis_logticks(long = 3)) +
          coord_cartesian(clip = "off")

        mid.gg <- ggplot(export.df, aes(x = Dp, y = `Residual dNdlogDp`)) +
          geom_point(shape = 1) +
          geom_segment(aes(x = Dp, xend = Dp, y = 0, yend = `Residual dNdlogDp`)) +
          scale_x_log10(breaks = 10^x.breaks, labels = x.labels, limits = c(min(x.limits), max(x.limits))) +
          scale_y_continuous(limits = range(export.df$`Residual dNdlogDp`)) +
          labs(title = paste0("Fraction Variance Unexplained: ", FVU)) +
          ylab(expression("Residual (dN/dlog"[10] * "D"[p] * ")")) +
          xlab(expression("D"[p])) +
          scale_color_manual(values = cm.palette) +
          base.theme +
          theme(
            plot.title  = element_text(hjust = 0.5, size = 12),
            plot.margin = margin(1, 1, 0.1, 1, unit = "cm")
          ) +
          guides(x = guide_axis_logticks(long = 3))

        bot.gg <- ggplot(export.df, aes(x = `Predicted dN`, y = `Actual dN`)) +
          geom_point(shape = 1) +
          geom_abline(slope = 1) +
          scale_y_continuous(n.breaks = 5) +
          xlab(expression("Predicted Concentration " ~ "n cm"^-3)) +
          ylab(expression("Actual Concentration " ~ "n cm"^-3)) +
          scale_color_manual(values = cm.palette) +
          labs(title = paste0("Pearson Correlation: ", stats.R2)) +
          base.theme +
          theme(
            plot.title    = element_text(hjust = 0.5, size = 12),
            plot.caption  = element_text(hjust = 0, vjust = -5),
            axis.title.x  = element_text(angle = 0, vjust = -1),
            plot.margin   = margin(1, 1, 1, 1, unit = "cm")
          )
      }

      # Merge plots using patchwork
      export.gg <- top.gg / mid.gg / bot.gg
    }

    # ------------------------------------------------------------------------ #
    # Exit Primary Loop  #####
    return(
      model.output(
        pass = flag.control,
        plot = export.gg,
        data = export.df,
        predict = export.ft,
        fits = export.lm,
        evaluation = export.pf,
        benchmark.start = benchmark.start
      )
    )
  })

  # ------------------------------------------------------------------------ #
  # END #####

  # Close log
  logr::log_close()

  return(export.list)
}

#' Lognormal particle size distribution
#'
#' Evaluates the theoretical single-mode lognormal particle size distribution
#' function, \eqn{dN/d\log_{10}D_p}, at one or more particle diameters. The
#' distribution is formally defined in natural-log space and internally
#' converted to base-10 logarithms via the \code{2.303} (\eqn{= \ln(10)})
#' scaling factor. Used internally by \code{\link{multimodal.fitting}} as the
#' model formula passed to \code{minpack.lm::nlsLM()}.
#'
#' @param dx Numeric vector of particle diameters (\eqn{D_p}) at which to
#'   evaluate the distribution. Must be strictly positive and in the same
#'   units as \code{Dpg}.
#' @param N Numeric. Total particle number concentration of the mode
#'   (integral of \eqn{dN/d\log_{10}D_p} over all diameters).
#' @param GSD Numeric. Geometric standard deviation of the mode (unitless,
#'   \eqn{> 1}).
#' @param Dpg Numeric. Geometric mean (modal) diameter, in the same units as
#'   \code{dx}.
#'
#' @returns A numeric vector, the same length as \code{dx}, of
#'   \eqn{dN/d\log_{10}D_p} values for the specified lognormal mode.
#'
#' @export
dNdlogDp.PSD <- function(dx, N, GSD, Dpg){

  A = N/(((2*pi)^(1/2))*log(GSD))
  B = -1*(log(dx)-log(Dpg))^2
  C = 2*(log(GSD)^2)

  result = 2.303*A*exp(B/C)

  return(result)
}
