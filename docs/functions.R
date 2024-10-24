## FilterSpec functions -------

#' Filter time series
#'
#' Apply a given filter to a time series using different endpoint constraints.
#'
#' Note that when passing objects of class \code{ts}, the time step provided is
#' not used; thus, for time series with a time step different from 1, the filter
#' has to be adapted accordingly.
#'
#' Leading and trailing NA values are automatically stripped from the input
#' vector so that they do not spread into the filtered data when applying the
#' endpoint constraints, but added in again after filtering so that the output
#' vector has the same length as the input. This does not apply to any internal
#' NA values, which instead are handled by \code{na.rm}.
#'
#' The function applies endpoint constrains following Mann et al., GRL, 2004;
#' available methods are:
#' \itemize{
#'   \item method = 0: no constraint (loss at both ends);
#'   \item method = 1: minimum norm constraint;
#'   \item method = 2: minimum slope constraint;
#'   \item method = 3: minimum roughness constraint;
#'   \item method = 4: circular filtering.
#' }
#'
#' @param data numeric vector with the input timeseries (standard or ts object).
#' @param filter numeric vector of filter weights.
#' @param method single integer for choosing an endpoint constraint method;
#'   available choices are integers 0-4, see details.
#' @param na.rm logical; control the handling of internal NA values in
#'   \code{data}. If set to \code{TRUE}, any internal NA values are removed by
#'   linear interpolation from the neighbouring values; defaults to
#'   \code{FALSE}.
#' @return a ts object with the filtered timeseries.
#' @author Thomas Laepple
#' @source The endpoint constraint methods are based on the study:\cr
#'   Michael E. Mann, On smoothing potentially non‐stationary climate time
#'   series, Geophys. Res. Lett., 31, L07214, doi:10.1029/2004GL019569, 2004.
#' @examples
#' # Simple running mean filter across three bins
#'
#' x <- 1 : 10
#' filter <- rep(1 / 3, 3)
#'
#' # no endpoint constraints lead to loss at both ends
#' ApplyFilter(x, filter, method = 0)
#'
#' # circular filtering avoids end losses, so as the other methods
#' ApplyFilter(x, filter, method = 4)
#'
#' # leading and trailing NA's are ignored but added in again afterwards
#' x <- c(NA, 1 : 10, NA)
#' ApplyFilter(x, filter, method = 4)
#'
#' # ... but not internal NA's
#' x <- c(1 : 5, NA, 7 : 10)
#' ApplyFilter(x, filter, method = 4)
#'
#' # if not explicitly removed by linear interpolation
#' ApplyFilter(x, filter, method = 4, na.rm = TRUE)
#'
#'
#' # Visual comparison of methods ----------
#' set.seed(20220302)
#'
#' x <- PaleoSpec::SimPowerlaw(1, 1e02)
#' x <- x + 0.1 * (1:length(x))
#'
#' filt <- rep(1/30, 30)
#'
#' plot(x, type = "l")
#'
#' x0 <- ApplyFilter(x, filt, method = 0)
#' lines(x0, col = "Green")
#'
#' x1 <- ApplyFilter(x, filt, method = 1)
#' lines(x1, col = "blue")
#'
#' x2 <- ApplyFilter(x, filt, method = 2)
#' lines(x2, col = "red")
#'
#' x3 <- ApplyFilter(x, filt, method = 3)
#' lines(x3, col = "orange")
#'
#' x4 <- ApplyFilter(x, filt, method = 4)
#' lines(x4, col = "Purple", lty = 2)
#'
#' lines(x0, col = "Green")
#'
#' legend(x = "topleft",
#'        legend = c("0: no constraint (lose ends)",
#'                   "1: min norm (pad with mean)",
#'                   "2: min slope (reflect ends)",
#'                   "3: min roughness (reflect and invert ends)",
#'                   "4: circular"),
#'        col = c("Green", "Blue", "Red", "Orange", "Purple"),
#'        lwd = 2, lty = c(1,1,1,1,2))
#'
#'
#' # Repeat with linear trend, no noise
#'
#' x <- 1:100
#' filt <- rep(1/30, 30)
#'
#' plot(x, type = "l")
#'
#' x0 <- ApplyFilter(x, filt, method = 0)
#' lines(x0, col = "Green")
#'
#' x1 <- ApplyFilter(x, filt, method = 1)
#' lines(x1, col = "blue")
#'
#' x2 <- ApplyFilter(x, filt, method = 2)
#' lines(x2, col = "red")
#'
#' x3 <- ApplyFilter(x, filt, method = 3)
#' lines(x3, col = "orange")
#'
#' x4 <- ApplyFilter(x, filt, method = 4)
#' lines(x4, col = "Purple", lty = 2)
#'
#' lines(x0, col = "Green")
#'
#' legend(x = "topleft",
#'        legend = c("0: no constraint (lose ends)",
#'                   "1: min norm (pad with mean)",
#'                   "2: min slope (reflect ends)",
#'                   "3: min roughness (reflect and invert ends)",
#'                   "4: circular"),
#'        col = c("Green", "Blue", "Red", "Orange", "Purple"),
#'        lwd = 2, lty = c(1,1,1,1,2))
#' @export
ApplyFilter <- function(data, filter, method = 0, na.rm = FALSE) {
  
  if (!method %in% (0 : 4))
    stop("Unknown method; only 0 : 4 available.")
  
  result <- rep(NA, length(data))
  
  # remove leading and trailing NA's
  x <- c(zoo::na.trim(data))
  n <- length(x)
  
  # linearly interpolate internal NA's if requested
  if (na.rm) {x <- stats::approx(1 : n, x, 1 : n)$y}
  
  circular = FALSE
  
  if (method == 0 | method == 4) {
    
    if (method == 4) {circular = TRUE}
    
    xf <- stats::filter(x, filter, circular = circular)
    
  } else {
    
    N <- floor(length(filter) / 2)
    
    if (method == 1) {
      
      before <- rep(mean(x), N)
      after  <- rep(mean(x), N)
      
    } else if (method == 2 | method == 3) {
      
      before <- x[N : 1]
      after  <- x[n : (n - N + 1)]
      
      if (method == 3) {
        
        before <- x[1] - (before - mean(before))
        after  <- x[n] - (after - mean(after))
        
      }
    }
    
    xf <- stats::filter(c(before, x, after), filter, circular = circular)
    xf <- xf[(N + 1) : (N + n)]
    
  }
  
  i <- seq(match(x[1], data), by = 1, length.out = n)
  result[i] <- xf
  
  return(ts(result, frequency = frequency(data)))
  
}


#' Filter a Power Spectrum Object
#'
#' @param spec A spec object
#' @param keep_low_f Keep filtered (smoothed) low frequencies or replace with unfiltered
#' @inheritParams stats::spec.pgram
#' @inheritParams ApplyFilter
#' @return A spec object (list)
#' @seealso [FilterSpecLog] for filtering with filter widths equal in log-space
#' @export
#'
#' @examples
#' ## Comparison of the four methods - for power spectra, methods 0, 2 or 3 make the most sense
#' library(PaleoSpec)
#'
#' a <- 100
#' b <- 1
#' N <- 1e03
#' set.seed(20230625)
#' ts1 <- SimPLS(N, beta = b, alpha = a)
#' sp1 <- SpecMTM(ts(ts1), bin.width = 1)
#' LPlot(sp1)
#' abline(log10(a), -b, col = "green")

#' fl <- seq(3, 9, by = 2)
#' sp1_f3_0 <- FilterSpec(sp1, spans = fl, method = 0)
#' sp1_f3_1 <- FilterSpec(sp1, spans = fl, method = 1)
#' sp1_f3_2 <- FilterSpec(sp1, spans = fl, method = 2)
#' sp1_f3_3 <- FilterSpec(sp1, spans = fl, method = 3)
#' sp1_f3_4 <- FilterSpec(sp1, spans = fl, method = 4)
#'
#' LPlot(sp1)
#' LLines(sp1_f3_0, col = "blue")
#' LLines(sp1_f3_1, col = "green", lty = 2)
#' LLines(sp1_f3_2, col = "red", lty = 3)
#' LLines(sp1_f3_3, col = "orange", lty = 4)
#' LLines(sp1_f3_4, col = "gold", lty = 5)
#'
#' ## Comparison of keeping the filtered values in the reflected end portions or not
#' sp1_f3_0T <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = TRUE)
#' sp1_f3_0F <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = FALSE)

#' LPlot(sp1_f3_0F)
#' LLines(sp1_f3_0T, col = "red")

#' sp1_f3_2T <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = TRUE)
#' sp1_f3_2F <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = FALSE)

#' LPlot(sp1_f3_2F)
#' LLines(sp1_f3_2T, col = "red")
FilterSpec <- function(spec, spans, method = 3, keep_low_f = TRUE) {
  if (length(spec$dof) == 1) {
    spec$dof <- rep(spec$dof, length(spec$freq))
  }
  
  dof0 <- spec$dof
  
  kernel <- stats::kernel("modified.daniell", spans %/% 2)
  filter <- kernel[-kernel$m:kernel$m]
  
  spec_filt <- ApplyFilter(spec$spec, filter = filter, method = method)
  
  if (keep_low_f == FALSE) {
    # replace filtered spec with original in area where freqs have been reflected
    i <- 1:ceiling(length(filter) / 2)
    spec_filt[i] <- spec$spec[i]
    
    iend <- length(spec$freq) - (i-1)
    
    spec_filt[iend] <- spec$spec[iend]
    
  }
  
  spec$spec <- as.numeric(spec_filt)
  
  # degrees of freedom of the kernel
  df.kern <- stats::df.kernel(kernel)
  
  spec$dof <- df.kern * spec$dof / 2
  
  if (keep_low_f == FALSE) {
    
    i <- 1:ceiling(length(filter) / 2)
    spec$dof[i] <- dof0[i]
    
    iend <- length(spec$freq) - (i-1)
    spec$dof[iend] <- dof0[iend]
    
  }
  
  # Adjust DOF in reflected filter region
  if (keep_low_f == TRUE){
    
    fl <- length(filter)
    i <- 1:ceiling(fl / 2)
    iend <- length(spec$freq) - (i-1)
    
    
    if (method %in% c(2,3)){
      scl <- 2 * (fl - (i - 1)) / fl
      spec$dof[i] <- spec$dof[i] / scl
      spec$dof[iend] <- spec$dof[iend] / scl
      
    }
    
    if (method == 0){
      # remove NA portion
      spec$freq <- spec$freq[is.na(spec$spec) == FALSE]
      spec$dof <- spec$dof[is.na(spec$spec) == FALSE]
      spec$shape <- spec$shape[is.na(spec$spec) == FALSE]
      spec$spec <- spec$spec[is.na(spec$spec) == FALSE]
    }
    
  }
  
  spec$shape <- spec$dof / 2
  
  spec <- AddConfInterval(spec)
  
  
  return(spec)
}



#' Smooth a Spectrum with Evenly Spaced Bins in Logspace
#'
#' @param spec A spec object
#' @inheritParams LogSmooth
#' @inheritParams ApplyFilter
#'
#' @return A spec object (list)
#' @seealso [LogSmooth()] for an alternative implementation of log spaced filtering
#' @export
#' @examples
#' library(PaleoSpec)
#'
#' # simulate a timeseries with powerlaw power spectrum
#' a <- 100
#' b <- 1
#' N <- 1e03
#'
#' set.seed(20230625)
#' ts1 <- SimPLS(N, beta = b, alpha = a)
#' sp1 <- SpecMTM(ts(ts1), bin.width = 1)
#' LPlot(sp1)
#' abline(log10(a), -b, col = "green")
#' #
#' sp1_f3_0 <- FilterSpecLog(sp1, method = 0)
#' sp1_f3_2 <- FilterSpecLog(sp1, method = 2)
#'
#' LPlot(sp1)
#' LLines(sp1_f3_0, col = "blue")
#' LLines(sp1_f3_2, col = "green", lty = 3)
#'
#' sp1_df0.05 <- FilterSpecLog(sp1)
#' sp1_df0.1 <- FilterSpecLog(sp1, df.log = 0.1)
#'
#' LPlot(sp1)
#' LLines(sp1_df0.05, col = "blue")
#' LLines(sp1_df0.1, col = "red")
#'
#' ## A combination of FilterSpec and FilterSpecLog
#'
#' sp1_FSL <- FilterSpecLog(sp1)
#' sp1_FSL_FS <- FilterSpec(FilterSpecLog(sp1), spans = c(3, 5))
#' sp1_FS_FSL <- FilterSpecLog(FilterSpec(sp1, spans = c(3, 5)))
#' LPlot(sp1)
#' LLines(sp1_FSL, col = "blue")
#' LLines(sp1_FSL_FS, col = "red")
#' LLines(sp1_FS_FSL, col = "green")
FilterSpecLog <- function(spec,
                          df.log = 0.05,
                          spans = NULL,
                          method = 3, f.res = 10){
  
  GetFW <- function(spec, df.log) {
    ((exp(df.log) - 1) * max(spec$freq)) / min(spec$freq)
  }
  
  if (length(spec$dof) == 1){
    spec$dof <- rep(spec$dof, length(spec$freq))
  }
  
  if (is.null(spans)){
    spans <- GetFW(spec, df.log = df.log)
  }
  
  # interpolate spectrum onto equal in log space freq axis
  delta_f <- min(spec$freq)
  logfreq <- log(spec$freq)
  
  freq_logspace <- (seq(min(logfreq), max(logfreq)+delta_f, length.out = f.res*length(spec$freq)))
  spec_loginterp <- stats::approx(logfreq, spec$spec, xout = freq_logspace, rule = 2)$y
  
  spans_adj <- spans * f.res
  
  # DOF of filter
  kernel <- stats::kernel("daniell", spans_adj %/% 2)
  filter <- kernel[-kernel$m:kernel$m]
  df.kern <- stats::df.kernel(kernel)
  
  # DOF of a boxcar filter the same width
  kernal.flat <- stats::kernel("daniell", length(filter) %/% 2)
  df.kern.flat <- stats::df.kernel(kernal.flat)
  
  # modify for non boxcar filters
  df.mod <- df.kern / df.kern.flat
  
  # smooth/filter in log space
  spec_filt <- ApplyFilter(spec_loginterp, filter = filter, method = method)
  
  # re-interpolate back to original freq axis
  spec3 <- stats::approx(freq_logspace, spec_filt, xout = logfreq)$y
  
  # overwrite spec with filtered spec
  spec$spec <- spec3
  
  # keep old DOF
  dof0 <- spec$dof
  
  # Gets the difference in delta_f for the log and standard freq axis
  NpF <- function(freq, fw, df){
    
    posdiff <- (exp(log(freq) + df) - freq)
    negdiff <- (freq - exp(log(freq) - df))
    
    fdiff <- rowMeans(cbind(negdiff, posdiff))
    
    2 * fw * (fdiff/df) * 1/(2*max(freq))
  }
  df.logkern <- NpF(spec$freq, length(filter), df = diff(freq_logspace[1:2]))
  
  spec$dof <- spec$dof + df.mod * df.logkern * spec$dof/2
  spec$shape <- spec$dof/2
  spec$spans <- paste(spans, collapse = ",")
  
  if (method == 0){
    # remove NA portion
    spec$freq <- spec$freq[is.na(spec$spec) == FALSE]
    spec$dof <- spec$dof[is.na(spec$spec) == FALSE]
    spec$shape <- spec$shape[is.na(spec$spec) == FALSE]
    spec$spec <- spec$spec[is.na(spec$spec) == FALSE]
  }
  
  spec <- AddConfInterval(spec)
  
  
  
  return(spec)
}


## SNRStack ------

#' Title
#'
#' @param x 
#' @param bin_width 
#' @param prefilter 
#' @param logsmooth 
#' @param equalise_var Equalise the variance of all records in a cluster. 
#' Preserves the mean variance of the cluster so that is can be compared with 
#' e.g. instrumental data, but avoids bias to SNR estimate if each recorder had
#'  a different sensitivity.
#' @inheritParams SpecMTM
#' @return
#' @export
#'
#' @examples
SNRStack <- function(x, bin_width,
                     prefilter = FALSE,
                     logsmooth = FALSE,
                     k = 1, nw = 0,
                     #equalise_mean = FALSE,
                     equalise_var = FALSE) {
  
  
  # equalise the variance of all records in a cluster
  if (equalise_var == TRUE){
    
    sd_x <- apply(x, MARGIN = 2, sd, na.rm = TRUE)
    mean_sd_x = mean(sd_x)
    
    x <- apply(x, MARGIN = 2, function(i) {
      mean_sd_x * i / sd(i, na.rm = TRUE)
    })
    
  }
  
  
  ## not sure if n should be the simple number of records
  # n <- ncol(x)
  
  # or the mean number of records per timepoint, which
  # will be less if there are a lot of gaps
  
  n <- apply(x, 1, function(i) sum(is.na(i) == FALSE))
  n <- mean(n, na.rm = TRUE)
  
  # mean record spectrum
  S_record <- SpecACF(x, bin.width = bin_width, k = k, nw = nw)
  # spectrum of the stack
  S_stack <- SpecACF(rowMeans(x, na.rm = TRUE),
                     bin.width = bin_width, k = k, nw = nw)
  
  
  
  # Spec of noise
  S_noise <- S_record
  S_noise$spec <- (S_record$spec - S_stack$spec) / (1 - 1 / n)
  
  # Spec of common signal (climate)
  S_clim <- S_record
  S_clim$spec <- S_record$spec - S_noise$spec
  
  
  if (isFALSE(prefilter) == FALSE) { #i.e. if TRUE
    S_record_f <- FilterSpec(S_record,
                             spans = prefilter,
                             method = 3, keep_low_f = TRUE
    )
    
    S_stack_f <- FilterSpec(S_stack,
                            spans = prefilter,
                            method = 3, keep_low_f = TRUE
    )
    
    S_clim_f <- FilterSpec(S_clim,
                           spans = prefilter,
                           method = 3, keep_low_f = TRUE
    )
    S_noise_f <- FilterSpec(S_noise,
                            spans = prefilter,
                            method = 3, keep_low_f = TRUE
    )
  }
  
  if (logsmooth) {
    S_record_f <- FilterSpecLog(
      S_record_f,
      df.log = 0.05
    )
    S_stack_f <- FilterSpecLog(
      S_stack_f,
      df.log = 0.05
    )
    S_clim_f <- FilterSpecLog(
      S_clim_f,
      df.log = 0.05
    )
    S_noise_f <- FilterSpecLog(
      S_noise_f,
      df.log = 0.05
    )
  }
  
  
  
  if (prefilter[1] == FALSE & logsmooth[1] == FALSE){
    S_SNR <- SpecRatio(S_clim, S_noise)
  } else {
    S_SNR <- SpecRatio(S_clim_f, S_noise_f)
  }
  
  
  # Add CI
  S_clim <- PaleoSpec::AddConfInterval(S_clim)
  S_noise <- PaleoSpec::AddConfInterval(S_noise)
  S_stack <- PaleoSpec::AddConfInterval(S_stack)
  S_proxy <- PaleoSpec::AddConfInterval(S_record)
  
  if (logsmooth == FALSE){
    spec_list <- list(
      S_noise = Spec2DF(S_noise),
      S_clim = Spec2DF(S_clim),
      S_proxy = Spec2DF(S_record),
      S_stack = Spec2DF(S_stack),
      SignalNoise = Spec2DF(S_SNR)
    )  
  } else {
    spec_list <- list(
      S_noise = Spec2DF(S_noise_f),
      S_clim = Spec2DF(S_clim_f),
      S_proxy = Spec2DF(S_record_f),
      S_stack = Spec2DF(S_stack_f),
      SignalNoise = Spec2DF(S_SNR)
    )
  }
  
  
  
  return(spec_list)
}

SpecRatio <- function(S1, S2){
  
  stopifnot(all(S1$freq == S2$freq))
  
  S3 <- S1
  
  S3$spec <- S1$spec / S2$spec
  S3$dof1 <- S1$dof
  S3$dof2 <- S2$dof
  
  S3 <- AddConfIntervalRatio2(S3)
  
  return(S3)
}



AddConfIntervalRatio <- function(spec, pval = 0.05, MINVALUE = 1e-10) {
  
  #is.spectrum(spec)
  
  spec$lim.1 <- qf(c(1 - pval / 2), spec$dof1, spec$dof2) * spec$spec
  spec$lim.2 <- qf(c(pval / 2), spec$dof1, spec$dof2) * spec$spec 
  
  spec$lim.1[spec$lim.1 < MINVALUE] <- MINVALUE
  spec$lim.2[spec$lim.2 < MINVALUE] <- MINVALUE
  
  return(spec)
  
}

AddConfIntervalRatio2 <- function(spec, pval = 0.05, MINVALUE = 1e-10) {
  
  #is.spectrum(spec)
  
  spec$lim.2 <- 1/qf(c(1 - pval / 2), spec$dof1, spec$dof2) * spec$spec
  spec$lim.1 <- 1/qf(c(pval / 2), spec$dof1, spec$dof2) * spec$spec 
  
  spec$lim.1[spec$lim.1 < MINVALUE] <- MINVALUE
  spec$lim.2[spec$lim.2 < MINVALUE] <- MINVALUE
  
  return(spec)
  
}


### plotting functions ------

#' Title
#'
#' @param x 
#' @param gg 
#' @param conf 
#' @param spec_id 
#' @param colour 
#' @param group 
#' @param alpha.line 
#' @param alpha.ribbon 
#' @param removeFirst 
#' @param removeLast 
#' @param min.colours 
#' @param force.lims 
#'
#' @return
#' @export
#'
#' @examples
gg_spec2 <- function(x, gg = NULL,
                     conf = TRUE,
                     spec_id = NULL,
                     colour = spec_id,
                     group = spec_id,
                     #linetype = spec_id,
                     alpha.line = 1,
                     alpha.ribbon = c(0.166, 0.333),
                     removeFirst = 0, removeLast = 0,
                     min.colours = 2,
                     force.lims = FALSE,
                     quantiles = FALSE) {
  
  gg_installed <- requireNamespace("ggplot2", quietly = TRUE)
  
  if (gg_installed == FALSE){
    stop("package ggplot2 is required to use gg_spec(). To install ggplot2, run install.packages(\"ggplot2\") from the console")
  }
  
  if (class(x)[1] != "spec_df" & "list" %in% class(x) == FALSE){
    x <- list(x)
    names(x) <- spec_id
  }
  
  if (class(x)[1] != "spec_df"){
    df <- Spec2DF(x)
  } else {
    df <- x
  }
  
  
  if (removeFirst > 0) {
    df <- df %>% 
      group_by({{group}}) %>% 
      filter(rank(df$freq) > removeFirst)
    
    #[rank(df$freq) > removeFirst, ]
  }
  
  if (removeLast > 0) {
    df <- df[rank(-df$freq) > removeLast,]
  }
  
  if (exists("spec_id", df) == FALSE){
    df$spec_id <- 1
  }
  
  
  if (is.numeric(df$spec_id)){
    df$spec_id <- as.character(df$spec_id)
  }
  
  
  
  if (is.null(gg)) {
    p <- ggplot2::ggplot(data = df, aes(group = {{group}}))
  } else {
    p <- gg
  }
  
  # rename to PSD so that y axis label can be overwritten later
  df$PSD <- df$spec
  df$Frequency <- df$freq
  
  df <- as.data.frame(df)
  
  if (conf == TRUE & exists("lim.1", df)){
    
    if (nrow(df) > 1e04 & force.lims == FALSE){
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
      
      p <- p +
        ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = lim.2,
                                                     #group = spec_id,
                                                     ymax = lim.1, fill = {{ colour }}),
                             alpha = alpha.ribbon[1], colour = NA)
    }
  }
  
  if (quantiles == TRUE & exists("X2.5.", df)){
    
    if (nrow(df) > 1e04 & force.lims == FALSE){
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
      
      p <- p +
        ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = `X2.5.`,
                                                     #group = spec_id,
                                                     ymax = `X97.5.`, fill = {{ colour }}),
                             alpha = alpha.ribbon[1], colour = NA)+
        ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = `X15.9.`,
                                                     #group = spec_id,
                                                     ymax = `X84.1.`, fill = {{ colour }}),
                             alpha = alpha.ribbon[2], colour = NA)
    }
  }
  
  p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Frequency, y = PSD,
                                                      #group = spec_id,
                                                      #linetype = {{ linetype }},
                                                      colour = {{ colour }}),
                              alpha = alpha.line
  ) +
    ggplot2::scale_x_continuous(trans = "log10",
                                sec.axis = ggplot2::sec_axis(~ 1/., "Timescale")) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::annotation_logticks(sides = "tlb") +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer("",
                                 type = "qual",
                                 palette = "Dark2",
                                 aesthetics = c("colour", "fill")) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_alpha()
  
  g <- ggplot2::ggplot_build(p)
  colrs <- unlist(unique(sapply(g$data, function(x) unique(x["colour"])$colour)))
  colrs <- colrs[is.na(colrs) == FALSE]
  ncolrs <- length(colrs)
  
  
  
  if (ncolrs <= min.colours){
    
    p <- p + ggplot2::scale_colour_manual("", values = "black",
                                          aesthetics = c("colour", "fill"))
    
    if (is.null({{ colour }})){
      p <- p + ggplot2::theme(legend.position = "none")
    }
    
  }
  
  p
}

## Utilities ----

summarise_q_2 <- function (dat, var, probs = 
                             c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)) 
{
  dat %>% dplyr::summarise(
    mean = mean({{ var }}, na.rm = TRUE),
    sd = stats::sd({{ var }}, na.rm = TRUE),
    n = sum(is.na({{ var }}) == FALSE),
    x = stats::quantile({{ var }}, probs, na.rm = TRUE), 
    q = paste0(round(100 * probs, 1), "%")
  ) %>% 
    tidyr::pivot_wider(names_from = .data$q, values_from = .data$x) %>% 
    dplyr::as_tibble()
}


