# Functions copied from the PaleoSpec package ------
# https://github.com/EarthSystemDiagnostics/paleospec
# Reproduced here for future reproducibility without dependency on PaleoSpec
# versioning.

## SpecACF ------
#' Estimate Power Spectra via the Autocovariance Function With Optional Slepian
#' Tapers
#'
#' @param x a vector or matrix of binned values, possibly with gaps
#' @param deltat,bin.width the time-step of the timeseries, equivalently the
#' width of the bins in a binned timeseries, set only one
#' @param pos.f.only return only positive frequencies, defaults to TRUE If TRUE,
#'   freq == 0, and frequencies higher than 1/(2*bin.width) which correspond to
#'   the negative frequencies are removed
#' @param demean remove the mean from each record (column) in x, defaults to
#'   TRUE. If detrend is TRUE, mean will be removed during detrending regardless
#'   of the value of demean
#' @param detrend remove the mean and any linear trend from each record (column)
#'   in x, defaults to FALSE
#' @description Estimates the power spectrum from a single time series, or the
#'   mean spectrum of a set of timeseries stored as the columns of a matrix.
#'   Timeseries can contain (some) gaps coded as NA values. Gaps results in
#'   additional estimation error so that the power estimates are no longer
#'   chi-square distributed and can contain additional additive error, to the
#'   extent that power at some frequencies can be negative. We do not have a
#'   full understanding of this estimation uncertainty, but simulation testing
#'   indicates that the estimates are unbiased such that smoothing across
#'   frequencies to remove negative estimates results in an unbiased power
#'   spectrum.
#' @inheritParams SpecMTM
#' @importFrom multitaper spec.mtm
#' @author Torben Kunz and Andrew Dolman <andrew.dolman@awi.de>
#' @return a spec object (list)
#' @family functions to estimate power spectra
#' @export
#'
#' @examples
#' set.seed(20230312)
#'
#' # Comparison with SpecMTM
#'
#' tsM <- replicate(2, SimPLS(1e03, 1, 0.1))
#' spMk3 <- SpecACF(tsM, bin.width = 1, k = 3, nw = 2)
#' spMk1 <- SpecACF(tsM, bin.width = 1, k = 1, nw = 0)
#'
#' spMTMa <- SpecMTM(tsM[,1], deltat = 1)
#' spMTMb <- SpecMTM(tsM[,2], deltat = 1)
#' spMTM <- spMTMa
#' spMTM$spec <- (spMTMa$spec + spMTMb$spec)/2
#'
#' gg_spec(list(
#'   `ACF k=1` = spMk1,
#'   `ACF k=3` = spMk3,
#'   `MTM k=3` = spMTM
#' ), alpha.line = 0.75) +
#'   ggplot2::facet_wrap(~spec_id)
#'
#' ## No gaps
#'
#' ts1 <- SimPLS(1000, 1, 0.1)
#'
#' sp_ACF1 <- SpecACF(ts1, 1, k = 1)
#' sp_MTM7 <- SpecMTM(ts1, nw = 4, k = 7, deltat = 1)
#' sp_ACF7 <- SpecACF(ts1, 1, k = 7, nw = 4)
#'
#' gg_spec(list(
#'   `ACF k=1` = sp_ACF1, `ACF k=7` = sp_ACF7, `MTM k=7` = sp_MTM7
#' ))
#'
#' # With Gaps
#'
#' gaps <- (arima.sim(list(ar = 0.5), n = length(ts1))) > 1
#' table(gaps)
#' ts1_g <- ts1
#' ts1_g[gaps] <- NA
#'
#' sp_ACF1_g <- SpecACF(ts1_g, 1)
#' sp_ACFMTM1_g <- SpecACF(ts1_g, bin.width = 1, nw = 4, k = 7)
#'
#' gg_spec(list(
#'   ACF_g = sp_ACF1_g,
#'   ACF_g_smoothed = FilterSpecLog(sp_ACF1_g),
#'   ACF_g_tapered = sp_ACFMTM1_g
#' ), conf = FALSE) +
#'   ggplot2::geom_abline(intercept = log10(0.1), slope = -1, lty = 2)
#'
#'
#'
#' ## AR4
#' arc_spring <- c(2.7607, -3.8106, 2.6535, -0.9238)
#'
#' tsAR4 <- arima.sim(list(ar = arc_spring), n = 1e03) + rnorm(1e03, 0, 10)
#' plot(tsAR4)
#' spAR4_ACF <- SpecACF(tsAR4, 1)
#' spAR4_MTACF <- SpecACF(as.numeric(tsAR4), 1, k = 15, nw = 8)
#'
#' gg_spec(list(#'
#'   `ACF k=1` = spAR4_ACF,
#'   `ACF k=15` = spAR4_MTACF)
#' )
#'
#' ## Add gaps to timeseries
#'
#' gaps <- (arima.sim(list(ar = 0.5), n = length(tsAR4))) > 2
#' table(gaps)
#' tsAR4_g <- tsAR4
#' tsAR4_g[gaps] <- NA
#'
#' plot(tsAR4, col = "green")
#' lines(tsAR4_g, col = "blue")
#'
#' table(tsAR4_g > 0, useNA = "always")
#'
#' spAR4_ACF_g <- SpecACF(as.numeric(tsAR4_g), 1)
#' spAR4_MTACF_g <- SpecACF(as.numeric(tsAR4_g), 1, nw = 8, k = 15)
#'
#' table(spAR4_ACF_g$spec < 0)
#' table(spAR4_MTACF_g$spec < 0)
#'
#' gg_spec(list(
#'   `ACF gaps k=1` = spAR4_ACF_g,
#'   `ACF gaps k = 15` = spAR4_MTACF_g,
#'   `ACF full k = 15` = spAR4_MTACF
#' )
#' )
SpecACF <- function(x,
                    deltat = NULL, bin.width = NULL,
                    k = 3, nw = 2,
                    demean = TRUE, detrend = TRUE,
                    TrimNA = TRUE,
                    pos.f.only = TRUE,
                    return.working = FALSE) {
  if (is.null(deltat) & is.null(bin.width) & is.ts(x) == FALSE) {
    stop("One of deltat or bin.width must be set")
  }
  
  # Convert ts to vector
  if (is.ts(x)) {
    d <- dim(x)
    dt_ts <- deltat(x)
    
    x <- as.vector(x)
    
    if (is.null(deltat) == FALSE) {
      if (dt_ts != deltat) stop("timeseries deltat does not match argument deltat")
    }
    
    if (is.null(bin.width) == FALSE) {
      if (dt_ts != bin.width) stop("timeseries deltat does not match argument bin.width")
    }
    
    if (is.null(deltat)) deltat <- dt_ts
    if (is.null(bin.width)) bin.width <- dt_ts
    
    if (is.null(d) == FALSE) {
      dim(x) <- d
    }
  }
  
  if (is.null(bin.width) & is.null(deltat) == FALSE) {
    bin.width <- deltat
  }
  
  if (is.null(bin.width) == FALSE & is.null(deltat)) {
    deltat <- bin.width
  }
  
  # Convert vector to matrix
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  
  if (is.data.frame(x) == TRUE) {
    x <- as.matrix(x)
  }
  
  if (TrimNA) {
    x <- TrimNA(x)
  }
  
  if (detrend) {
    i <- seq_along(x[, 1])
    x <- apply(x, 2, function(y) {
      stats::residuals(stats::lm(y ~ i, na.action = "na.exclude"))
    })
  }
  
  # remove mean from each record
  if (demean) {
    x <- x - colMeans(x, na.rm = TRUE)
  }
  
  lag <- (0:(nrow(x) - 1))
  
  ncolx <- ncol(x)
  
  ## Tapering
  if (k > 1) {
    n <- nrow(x)
    
    dpssIN <- multitaper:::dpss(n,
                                k = k, nw = nw,
                                returnEigenvalues = TRUE
    )
    dw <- dpssIN$v #* sqrt(bin.width)
    ev <- dpssIN$eigen
    
    x <- matrix(unlist(lapply(1:ncol(x), function(i) {
      lapply(1:ncol(dw), function(j) {
        x[, i] * dw[, j]
      })
    })), nrow = nrow(x), byrow = FALSE)
  }
  
  if (return.working == TRUE) {
    working <- PaleoSpec:::mean.ACF(x, return.working = TRUE)
    acf <- working[["acf"]]
  } else {
    acf <- PaleoSpec:::mean.ACF(x)
  }
  
  freq <- (lag / (nrow(x))) / bin.width
  
  if (k > 1){
    spec <- Re(stats::fft(acf)) * bin.width * nrow(x)
  } else {
    spec <- Re(stats::fft(acf)) * bin.width
  }
  
  
  
  
  rfreq <- max(freq)
  
  if (pos.f.only) {
    spec <- spec[freq > 0 & freq <= 1 / (2 * bin.width)]
    freq <- freq[freq > 0 & freq <= 1 / (2 * bin.width)]
  }
  
  out <- list(
    bin.width = bin.width,
    rfreq = rfreq,
    nrec = ncolx,
    lag = lag, acf = acf,
    freq = freq, spec = spec, f.length = 1,
    dof = rep(2 * ncolx * k, length(freq))
  )
  
  class(out) <- c("SpecACF", "spec")
  
  if (return.working) {
    out <- list(working = working, spec = out)
  }
  return(out)
}


#' Bin a Timeseries Preserving Empty Bins
#'
#' @param time.var time variable, a vector
#' @param value.var value variable, a vector
#' @param bin.width with of time bins in output
#' @param strt.time first time value in output, defaults to min(time.var) - bin.width/2
#' @param end.time last time values in output, defaults to max(time.var) + bin.width/2
#' @description Convert a continuous time series to a discrete (regular)
#'   timeseries by binning, preserving any empty bins. Additionally returns the
#'   number of data points in each bin
#' @return a dataframe
#' @export
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' set.seed(20230312)
#' x <- cumsum(rgamma(20, shape = 1.5, rate = 1.5/10))
#' y <- SimProxySeries(a = 0.1, b = 1, t.smpl = x)
#' BinTimeseries(x, y, bin.width = 15)
BinTimeseries <- function(time.var, value.var, bin.width,
                          strt.time = NULL, end.time = NULL){
  
  if (is.null(strt.time)){
    strt.time <-  round((min(time.var) - bin.width/2))
  }
  
  if (is.null(end.time)){
    end.time <-   round((max(time.var) + bin.width/2))
  }
  
  breaks <- seq(strt.time, end.time + bin.width/2, by = bin.width)
  
  time <- seq(strt.time + bin.width/2, end.time, by = bin.width)
  
  f <- cut(time.var, breaks = breaks)
  
  mean.value <- tapply(value.var, f, mean)
  mean.time <- tapply(time.var, f, mean)
  n.bin <- tapply(value.var, f, length)
  
  data.frame(time, bin = levels(f), mean.time, mean.value, n.bin, bin.width)
}


#' Title
#'
#' @param x a matrix
#' @param return.working return intermediate steps, defaults to FALSE
#'
#' @return a matrix
#' @keywords internal
mean.ACF <- function(x, return.working = FALSE){
  # Arguments:
  # x is a matrix, each column of which contains a timeseries (with missing values indicated by NA)
  
  # Value:
  # the circular ACF, with length(ACF)==dim(x)[1], being the length of the timeseries
  
  ix <- !is.na(x)
  x[is.na(x)] <- 0
  
  mvacfx <- mvacf.by.fft(x)
  mvacfix <- mvacf.by.fft(ix)
  
  sum.x.mvacf <- rowSums(mvacfx)
  sum.ix.mvacf <- rowSums(mvacfix)
  
  if (min(round(sum.ix.mvacf * dim(x)[1]))==0){
    warning("Some lags cannot be computed!")
    sum.ix.mvacf[round(sum.ix.mvacf * dim(x)[1])==0] <- NA
  }
  # This means the ACF could not be computed for at least one lag.
  # One could set ix.mvacf[round(ix.mvacf)==0] <- NA
  # but from an ACF with NAs one cannot compute the PSD via the FFT
  # unless one first interpolates the ACF across the missing lags.
  
  if (return.working == TRUE){
    out <- list(x = x, ix = ix, mvacfx = mvacfx, mvacfix = mvacfix,
                sum.x.mvacf = sum.x.mvacf, sum.ix.mvacf = sum.ix.mvacf,
                acf = sum.x.mvacf / sum.ix.mvacf)
  } else {
    out <- sum.x.mvacf / sum.ix.mvacf
  }
  
  return(out)
  
}


#' Title
#'
#' @param x a matrix
#'
#' @return a matrix
#' @keywords internal
mvacf.by.fft <- function(x){
  x.fft <- stats::mvfft(x)
  (Re(stats::mvfft(x.fft * Conj(x.fft), inverse=TRUE))) / dim(x)[1]^2
}

#' Remove leading and trailing rows of all NA
#'
#' @param m a numeric matrix, data.frame or vector
#' @param trim trim leading and trailing rows of "all" NA or containing "any" NA values
#' @return a numeric matrix, data.frame or vector
#' @export
#' @examples
#' m <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 1, 1, NA, NA, NA, 1:9, NA,NA,NA, 10:12, NA, 1, NA, NA,NA,NA), ncol = 3, byrow = TRUE)
#' m
#' TrimNA(m)
#' TrimNA(m, trim = "any")
TrimNA <- function(m, trim = c("all", "any")) {
  trim <- match.arg(trim)
  
  # make it work on vectors
  class_m <- class(m)
  if (class_m[1] == "numeric") {
    m <- cbind(m)
  }
  
  if (trim == "all") {
    empty.row <- is.nan(rowMeans(m, na.rm = TRUE))
    rank.good <- (empty.row == FALSE) * 1:length(empty.row)
  } else if (trim == "any") {
    empty.row <- is.na(rowMeans(m))
    rank.good <- (empty.row == FALSE) * 1:length(empty.row)
  }
  
  first.good <- which.min(empty.row * 1:length(empty.row))
  last.good <- which.max(rank.good)
  
  m <- m[first.good:last.good, , drop = FALSE]
  
  # return to a vector
  if (class_m[1] == "numeric") {
    m <- as.numeric(m)
  }
  
  return(m)
}

## Interpolated spectra ----

#' Interpolate spectrum
#'
#' Interpolate a spectrum to a given frequency resolution. The interpolation
#' includes the spectral estimates themselves as well as their degrees of
#' freedom.
#'
#' @param spec a spectral object of class \code{"spec"} or a list with the
#'   minimum components \code{freq}, \code{spec} and \code{dof} which are
#'   vectors of the same length giving the original frequency resolution and the
#'   corresponding spectral estimates and degrees of freedom.
#' @param freqRef numeric vector with the target frequency resolution.
#' @param check logical; if \code{TRUE} (the default) \code{spec} is checked for
#'   being a proper spectral object. This can be turned off in programmatic use
#'   when the respective check has already been performed.
#' @return a spectral object of class \code{"spec"} with the spectrum
#'   interpolated to the target resolution.
#' @author Thomas Laepple
#' @examples
#'
#' freqRef <- seq(0.1, 0.5, 0.1)
#' spec <- list(freq = c(0.1, 0.2, 0.4, 0.5), spec = c(1, 2, 4, 5),
#'              dof = rep(1, 4))
#'
#' SpecInterpolate(spec, freqRef)
#' @export
SpecInterpolate <- function(spec, freqRef, check = TRUE) {
  
  if (check) is.spectrum(spec)
  
  result <- list(
    freq = freqRef,
    spec = approx(spec$freq, spec$spec, freqRef)$y,
    dof  = approx(spec$freq, spec$dof, freqRef)$y
  )
  
  class(result) <- "spec"
  
  return(result)
  
}

## Integrate spectra -------
#' @title Variance estimate by integrating a part of the spectrum
#' @param spec spectrum (list of spec,freq,dof) to be analysed
#' @param f f[1],f[2]: frequency interval to be analysed
#' @param dfreq frequency discretisation used in the temporary  interpolation
#' @param df.log if > 0, smooth the spectra prior to integrating
#' @param bw the bandwidth assumed for the confinterval calculation (from the multitaper spectral estimate)
#' @return list(var,dof)  variance and corresponding dof
#' @author Thomas Laepple
#' @examples
#' x<-ts(rnorm(100))
#' spec<-SpecMTM(x)
#' var(x) #Sample variance of the timeseries
#' GetVarFromSpectra(spec,c(1/100,0.5))
#' GetVarFromSpectra(spec,c(0.25,0.5))
#' @export
GetVarFromSpectra <- function(spec,f,dfreq=NULL,df.log = 0,bw = 3)
{
  
  if (f[1] >= f[2]) stop("f1 must be < f2")
  freqVector <- spec$freq
  
  ## Test it both frequencies are included
  if (f[1] < FirstElement(freqVector)) {
    warning("f[1] is smaller than the lowest frequency in the spectrum, set to the lowest frequency")
    f[1] <- FirstElement(freqVector)
  }
  if (f[2] > LastElement(freqVector))  {
    warning("f[2] is larger than the highest frequency in the spectrum, set to the highest frequency")
    f[2] <- LastElement(freqVector)
  }
  
  if (is.null(dfreq)) dfreq <- min(diff(spec$freq)[1]/5, (f[2]-f[1])/100)
  newFreq <- seq(from = f[1],to = f[2],by = dfreq) # temporary frequency vector
  
  ## For spectra fromn the periodogram, the dof are supplied as a scalar named df; for MTN as a vector called DOF
  if (is.null(spec$dof)) spec$dof <- rep(spec$df,length(spec$freq))
  
  
  dof.original <- mean(SpecInterpolate(spec,newFreq)$dof)
  ## DOF before smoothing
  
  if (df.log > 0) spec <- LogSmooth(spec,removeLast = 0,df.log = df.log)
  
  spec.int <- SpecInterpolate(spec, newFreq)
  vars <- mean(spec.int$spec)*diff(f)*2
  dof <- mean(spec.int$dof)
  
  ## Estimate the DOF by calculating how many independent
  ## spectral estimates contribute to the calculated mean value
  dfreq <- mean(diff(spec$freq))
  nSpecEstimate <- (f[2]-f[1])/dfreq
  dof <- dof*nSpecEstimate/bw
  
  if (dof < dof.original) dof = dof.original
  
  return(list(var = vars,dof = dof))
}

## TidySpec -----
# tidy spectra classes, methods and functions

#' Make a spec Object
#'
#' @param x A named list with (at least) the elements freq, spec, and dof
#'
#' @return Object of class spec
#' @export
#' @family TidySpec
#' @author Andrew Dolman <andrew.dolman@awi.de>
as.spec <- function(x){
  class(x) <- unique(append(c("spec", "list"), class(x)) )
  x
}


#' Make a spec_df Object
#'
#' @param x A dataframe or tibble with (at least) columns freq, spec, and dof
#'
#' @return Object of class spec_df
#' @export
#' @family TidySpec
#' @author Andrew Dolman <andrew.dolman@awi.de>
as_spec_df <- function(x) {
  class(x) <- unique(append(c("spec_df", "list"), class(x)))
  x
}



#' Transform a spec Object Into a Dataframe
#'
#' @param x A spec object
#'
#' @return A dataframe or tibble (if package tibble is installed)
#' @export
#' @family TidySpec
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(rnorm(100))
#' sp1 <- SpecMTM(ts1)
#' sp1_df <- as.data.frame(sp1)
as.data.frame.spec <- function(x){
  
  df <- data.frame(
    freq = x$freq,
    spec = x$spec
  )
  
  if (exists("dof", x)){
    df$dof <- x$dof
  }
  
  if (exists("lim.1", x)){
    
    df$lim.1 <- x$lim.1
    df$lim.2 <- x$lim.2
    
  }
  
  tibble_installed <- requireNamespace("tibble", quietly = TRUE)
  
  if (tibble_installed){
    df <- tibble::as_tibble(df)
  }
  
  df <- as_spec_df(df)
  
  return(df)
  
}


#' Transform Spec Object(s) Into a Dataframe
#'
#' @param x A spec object or (optionally named) list of spec objects.
#'
#' @return A spec_df object
#' @export
#' @family TidySpec
#' @importFrom data.table rbindlist
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(SimPLS(1000, beta = 1))
#' sp1 <- SpecMTM(ts1)
#' sp1 <- AddConfInterval(sp1)
#' ts2 <- ts(rnorm(1000))
#' sp2 <- SpecMTM(ts2)
#' sp_lst <- list(sp1 = sp1, sp2 = sp2)
#' sp_df <- Spec2DF(sp_lst)
#' sp_df
Spec2DF <- function(x){
  
  if ("spec" %in% class(x)){
    x <- list(x)
  }
  
  df.lst <- lapply(x, as.data.frame.spec)
  
  df <- data.table::rbindlist(df.lst, fill = TRUE, idcol = "spec_id")
  
  tibble_installed <- requireNamespace("tibble", quietly = TRUE)
  
  if (tibble_installed){
    df <- tibble::as_tibble(df)
  }
  
  class(df) <- c("spec_df", class(df))
  
  return(df)
}


#' Transform a spec_df Object into a spec Object or List of spec Objects
#'
#' @param spec_df A spec_df object
#'
#' @return A spec object or a list of spec objects
#' @export
#' @family TidySpec
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(SimPLS(1000, beta = 1))
#' sp1 <- SpecMTM(ts1)
#' ts2 <- ts(rnorm(1000))
#' sp2 <- SpecMTM(ts2)
#' sp_lst <- list(sp1 = sp1, sp2 = sp2)
#' sp_df <- Spec2DF(sp_lst)
#' sp_df
#' str(DF2Spec(sp_df))
DF2Spec <- function(spec_df){
  
  stopifnot(c("freq", "spec") %in% names(spec_df))
  
  if ("spec_id" %in% names(spec_df)) {
    lst <- split(spec_df, spec_df$spec_id)
  } else {
    lst <- list(spec_df)
  }
  
  lst <- lapply(lst, function(x) as.spec(as.list(x)))
  
  
  if (length(lst) == 1) {
    lst <- lst[[1]]
    lst <- as.spec(lst)
  } else {
    class(lst) <- c("list", "spec_list", class(lst))
  }
  
  class(lst) <- unique(class(lst))
  
  
  return(lst)
}



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
#' @family functions to filter / smooth spectra
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
#' @family functions to filter / smooth spectra
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

## Plotting functions -------
#' Plot One or More Spectra with ggplot2
#'
#' @param x An object of class "spec" or "spec_df"
#' @param gg An existing ggplot object on which to add a new spec layer
#' @param spec_id Name for this plot layer
#' @param conf Plot shaded confidence interval if it exists in the spec object
#' @param alpha.line Alpha level for the spectra line(s)
#' @param alpha.ribbon Alpha level for the confidence region(s)
#' @param colour Name of variable map to colour, unquoted
#' @param linetype Name of variable map to linetype, unquoted
#' @param group Name of variable to map to group, unquoted
#' @param min.colours Minimum number of spectra before starting to colour them
#' separately
#' @param removeFirst,removeLast Remove first or last "n" values on the low or
#'   high frequency side respectively. Operates on a per group basis.
#' @param force.lims,force.CI Force the plotting of confidence regions when the total number
#' of frequencies exceeds 10000. force.CI is deprecated, use force.lims. Defaults to FALSE
#' @param quantiles Plot uncertainty quantiles if they exist in the spec object
#' @param time_unit Optional string giving the time unit for axes labels, e.g. "years"
#'
#' @return a ggplot object
#' @family functions to plot power spectra
#' @export
#'
#' @examples
#' library(PaleoSpec)
#' N <- 1e03
#' beta <- 1
#' alpha <- 0.1
#'
#' ts1 <- SimPLS(N = N, b = beta, a = alpha)
#' ts2 <- SimPLS(N = N, b = beta, a = alpha)
#' sp1 <- SpecMTM(ts1, deltat = 1)
#' sp1 <- AddConfInterval(sp1)
#' sp2 <- SpecMTM(ts2, deltat = 1)
#'
#' # plot single spectrum
#' p <- gg_spec(sp1, spec_id = "df1")
#' p
#'
#' # Add additional second spectra
#' p <- gg_spec(sp2, p, spec_id = "df2", removeFirst = 2)
#' p
#'
#' p <- gg_spec(sp1, p, spec_id = "df3", removeLast = 200)
#' p
#'
#' sp2 <- LogSmooth(sp1)
#' p <- gg_spec(sp1, spec_id = "df1")
#' p <- gg_spec(sp2, p, spec_id = "df2")
#' p <- p + ggplot2::geom_abline(intercept = log10(alpha), slope = -beta, colour = "red")
#' p
#'
#' # Or directly plot named or unnamed list
#'
#' gg_spec(list(sp1, sp2))
#' gg_spec(list(raw = sp1, smoothed = sp2))
#'
#' # without setting any names all spectra will be black
#' p <- gg_spec(sp1)
#' sp2 <- LogSmooth(sp1)
#' p <- gg_spec(sp2, p)
#' p
gg_spec <- function(x, gg = NULL,
                    conf = TRUE,
                    spec_id = NULL,
                    colour = spec_id,
                    group = spec_id,
                    linetype = NULL,
                    alpha.line = 1,
                    alpha.ribbon = c(0.166, 0.333),
                    removeFirst = 0, removeLast = 0,
                    min.colours = 2,
                    force.lims = FALSE,
                    force.CI = NULL,
                    quantiles = FALSE,
                    time_unit = NULL) {
  if (!missing("force.CI")) {
    
    warning("the force.CI argument is deprecated please use force.lims")
    force.lims <- force.CI
  }
  
  gg_installed <- requireNamespace("ggplot2", quietly = TRUE)
  
  if (gg_installed == FALSE) {
    stop("package ggplot2 is required to use gg_spec(). To install ggplot2, run install.packages(\"ggplot2\") from the console")
  }
  
  if (class(x)[1] != "spec_df" & "list" %in% class(x) == FALSE) {
    x <- list(x)
    names(x) <- spec_id
  }
  
  if (class(x)[1] != "spec_df") {
    df <- Spec2DF(x)
  } else {
    df <- x
  }
  
  
  if (is.numeric(df$spec_id)){
    df$spec_id <- as.character(df$spec_id)
  }
  
  
  if (removeFirst > 0) {
    # Convert 'df' to a data.table
    data.table::setDT(df)
    
    # Group by 'spec_id' and filter rows where 'rank(freq)' is greater than 'removeFirst'
    df <- df[, .SD[rank(freq) > removeFirst], by = spec_id]
    
    # return to data.frame
    df <- as.data.frame(df)
  }
  
  if (removeLast > 0) {
    # Convert 'df' to a data.table
    data.table::setDT(df)
    
    # Group by 'spec_id' and filter rows where 'rank(-freq)' is greater than 'removeLast'
    df <- df[, .SD[rank(-freq) > removeLast], by = spec_id]
    
    # return to data.frame
    df <- as.data.frame(df)
    
  }
  
  if (is.null(time_unit) == FALSE) {
    lab_timescale <- paste0("Timescale [", time_unit, "]")
    lab_freq <- paste0("Frequency [1/", time_unit, "]")
  } else {
    lab_timescale <- "Timescale"
    lab_freq <- "Frequency"
  }
  
  
  if (is.null(gg)) {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(group = {{ group }}))
  } else {
    p <- gg
  }
  
  # rename to PSD so that y axis label can be overwritten later
  df$PSD <- df$spec
  df$Frequency <- df$freq
  
  df <- as.data.frame(df)
  
  if (conf == TRUE & exists("lim.1", df)) {
    if (nrow(df) > 1e04 & force.lims == FALSE) {
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
      p <- p +
        ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = lim.2,
            ymax = lim.1, fill = {{ colour }}
          ),
          alpha = alpha.ribbon[1], colour = NA
        )
    }
  }
  
  
  if (quantiles == TRUE & exists("X2.5.", df)) {
    if (nrow(df) > 1e04 & force.lims == FALSE) {
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
      p <- p +
        ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = `X2.5.`,
            ymax = `X97.5.`, fill = {{ colour }}),
          alpha = alpha.ribbon[1], colour = NA
        ) +
        ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = `X15.9.`,
            ymax = `X84.1.`, fill = {{ colour }}
          ),
          alpha = alpha.ribbon[2], colour = NA
        )
    }
  }
  
  p <- p + ggplot2::geom_line(
    data = df, ggplot2::aes(
      x = Frequency, y = PSD,
      linetype = {{ linetype }},
      colour = {{ colour }}
    ),
    alpha = alpha.line
  ) +
    ggplot2::scale_x_continuous(
      name = lab_freq, trans = "log10",
      sec.axis = ggplot2::sec_axis(~ 1 / ., lab_timescale, labels = scales::comma)
    ) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::annotation_logticks(sides = "tlb") +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer("",
                                 type = "qual",
                                 palette = "Dark2",
                                 aesthetics = c("colour", "fill")
    ) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_alpha() +
    ggplot2::labs(linetype = "")
  
  g <- ggplot2::ggplot_build(p)
  colrs <- unlist(unique(sapply(g$data, function(x) unique(x["colour"])$colour)))
  colrs <- colrs[is.na(colrs) == FALSE]
  ncolrs <- length(colrs)
  
  
  if (ncolrs <= min.colours) {
    p <- p + ggplot2::scale_colour_manual("",
                                          values = "black",
                                          aesthetics = c("colour", "fill")
    )
    
    # remove legend if only 1 spec and no name given
    if (({{ deparse(substitute(colour))[[1]] }}) == "spec_id" &
        ({{ deparse(substitute(linetype))[[1]] }}) == "NULL"){
      if (is.null({{ colour }}) & is.null({{ linetype }})) {
        p <- p + ggplot2::theme(legend.position = "none")
      }
    }
  }
  p
}


## Utility functions ---------
# Number of spectral observations
#
# @param x spectral object
# @return numerical value.
# @author Thomas Münch
get.length <- function(x) {
  length(x$freq)
}
# Extract frequency vector
#
# @param x spectral object
# @return numerical vector of frequency axis.
# @author Thomas Münch
get.freq <- function(x) {
  x$freq
}
# Extract PSD vector
#
# @param x spectral object
# @return numerical vector of power spectral density values.
# @author Thomas Münch
get.spec <- function(x) {
  x$spec
}
# Extract DOF vector
#
# @param x spectral object
# @return numerical vector of degrees of freedom.
# @author Thomas Münch
get.dofs <- function(x) {
  x$dof
}

FirstElement<-function (x)
{
  return(x[1])
}

LastElement<-function (x)
{
  return(x[length(x)])
}


AddConfInterval <- function(spec, pval = 0.05, MINVALUE = 1e-10) {
  
  is.spectrum(spec)
  
  spec$lim.1 <- spec$spec * qchisq(c(1 - pval / 2), spec$dof) / (spec$dof)
  spec$lim.2 <- spec$spec * qchisq(c(pval / 2), spec$dof) / (spec$dof)
  spec$lim.1[spec$lim.1 < MINVALUE] <- MINVALUE
  spec$lim.2[spec$lim.2 < MINVALUE] <- MINVALUE
  
  spec$pval <- pval
  
  return(spec)
  
}

#' Check for spectral object
#'
#' In \code{PaleoSpec} a "spectral object" is defined to be a list consisting as
#' the minimum requirement of the three elements \code{freq}, \code{spec}, and
#' \code{dof}, which are numeric vectors of the same length giving the frequency
#' axis, the power spectral density, and the degrees of freedom of the spectral
#' estimate. Most functions in this package require such an object as their main
#' input. This function checks an object for the minimum requirements.
#'
#' @param x an object to be checked.
#' @param check.only logical; the default means to issue an error when a
#'   requirement is not met. If set to \code{TRUE}, all requirements are
#'   checked, and the function returns \code{TRUE} if they all pass, or
#'   \code{FALSE} if one or more checks fail.
#' @param dof logical; per default, \code{x} is also checked for containing a
#'   \code{dof} vector, which can be turned off here when unnecessary.
#' @return a logical value for \code{check.only = TRUE}.
#' @author Thomas Münch
is.spectrum <- function(x, check.only = FALSE, dof = TRUE) {
  
  isList <- is.list(x)
  
  if (isList) {
    
    hasFreq <- !is.null(get.freq(x))
    hasSpec <- !is.null(get.spec(x))
    hasDofs <- ifelse(dof, !is.null(get.dofs(x)), TRUE)
    
    hasAll <- all(hasFreq, hasSpec, hasDofs)
    
    if (hasAll) {
      lengths <- c(get.length(x), length(get.spec(x)))
      if (dof) lengths <- c(lengths, length(get.dofs(x)))
      hasEqualLength <- stats::var(lengths) == 0
      
    } else {
      hasEqualLength <- FALSE
    }
  }
  
  if (check.only) {
    
    if (isList) {
      hasAll & hasEqualLength
    } else {
      isList
    }
    
  } else {
    
    if (!isList) {
      stop("Passed argument is not a spectral list object.", call. = FALSE)
    }
    if (!hasFreq) {
      stop("Passed object has no frequency vector.", call. = FALSE)
    }
    if (!hasSpec) {
      stop("Passed object has no spectral density vector.", call. = FALSE)
    }
    if (!hasDofs) {
      stop("Passed object has no dof vector.", call. = FALSE)
    }
    if (!hasEqualLength) {
      stop("Frequency, PSD and DOF vectors have different lengths.",
           call. = FALSE)
    }
  }
  
}


# Additional custom functions -------
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

AddConfIntervalRatio2 <- function(spec, pval = 0.05, MINVALUE = 1e-10) {
  
  #is.spectrum(spec)
  
  spec$lim.2 <- 1/qf(c(1 - pval / 2), spec$dof1, spec$dof2) * spec$spec
  spec$lim.1 <- 1/qf(c(pval / 2), spec$dof1, spec$dof2) * spec$spec 
  
  spec$lim.1[spec$lim.1 < MINVALUE] <- MINVALUE
  spec$lim.2[spec$lim.2 < MINVALUE] <- MINVALUE
  
  return(spec)
  
}


#' ### plotting functions ------
#' 
#' #' Title
#' #'
#' #' @param x 
#' #' @param gg 
#' #' @param conf 
#' #' @param spec_id 
#' #' @param colour 
#' #' @param group 
#' #' @param alpha.line 
#' #' @param alpha.ribbon 
#' #' @param removeFirst 
#' #' @param removeLast 
#' #' @param min.colours 
#' #' @param force.lims 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' gg_spec2 <- function(x, gg = NULL,
#'                      conf = TRUE,
#'                      spec_id = NULL,
#'                      colour = spec_id,
#'                      group = spec_id,
#'                      #linetype = spec_id,
#'                      alpha.line = 1,
#'                      alpha.ribbon = c(0.166, 0.333),
#'                      removeFirst = 0, removeLast = 0,
#'                      min.colours = 2,
#'                      force.lims = FALSE,
#'                      quantiles = FALSE) {
#'   
#'   gg_installed <- requireNamespace("ggplot2", quietly = TRUE)
#'   
#'   if (gg_installed == FALSE){
#'     stop("package ggplot2 is required to use gg_spec(). To install ggplot2, run install.packages(\"ggplot2\") from the console")
#'   }
#'   
#'   if (class(x)[1] != "spec_df" & "list" %in% class(x) == FALSE){
#'     x <- list(x)
#'     names(x) <- spec_id
#'   }
#'   
#'   if (class(x)[1] != "spec_df"){
#'     df <- Spec2DF(x)
#'   } else {
#'     df <- x
#'   }
#'   
#'   
#'   if (removeFirst > 0) {
#'     df <- df %>% 
#'       group_by({{group}}) %>% 
#'       filter(rank(df$freq) > removeFirst)
#'     
#'     #[rank(df$freq) > removeFirst, ]
#'   }
#'   
#'   if (removeLast > 0) {
#'     df <- df[rank(-df$freq) > removeLast,]
#'   }
#'   
#'   if (exists("spec_id", df) == FALSE){
#'     df$spec_id <- 1
#'   }
#'   
#'   
#'   if (is.numeric(df$spec_id)){
#'     df$spec_id <- as.character(df$spec_id)
#'   }
#'   
#'   
#'   
#'   if (is.null(gg)) {
#'     p <- ggplot2::ggplot(data = df, aes(group = {{group}}))
#'   } else {
#'     p <- gg
#'   }
#'   
#'   # rename to PSD so that y axis label can be overwritten later
#'   df$PSD <- df$spec
#'   df$Frequency <- df$freq
#'   
#'   df <- as.data.frame(df)
#'   
#'   if (conf == TRUE & exists("lim.1", df)){
#'     
#'     if (nrow(df) > 1e04 & force.lims == FALSE){
#'       warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
#'     } else {
#'       
#'       p <- p +
#'         ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = lim.2,
#'                                                      #group = spec_id,
#'                                                      ymax = lim.1, fill = {{ colour }}),
#'                              alpha = alpha.ribbon[1], colour = NA)
#'     }
#'   }
#'   
#'   if (quantiles == TRUE & exists("X2.5.", df)){
#'     
#'     if (nrow(df) > 1e04 & force.lims == FALSE){
#'       warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
#'     } else {
#'       
#'       p <- p +
#'         ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = `X2.5.`,
#'                                                      #group = spec_id,
#'                                                      ymax = `X97.5.`, fill = {{ colour }}),
#'                              alpha = alpha.ribbon[1], colour = NA)+
#'         ggplot2::geom_ribbon(data = df, ggplot2::aes(x = Frequency, ymin = `X15.9.`,
#'                                                      #group = spec_id,
#'                                                      ymax = `X84.1.`, fill = {{ colour }}),
#'                              alpha = alpha.ribbon[2], colour = NA)
#'     }
#'   }
#'   
#'   p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Frequency, y = PSD,
#'                                                       #group = spec_id,
#'                                                       #linetype = {{ linetype }},
#'                                                       colour = {{ colour }}),
#'                               alpha = alpha.line
#'   ) +
#'     ggplot2::scale_x_continuous(trans = "log10",
#'                                 sec.axis = ggplot2::sec_axis(~ 1/., "Timescale")) +
#'     ggplot2::scale_y_continuous(trans = "log10") +
#'     ggplot2::annotation_logticks(sides = "tlb") +
#'     ggplot2::theme_bw() +
#'     ggplot2::scale_colour_brewer("",
#'                                  type = "qual",
#'                                  palette = "Dark2",
#'                                  aesthetics = c("colour", "fill")) +
#'     ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
#'     ggplot2::scale_alpha()
#'   
#'   g <- ggplot2::ggplot_build(p)
#'   colrs <- unlist(unique(sapply(g$data, function(x) unique(x["colour"])$colour)))
#'   colrs <- colrs[is.na(colrs) == FALSE]
#'   ncolrs <- length(colrs)
#'   
#'   
#'   
#'   if (ncolrs <= min.colours){
#'     
#'     p <- p + ggplot2::scale_colour_manual("", values = "black",
#'                                           aesthetics = c("colour", "fill"))
#'     
#'     if (is.null({{ colour }})){
#'       p <- p + ggplot2::theme(legend.position = "none")
#'     }
#'     
#'   }
#'   
#'   p
#' }

## Utilities ----

summarise_q_2 <- function (dat, var, probs = 
                             c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)) 
{
  dat %>% dplyr::reframe(
    mean = mean({{ var }}, na.rm = TRUE),
    sd = stats::sd({{ var }}, na.rm = TRUE),
    n = sum(is.na({{ var }}) == FALSE),
    x = stats::quantile({{ var }}, probs, na.rm = TRUE), 
    q = paste0(round(100 * probs, 1), "%")
  ) %>% 
    tidyr::pivot_wider(names_from = .data$q, values_from = .data$x) %>% 
    dplyr::as_tibble()
}

#' Interpolate Over a Limited Distance Leaving Gaps
#'
#' @param x
#' @param y
#' @param xout
#' @param method
#' @param rule
#' @param max_dx
#' @param return
#'
#' @return
#' @export
#'
#' @examples
ApproxABit <- function(x, y, xout,
                       method = "linear",
                       rule = 2, #ties = mean,
                       max_dx = NULL,
                       return = c("list", "dataframe")){
  
  return <- match.arg(return)
  
  o <- order(x)
  x <- x[o]
  y <- y[o]
  
  if (is.null(max_dx)){
    max_dx <- median(diff(x))
  }
  
  xout <- sort(xout)
  
  yout <- approx(x, y, xout, method = method, rule = rule#, ties = ties
  )$y
  
  # get nearest x to xout
  x_exists <- x[is.na(y) == FALSE]
  dxout <- abs(xout - sapply(xout, function(i) x_exists[which.min(abs(x_exists-i))]))
  yout[dxout > max_dx] <- NA
  
  if (return == "list"){
    out <- list(x = xout, y = yout)
  } else if (return == "dataframe"){
    out <- data.frame(x = xout, y = yout)
  }
  return(out)
}

