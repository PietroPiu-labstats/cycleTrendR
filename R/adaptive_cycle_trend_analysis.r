ACF <- NULL

#' Adaptive Trend and Cycle Analysis for Time Series
#'
#' @description
#' Performs adaptive trend estimation, cycle detection, Fourier harmonic
#' selection, bootstrap confidence intervals, change points detection, and
#' rolling-origin forecasting. Supports LOESS, GAM, and GAMM models, and
#' automatically handles irregular sampling using the Lomb Scargle periodogram.
#'
#' @importFrom mgcv gam
#' @importFrom lomb lsp
#' @importFrom stats coef Box.test acf as.formula gaussian lm loess median predict qchisq quantile setNames spec.pgram stl ts
#' @importFrom changepoint cpts cpt.meanvar
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs geom_ribbon geom_vline theme_minimal annotate geom_col geom_hline scale_color_manual ggsave
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom officer body_add_img body_add_par read_docx
#' @importFrom flextable body_add_flextable flextable
#' @importFrom tseries adf.test kpss.test
#' @importFrom utils tail
#'
#' @param signal Numeric vector of observed values.
#' @param dates Date vector of the same length as `signal`.
#' @param normalize Logical; if TRUE, Z score normalization is applied.
#' @param trendmethod Character; `"loess"` or `"gam"`.
#' @param usefourier Logical; whether to include Fourier harmonics.
#' @param fourierK Integer; fixed number of harmonics if auto selection disabled.
#' @param auto_fourier_select Logical; if TRUE, selects K via AICc/BIC.
#' @param fourier_selection_criterion `"AICc"` or `"BIC"`.
#' @param fourierK_max Maximum K to consider during selection.
#' @param cimethod `"model"`, `"bootstrapiid"`, or `"bootstrapmbb"`.
#' @param nboot Number of bootstrap samples.
#' @param blocksize Block size for MBB bootstrap.
#' @param seasonalfrequency Seasonal frequency for STL (regular sampling).
#' @param stlrobust Logical; robust STL decomposition.
#' @param specspans Smoothing spans for spectral estimation.
#' @param auto_seasonality Logical; if TRUE, uses dominant period.
#' @param lagmax Maximum lag for ACF and Ljung Box tests.
#' @param loess_span_mode `"auto_aicc"`, `"auto_gcv"`, `"cv"`, `"fixed"`.
#' @param loess_span_fixed Numeric; fixed LOESS span.
#' @param loess_span_grid Grid of spans for CV.
#' @param loess_cv_k Number of folds for blocked CV.
#' @param blocklength_mode `"auto_pwsd"`, `"heuristic"`, `"fixed"`.
#' @param blocklength_fixed Fixed block length.
#' @param robust Logical; robust LOESS or robust GAM family.
#' @param use_gamm Logical; fit GAMM instead of GAM.
#' @param group_var Character; grouping variable for random intercepts.
#' @param group_values Optional vector to attach as grouping variable.
#' @param random_effect Optional random effects list for `mgcv::gamm`.
#' @param cor_struct `"none"`, `"ar1"`, `"arma"`.
#' @param arma_p,arma_q ARMA orders.
#' @param forecast_holdout_h Holdout horizon for forecasting.
#' @param forecast_origin_mode `"expanding"` or `"sliding"`.
#' @param train_window Training window for sliding origin.
#' @param forecast_lock_K Logical; lock Fourier K across origins.
#' @param exportdocx Logical; export DOCX report.
#' @param exportplot Logical; export plots.
#' @param outputpath Path for DOCX export.
#' @param logopath Optional logo for DOCX.
#' @param logowidth,logoheight Logo dimensions.
#' @param project_id,cohort_id,assay_version,analyst,run_date,notes Metadata.
#' @param include_parameters_appendix Logical; include appendix in DOCX.
#'
#' @return A list containing:
#' \itemize{
#'   \item Trend estimates
#'   \item Confidence intervals
#'   \item Residuals and diagnostics
#'   \item Fourier selection results
#'   \item Change-point locations
#'   \item Spectral analysis
#'   \item Forecast results (if enabled)
#'   \item ggplot2 objects for visualization
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dates <- as.Date("2020-01-01") + cumsum(sample(1:3, 300, replace = TRUE))
#' signal <- sin(2*pi*as.numeric(dates)/20) + rnorm(300, 0, 0.3)
#'
#' res_gam <- adaptive_cycle_trend_analysis(
#'   signal = signal,
#'   dates = dates,
#'   usefourier = TRUE,
#'   trendmethod = "gam"
#' )
#'
#' dates <- as.Date("2020-01-01") + cumsum(sample(1:3, 150, replace = TRUE))
#' signal <- sin(2*pi*as.numeric(dates)/25) + rnorm(150, 0, 0.3)
#' group  <- rep(letters[1:4], length.out = length(signal))
#'
#' res_gamm <- adaptive_cycle_trend_analysis(
#'   signal = signal,
#'   dates = dates,
#'   trendmethod = "gam",
#'   use_gamm = TRUE,
#'   group_var = "subject",
#'   group_values = group,
#'   usefourier = FALSE,
#'   nboot = 20
#' )
#'
#' plot(res_gamm$Plot$Trend)
#'
#' dates <- as.Date("2020-01-01") + 1:120
#' signal <- sin(2*pi*(1:120)/20) + rnorm(120, 0, 0.2)
#'
#' res_loess <- adaptive_cycle_trend_analysis(
#'   signal = signal,
#'   dates = dates,
#'   trendmethod = "loess",
#'   usefourier = TRUE,
#'   auto_fourier_select = TRUE,
#'   nboot = 50
#' )
#'
#' plot(res_loess$Plot$Trend)
#' }
#'
#' @export

adaptive_cycle_trend_analysis <- function(

  signal,
  dates,
  normalize = FALSE,
  trendmethod = c("loess", "gam"),
  # Fourier options
  usefourier = FALSE,
  fourierK = 2,                       # legacy default if auto selection disabled
  auto_fourier_select = TRUE,
  fourier_selection_criterion = c("AICc", "BIC"),
  fourierK_max = 6,
  # CI options
  cimethod = c("model", "bootstrapiid", "bootstrapmbb"),
  nboot = 1000,
  blocksize = NULL,
  # Seasonality and spectra
  seasonalfrequency = 7,
  stlrobust = TRUE,
  specspans = c(7, 7),
  auto_seasonality = TRUE,
  # ACF/Ljung-Box
  lagmax = NULL,
  # LOESS span selection controls
  loess_span_mode = c("auto_aicc", "auto_gcv", "cv", "fixed"),
  loess_span_fixed = NULL,
  loess_span_grid = seq(0.15, 0.60, by = 0.05),
  loess_cv_k = 5,
  # Block-length selection controls for MBB
  blocklength_mode = c("auto_pwsd", "heuristic", "fixed"),
  blocklength_fixed = NULL,
  # Robust trend
  robust = TRUE,
  # ---- GAMM controls ----
  use_gamm = FALSE,                  # enable GAMM via mgcv::gamm
  group_var = NULL,                  # convenience: random intercept per group (character scalar)
  group_values = NULL,               # NEW: vector aligned to signal to attach as df[[group_var]] if missing
  random_effect = NULL,              # explicit random = list(group = ~1, ...)
  cor_struct = c("none", "ar1", "arma"),
  arma_p = 1, arma_q = 0,           # ARMA orders if cor_struct == "arma"
  # ---- Forecast hold-out & rolling-origin ----
  forecast_holdout_h = 0,
  forecast_origin_mode = c("expanding", "sliding"),
  train_window = NULL,               # for sliding; default to N - H
  forecast_lock_K = TRUE,            # lock selected K across origins
  # Export
  exportdocx = FALSE,
  exportplot = FALSE,
  outputpath = "analysisreport.docx",
  # DOCX logo
  logopath = NULL,
  logowidth = 1.5,
  logoheight = 0.5,
  # Report metadata
  project_id = NULL,
  cohort_id = NULL,
  assay_version = NULL,
  analyst = NULL,
  run_date = Sys.Date(),
  notes = NULL,
  include_parameters_appendix = TRUE
) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  trendmethod <- match.arg(trendmethod)
  cimethod    <- match.arg(cimethod)
  loess_span_mode <- match.arg(loess_span_mode)
  blocklength_mode <- match.arg(blocklength_mode)
  forecast_origin_mode <- match.arg(forecast_origin_mode)
  fourier_selection_criterion <- match.arg(fourier_selection_criterion)
  cor_struct <- match.arg(cor_struct)

  # ---- Validate ----
  if (!is.numeric(signal)) stop("signal must be numeric.")
  if (!inherits(dates, "Date")) stop("dates must be of class Date.")
  if (length(signal) != length(dates)) stop("signal and dates must have same length.")
  if (usefourier && fourierK < 0) stop("fourierK must be >= 0.")
  N <- length(signal)
  if (forecast_holdout_h < 0 || forecast_holdout_h > max(0, N - 5)) {
    stop("forecast_holdout_h must be between 0 and N-5.")
  }
  if (use_gamm && !requireNamespace("nlme", quietly = TRUE)) {
    stop("use_gamm=TRUE requires the 'nlme' package. Please install.packages('nlme').")
  }
  if (use_gamm && identical(trendmethod, "loess")) {
    message("Note: use_gamm=TRUE is ignored because trendmethod='loess'. Set trendmethod='gam' to enable GAMM.")
  }

  # ---- Normalize ----
  if (normalize) {
    signal <- as.numeric(scale(signal))
    message("Data normalized using Z-score standardization.")
  } else {
    message("Using original scale.")
  }

  # ---- Regularity ----
  datediff <- diff(as.numeric(dates))
  irregular <- if (length(datediff) > 0) any(abs(datediff - median(datediff)) > 1e-6) else FALSE
  message(ifelse(irregular, "Irregular dates detected: Lomb Scargle + LOESS/GAM(M)",
                 "Regular dates: STL + LOESS/GAM(M)"))

  # ---- Data & time index ----
  df <- data.frame(Date = dates, Value = signal)
  timenum <- as.numeric(df$Date)
  df$timenum <- timenum
  avgintervaldays <- if (length(datediff)) mean(datediff) else 1

  # ---- Attach/convert group_var ----
  if (use_gamm && !is.null(group_var)) {
    # Ensure group_var is a proper character scalar
    if (!is.character(group_var) || length(group_var) != 1L || nchar(group_var) == 0L) {
      stop("`group_var` must be a non-empty character scalar naming a column.")
    }
    # Attach group_values if df lacks the column
    if (!group_var %in% names(df)) {
      if (!is.null(group_values)) {
        if (length(group_values) != nrow(df)) {
          stop(sprintf("`group_values` length (%d) must match data length (%d).",
                       length(group_values), nrow(df)))
        }
        df[[group_var]] <- group_values
        message(sprintf("Attached grouping column '%s' from `group_values`.", group_var))
      }
    }
    # Column must exist now
    if (!group_var %in% names(df)) {
      stop(sprintf("group_var '%s' not found in the data.", group_var))
    }
    # Convert to factor if needed
    if (!is.factor(df[[group_var]])) {
      df[[group_var]] <- as.factor(df[[group_var]])
      message(sprintf("Converted group_var '%s' to factor (%d levels).",
                      group_var, nlevels(df[[group_var]])))
    }
  }

  # ---- Dominant period estimation (fixed irregular branch) ----
  estimate_dominant <- function(x, dates, irregular, specspans, avgintervaldays) {
    if (!irregular) {
      spec <- spec.pgram(x, spans = specspans, log = "no", taper = 0.1, detrend = TRUE)
      freq <- spec$freq; power <- spec$spec
      if (length(power) > 0) {
        ix <- which.max(power); dfreq <- freq[ix]
        list(freq = dfreq,
             period_samples = if (!is.na(dfreq) && dfreq > 0) 1 / dfreq else NA_real_,
             period_days = if (!is.na(dfreq) && dfreq > 0) (1 / dfreq) * avgintervaldays else NA_real_)
      } else {
        list(freq = NA_real_, period_samples = NA_real_, period_days = NA_real_)
      }
    } else {
      tnum <- as.numeric(dates - min(dates))
      lombres <- try(lsp(x, tnum, type = "frequency", plot = FALSE), silent = TRUE)
      if (inherits(lombres, "try-error") ||
          is.null(lombres$peak.at) || length(lombres$peak.at) < 2 ||
          any(!is.finite(lombres$peak.at))) {
        return(list(freq = NA_real_, period_samples = NA_real_, period_days = NA_real_))
      }
      dfreq <- lombres$peak.at[1]
      p_days <- lombres$peak.at[2]
      list(
        freq = dfreq,
        period_samples = if (!is.na(dfreq) && dfreq > 0) 1 / dfreq else NA_real_,
        period_days = p_days
      )
    }
  }
  dom <- estimate_dominant(signal, dates, irregular, specspans, avgintervaldays)

  # ---- Seasonality from dominant period ----
  if (auto_seasonality && !is.na(dom$period_samples) && dom$period_samples > 1) {
    seasonalfrequency <- max(2L, round(dom$period_samples))
    message(sprintf("Auto seasonality: dominant period approx %.2f samples -> seasonalfrequency=%d",
                    dom$period_samples, seasonalfrequency))
  } else {
    message(sprintf("Using provided seasonalfrequency=%d", seasonalfrequency))
  }
  perioddays <- if (!is.na(dom$period_days)) dom$period_days else seasonalfrequency * avgintervaldays

  # ---- Fourier design ----
  makefourier <- function(timenum, perioddays, K) {
    phi <- 2 * pi * timenum / perioddays
    out <- data.frame(row.names = seq_along(timenum))
    if (K >= 1) for (k in seq_len(K)) {
      out[[paste0("sin", k)]] <- sin(k * phi)
      out[[paste0("cos", k)]] <- cos(k * phi)
    }
    out
  }
  maxK <- if (usefourier) max(fourierK, fourierK_max) else 0
  full_fourier_design <- if (usefourier && maxK > 0) makefourier(timenum, perioddays, maxK) else NULL

  # ---- Helper: Blocked k-fold CV for LOESS span ----
  blocked_loess_cv <- function(y, x, spans, k = 5, robust = TRUE) {
    n <- length(y)
    k <- min(k, max(2, floor(n/2)))
    idx <- seq_len(n)
    folds <- split(idx, cut(idx, breaks = k, labels = FALSE))
    mse <- sapply(spans, function(sp) {
      errs <- numeric(0)
      for (fi in seq_along(folds)) {
        val <- folds[[fi]]; train <- setdiff(idx, val)
        fit <- try(loess(y ~ x, span = sp, degree = 2,
                         family = if (robust) "symmetric" else "gaussian",
                         subset = train), silent = TRUE)
        if (inherits(fit, "try-error")) next
        pred <- try(predict(fit, newdata = data.frame(x = x[val])), silent = TRUE)
        if (inherits(pred, "try-error")) next
        errs <- c(errs, mean((y[val] - pred)^2, na.rm = TRUE))
      }
      mean(errs, na.rm = TRUE)
    })
    spans[which.min(mse)]
  }

  # ---- Trend (base, without Fourier) ----
  df$FourierComponent <- 0
  trendlabel <- NULL
  trendobj <- NULL

  # Decide LOESS span if needed
  span_opt <- NULL
  if (trendmethod == "loess") {
    if (loess_span_mode == "fixed") {
      if (is.null(loess_span_fixed)) stop("loess_span_fixed must be provided when loess_span_mode='fixed'.")
      span_opt <- loess_span_fixed
      message(sprintf("Using fixed LOESS span = %.3f", span_opt))
    } else {
      have_fancova <- requireNamespace("fANCOVA", quietly = TRUE)
      if (have_fancova && loess_span_mode %in% c("auto_aicc", "auto_gcv")) {
        crit <- if (loess_span_mode == "auto_aicc") "aicc" else "gcv"
        as_fit <- try(
          fANCOVA::loess.as(
            x = df$timenum, y = df$Value,
            criterion = crit,
            family = if (robust) "symmetric" else "gaussian",
            plot = FALSE
          ),
          silent = TRUE
        )
        if (!inherits(as_fit, "try-error") && !is.null(as_fit$pars$span)) {
          span_opt <- as_fit$pars$span
          message(sprintf("loess.as proposed span (criterion=%s) = %.3f", crit, span_opt))
        }
      }
      if (is.null(span_opt) || loess_span_mode == "cv") {
        span_opt <- blocked_loess_cv(
          y = df$Value, x = df$timenum,
          spans = loess_span_grid, k = loess_cv_k, robust = robust
        )
        message(sprintf("CV-selected LOESS span = %.3f", span_opt))
      }
    }
    loessfit <- loess(Value ~ timenum, data = df,
                      span = span_opt, degree = 2,
                      family = if (robust) "symmetric" else "gaussian")
    base_trend <- as.numeric(predict(loessfit, newdata = data.frame(timenum = timenum)))
    trendobj <- loessfit
    trendlabel <- sprintf("LOESS (span=%.2f%s)", span_opt, if (robust) ", robust" else "")
  } else {
    # GAM or GAMM
    family_for_fit <- if (use_gamm) gaussian() else if (robust) mgcv::scat() else gaussian()
    if (use_gamm) {
      # Construct random and correlation structures
      rand_effect <- if (!is.null(random_effect)) random_effect else {
        if (!is.null(group_var)) setNames(list(~1), group_var) else NULL
      }
      correlation <- NULL
      if (cor_struct != "none") {
        if (!requireNamespace("nlme", quietly = TRUE)) stop("nlme is required for cor_struct.")
        if (cor_struct == "ar1") {
          form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
          correlation <- nlme::corAR1(form = form)
        } else if (cor_struct == "arma") {
          form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
          correlation <- nlme::corARMA(p = arma_p, q = arma_q, form = form)
        }
      }
      gammfit <- mgcv::gamm(
        Value ~ s(timenum, bs = "cs"),
        data = df,
        family = family_for_fit,
        method = "REML",
        random = rand_effect,
        correlation = correlation
      )
      # In-sample predictions including random effects if available
      base_trend <- try(as.numeric(stats::predict(gammfit$lme,
                                                  newdata = df,
                                                  level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                        silent = TRUE)
      if (inherits(base_trend, "try-error") || any(!is.finite(base_trend))) {
        # Fallback to fixed-effect smooth only
        base_trend <- as.numeric(predict(gammfit$gam, newdata = df))
      }
      trendobj <- gammfit
      edf <- tryCatch(sum(summary(gammfit$gam)$s.table[, "edf"]), error = function(e) NA)
      trendlabel <- paste0("GAMM (REML",
                           if (!is.na(edf)) paste0(", edf approx", round(edf, 2)) else "",
                           if (!is.null(group_var) || !is.null(random_effect)) ", random effects" else "",
                           if (cor_struct != "none") paste0(", cor=", cor_struct) else "",
                           ")")
    } else {
      gamfit <- gam(Value ~ s(timenum, bs = "cs"), data = df, method = "REML",
                    family = family_for_fit)
      base_trend <- as.numeric(predict(gamfit, newdata = df))
      trendobj <- gamfit
      edf <- tryCatch(sum(summary(gamfit)$s.table[, "edf"]), error = function(e) NA)
      trendlabel <- paste0("GAM (REML",
                           if (!is.na(edf)) paste0(", edf approx", round(edf, 2)) else "",
                           if (robust && !use_gamm) ", robust=scat" else "",
                           ")")
    }
  }

  # ---- Automatic Fourier K selection via AICc/BIC ----
  selectedK <- 0
  fourfit <- NULL
  if (usefourier) {
    resid_base <- df$Value - base_trend
    crit_vals <- rep(Inf, fourierK_max + 1)   # index 1 corresponds to K=0
    fit_list <- vector("list", length = fourierK_max + 1)
    # K = 0 case (no Fourier)
    crit_vals[1] <- {
      rss <- sum((resid_base - mean(resid_base))^2, na.rm = TRUE)
      n <- sum(is.finite(resid_base))
      kparams <- 1
      if (fourier_selection_criterion == "BIC") {
        n*log(rss/n) + kparams*log(n)
      } else {
        aic <- n*log(rss/n) + 2*kparams
        if (n - kparams - 1 > 0) aic + (2*kparams*(kparams+1))/(n - kparams - 1) else aic
      }
    }
    if (!is.null(full_fourier_design)) {
      for (K in 1:fourierK_max) {
        terms <- paste(
          unlist(lapply(seq_len(K), function(k) c(paste0("sin", k), paste0("cos", k)))),
          collapse = " + "
        )
        fml <- as.formula(paste0("resid_base ~ ", terms))
        fitK <- try(lm(fml, data = cbind(df, full_fourier_design)), silent = TRUE)
        if (inherits(fitK, "try-error")) next
        fit_list[[K + 1]] <- fitK
        n <- nrow(df)
        kparams <- length(coef(fitK))
        rss <- sum(residuals(fitK)^2)
        if (fourier_selection_criterion == "BIC") {
          crit_vals[K + 1] <- n*log(rss/n) + kparams*log(n)
        } else {
          aic <- n*log(rss/n) + 2*kparams
          crit_vals[K + 1] <- if (n - kparams - 1 > 0) aic + (2*kparams*(kparams+1))/(n - kparams - 1) else aic
        }
      }
      selected_ix <- which.min(crit_vals)
      selectedK <- selected_ix - 1
      if (selectedK > 0) {
        fourfit <- fit_list[[selected_ix]]
        df$FourierComponent <- as.numeric(predict(fourfit, newdata = cbind(df, full_fourier_design)))
        message(sprintf("Fourier selection: criterion=%s -> K=%d (from 0..%d)",
                        fourier_selection_criterion, selectedK, fourierK_max))
      } else {
        df$FourierComponent <- 0
        message(sprintf("Fourier selection: criterion=%s -> K=0 (no harmonics)", fourier_selection_criterion))
      }
    } else {
      df$FourierComponent <- 0
      selectedK <- 0
      message("Fourier design not available; skipping Fourier component.")
    }
  } else {
    selectedK <- fourierK
    df$FourierComponent <- 0
  }

  # ---- Final Trend (base + Fourier) ----
  df$Trend <- base_trend + df$FourierComponent
  if (usefourier) {
    trendlabel <- sprintf("%s + Fourier (K=%d; period approx %.1f days)", trendlabel, selectedK, perioddays)
  }

  # ---- Residuals & diagnostics ----
  residuals <- df$Value - df$Trend
  if (is.null(lagmax)) lagmax <- max(5L, floor(N/10))
  normalitytest <- nortest::sf.test(residuals)
  acfvalues <- acf(residuals, plot = FALSE, lag.max = lagmax)
  acfsig <- any(abs(acfvalues$acf[-1]) > 1.96 / sqrt(length(residuals)))
  ljungboxtest <- Box.test(residuals, lag = lagmax, type = "Ljung-Box")

  # ---- Confidence intervals ----
  if (missing(cimethod) || is.null(cimethod)) {
    cimethod <- if (ljungboxtest$p.value < 0.05) "bootstrapmbb" else "bootstrapiid"
  }
  cilower <- ciupper <- rep(NA_real_, N)

  # helper for bootstrap refits (lock selectedK)
  refitpredict <- function(yvec) {
    if (trendmethod == "loess") {
      fit1 <- loess(yvec ~ timenum, span = if (exists("span_opt")) span_opt else 0.3,
                    degree = 2, family = if (robust) "symmetric" else "gaussian")
      tr1 <- as.numeric(predict(fit1, newdata = data.frame(timenum = timenum)))
    } else {
      if (use_gamm) {
        rand_effect <- if (!is.null(random_effect)) random_effect else {
          if (!is.null(group_var)) setNames(list(~1), group_var) else NULL
        }
        correlation <- NULL
        if (cor_struct != "none") {
          if (cor_struct == "ar1") {
            form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
            correlation <- nlme::corAR1(form = form)
          } else if (cor_struct == "arma") {
            form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
            correlation <- nlme::corARMA(p = arma_p, q = arma_q, form = form)
          }
        }
        fit <- mgcv::gamm(
          Value ~ s(timenum, bs = "cs"),
          data = transform(df, Value = yvec),
          family = gaussian(), method = "REML",
          random = rand_effect, correlation = correlation
        )
        tr1 <- try(as.numeric(stats::predict(fit$lme,
                                             newdata = df,
                                             level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                   silent = TRUE)
        if (inherits(tr1, "try-error") || any(!is.finite(tr1))) {
          tr1 <- as.numeric(predict(fit$gam, newdata = df))
        }
      } else {
        fit <- gam(Value ~ s(timenum, bs = "cs"),
                   data = transform(df, Value = yvec),
                   method = "REML",
                   family = if (robust) mgcv::scat() else gaussian())
        tr1 <- as.numeric(predict(fit, newdata = df))
      }
    }
    if (usefourier && selectedK > 0 && !is.null(full_fourier_design)) {
      resid1 <- yvec - tr1
      terms <- paste(unlist(lapply(seq_len(selectedK), function(k) c(paste0("sin", k), paste0("cos", k)))), collapse = " + ")
      fml <- as.formula(paste0("resid1 ~ ", terms))
      fit2 <- lm(fml, data = cbind(df, full_fourier_design))
      four <- as.numeric(predict(fit2, newdata = cbind(df, full_fourier_design)))
      tr1 + four
    } else tr1
  }

  if (cimethod == "model") {
    # Model-based intervals implemented only for base LOESS/GAM without Fourier
    if (trendmethod == "loess" && (!usefourier || selectedK == 0)) {
      pred <- predict(trendobj, newdata = data.frame(timenum = timenum), se = TRUE)
      ciupper <- pred$fit + 1.96 * pred$se.fit
      cilower <- pred$fit - 1.96 * pred$se.fit
    } else if (trendmethod == "gam" && !use_gamm && (!usefourier || selectedK == 0)) {
      pred <- predict(trendobj, newdata = df, se.fit = TRUE)
      ciupper <- pred$fit + 1.96 * pred$se.fit
      cilower <- pred$fit - 1.96 * pred$se.fit
    } else {
      warning("Model-based CI not implemented for this configuration; switching to bootstrapiid.")
      cimethod <- "bootstrapiid"
    }
  }

  if (cimethod %in% c("bootstrapiid", "bootstrapmbb")) {
    if (cimethod == "bootstrapmbb") {
      if (!is.null(blocksize)) {
        message(sprintf("Using user-specified MBB blocksize = %d", as.integer(blocksize)))
      } else if (blocklength_mode == "fixed") {
        if (is.null(blocklength_fixed)) stop("blocklength_fixed must be provided when blocklength_mode='fixed'.")
        blocksize <- as.integer(blocklength_fixed)
        message(sprintf("Using fixed MBB blocksize = %d", blocksize))
      } else if (blocklength_mode == "auto_pwsd") {
        have_blocklength <- requireNamespace("blocklength", quietly = TRUE)
        if (have_blocklength) {
          pw <- try(blocklength::pwsd(series = residuals, plots = FALSE), silent = TRUE)
          if (!inherits(pw, "try-error") && !is.null(pw$block_length)) {
            blocksize <- as.integer(round(pw$block_length))
            message(sprintf("PWSD-selected MBB blocksize = %d", blocksize))
          }
        }
        if (is.null(blocksize)) {
          n <- length(residuals)
          crit <- 1.96 / sqrt(n)
          acfobj <- acf(residuals, plot = FALSE, lag.max = max(5L, floor(n/10)))
          siglags <- which(abs(acfobj$acf[-1]) > crit)
          corr_len <- if (length(siglags)) min(siglags) else 1L
          blocksize <- as.integer(max(5L, floor(max(n^(1/3), 1.5 * corr_len))))
          message(sprintf("Heuristic MBB blocksize = %d (fallback)", blocksize))
        }
      } else {
        n <- length(residuals)
        crit <- 1.96 / sqrt(n)
        acfobj <- acf(residuals, plot = FALSE, lag.max = max(5L, floor(n/10)))
        siglags <- which(abs(acfobj$acf[-1]) > crit)
        corr_len <- if (length(siglags)) min(siglags) else 1L
        blocksize <- as.integer(max(5L, floor(max(n^(1/3), 1.5 * corr_len))))
        message(sprintf("Heuristic MBB blocksize = %d", blocksize))
      }
    }

    mbbresample <- function(residvec, B, blocksize) {
      n  <- length(residvec)
      bs <- as.integer(blocksize)
      out <- matrix(NA_real_, nrow = n, ncol = B)
      starts <- 1:(n - bs + 1)
      for (b in seq_len(B)) {
        blocks <- integer(0)
        while (length(blocks) < n) {
          s <- sample(starts, 1)
          blocks <- c(blocks, residvec[s:(s + bs - 1)])
        }
        out[, b] <- blocks[seq_len(n)]
      }
      out
    }

    boottrends <- matrix(NA_real_, nrow = N, ncol = nboot)
    if (cimethod == "bootstrapiid") {
      for (b in seq_len(nboot)) {
        epsb <- sample(residuals, replace = TRUE)
        yb   <- df$Trend + epsb
        boottrends[, b] <- refitpredict(yb)
      }
    } else {
      mbbmat <- mbbresample(residuals, nboot, blocksize)
      for (b in seq_len(nboot)) {
        yb <- df$Trend + mbbmat[, b]
        boottrends[, b] <- refitpredict(yb)
      }
    }
    cilower <- apply(boottrends, 1, quantile, probs = 0.025, na.rm = TRUE)
    ciupper <- apply(boottrends, 1, quantile, probs = 0.975, na.rm = TRUE)
    message(sprintf("Bootstrap CI (%s): %d resamples%s",
                    cimethod, nboot,
                    if (cimethod == "bootstrapmbb") paste0("; blocksize=", blocksize) else ""))
  }

  # ---- Outliers ----
  df$Outlier <- df$Trend < cilower | df$Trend > ciupper

  # ---- Stationarity tests, change-points, spectra ----
  adftest  <- tseries::adf.test(signal, alternative = "stationary")
  kpsstest <- tseries::kpss.test(signal, null = "Level")
  adfinterpret  <- if (adftest$p.value < 0.05) "ADF: Reject H0 -> likely stationary" else "ADF: Fail to reject H0 -> likely non stationary"
  kpssinterpret <- if (kpsstest$p.value < 0.05) "KPSS: Reject H0 -> likely non stationary" else "KPSS: Fail to reject H0 -> likely stationary"

  cpobj    <- cpt.meanvar(df$Trend, method = "PELT", penalty = "BIC")
  cppoints <- cpts(cpobj)
  cpdates  <- if (length(cppoints)) df$Date[cppoints] else as.Date(character())

  # ---- Spectra for plotting ----
  spectrumdf <- significantspectrum <- data.frame()
  significantperiods <- numeric(0)
  dominantperioddays <- dom$period_days
  dominantfrequency  <- dom$freq
  if (!irregular) {
    tsdata <- ts(signal, frequency = seasonalfrequency)
    decomp <- stl(tsdata, s.window = "periodic", robust = stlrobust)
    df$Seasonal <- as.numeric(decomp$time.series[, "seasonal"])
    spec  <- spec.pgram(signal, spans = specspans, log = "no", taper = 0.1, detrend = TRUE)
    freq  <- spec$freq; power <- spec$spec
    spectrumdf <- data.frame(Frequency = freq, Power = power)
    k      <- spec$df; alpha <- 0.01
    upper  <- power * qchisq(1 - alpha/2, 2 * k) / (2 * k)
    sigix  <- which(power > upper)
    significantspectrum <- spectrumdf[sigix, , drop = FALSE]
    significantperiods <- round(1 / significantspectrum$Frequency * avgintervaldays, 2)
  } else {
    tnum <- as.numeric(df$Date - min(df$Date))
    lombres  <- lsp(signal, tnum, type = "frequency", plot = FALSE)
    lombdf   <- data.frame(Frequency = lombres$scanned, Power = lombres$power)
    spectrumdf <- lombdf
    signiflevel <- lombres$sig.level
    sigix <- which(lombres$power > signiflevel)
    significantspectrum <- lombdf[sigix, , drop = FALSE]
    significantperiods  <- round(1 / significantspectrum$Frequency, 2)
  }

  # ---- Plots ----
  ptrend <- ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Value), color = "steelblue", linewidth = 0.8) +
    geom_line(aes(y = Trend), color = "darkred", linewidth = 0.9) +
    geom_ribbon(aes(ymin = cilower, ymax = ciupper), fill = "grey80", alpha = 0.35) +
    geom_point(data = subset(df, Outlier), aes(y = Trend), color = "orange", size = 2) +
    { if (length(cpdates)) geom_vline(xintercept = cpdates, linetype = "dotted", color = "orange") } +
    labs(title = sprintf("Trend (%s) with 95%% CI (%s)", trendlabel, cimethod),
         y = "Signal", x = "Date") +
    theme_minimal()

  pperiod <- ggplot(spectrumdf, aes(x = Frequency, y = Power)) +
    geom_line(color = "black", linewidth = 0.6) +
    geom_point(data = significantspectrum, aes(x = Frequency, y = Power), color = "red", size = 1.6) +
    annotate("text",
             x = dominantfrequency,
             y = max(spectrumdf$Power, na.rm = TRUE),
             label = paste0("Dominant period approx ", round(dominantperioddays, 1), " days"),
             hjust = -0.1, vjust = -0.5, color = "red") +
    labs(title = if (!irregular) "Smoothed periodogram (FFT)" else "Lomb Scargle periodogram",
         x = "Frequency (cycles per time unit)", y = "Power") +
    theme_minimal()

  acfdf <- data.frame(Lag = acfvalues$lag[-1], ACF = acfvalues$acf[-1])
  pacf <- ggplot(acfdf, aes(x = Lag, y = ACF)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = c(-1.96/sqrt(length(residuals)), 1.96/sqrt(length(residuals))),
               linetype = "dashed", color = "red") +
    labs(title = sprintf("Residual ACF (lag.max=%d)", lagmax),
         x = "Lag", y = "ACF") +
    theme_minimal()

  combinedplot <- arrangeGrob(ptrend, pperiod, pacf, ncol = 1)
  grid.arrange(combinedplot)

  # ---- Helper plot ----
  fourierplot <- NULL
  if (usefourier && selectedK > 0) {
    plotdf <- df |>
      mutate(Residual = Value - base_trend) |>
      dplyr::select(Date, FourierComponent, Residual)
    fourierplot <- ggplot(plotdf, aes(x = Date)) +
      geom_line(aes(y = Residual, color = "Residual (Value base trend)"), linewidth = 0.7, alpha = 0.65) +
      geom_line(aes(y = FourierComponent, color = "Fourier component"), linewidth = 1.0) +
      scale_color_manual(values = c("Residual (Value base trend)" = "grey40",
                                    "Fourier component" = "#7B1FA2")) +
      labs(title = "Helper plot: Fourier (seasonal) component vs Residual",
           y = "Contribution (signal units)", x = "Date", color = NULL) +
      theme_minimal()
  }

  # ---- Interpretations (recommend GAMM when cycles + autocorr) ----
  acfcomment <- if (ljungboxtest$p.value > 0.05 && acfsig) {
    "ACF has some significant lags, but Ljung Box does not reject -> dependence weak/non-systematic."
  } else if (ljungboxtest$p.value < 0.05) {
    "Ljung Box indicates significant serial correlation in residuals."
  } else {
    "Residuals show no meaningful autocorrelation."
  }
  cyclecomment <- if (ljungboxtest$p.value < 0.05 && length(significantperiods) == 0) {
    "Autocorrelation without significant cycles -> short run dependence; consider ARIMA or GAMM with correlation."
  } else if (ljungboxtest$p.value < 0.05 && length(significantperiods) > 0) {
    "Autocorrelation and significant cycles -> consider SARIMA or GAMM with Fourier terms."
  } else {
    "No strong autocorrelation/cycles -> trend model likely sufficient."
  }
  interpretationsummary <- paste0(
    "Trend: ", trendlabel, ". ",
    if (normalitytest$p.value > 0.05) "Residuals approximately normal. " else "Residuals deviate from normality. ",
    acfcomment, " ",
    if (length(significantperiods) > 0) {
      paste0("Detected cycles (periods approx ", paste(significantperiods, collapse = ", "), " days). ")
    } else {
      "No significant cycles detected. "
    },
    if (!is.na(dominantperioddays)) paste0("Dominant cycle approx ", round(dominantperioddays, 2), " days. ") else "",
    adfinterpret, ". ",
    kpssinterpret, ". ",
    cyclecomment
  )

  trendpresence <- if (adftest$p.value >= 0.05 && kpsstest$p.value < 0.05) "Non stationary (trend present)"
                   else if (adftest$p.value < 0.05 && kpsstest$p.value >= 0.05) "Stationary (trend limited)"
                   else "Ambiguous (check plot/tests)"

  # ---- Rolling-origin forecast evaluation ----
  forecast <- NULL
  if (forecast_holdout_h > 0) {
    H <- forecast_holdout_h
    yhat <- rep(NA_real_, H)
    yobs <- tail(df$Value, H)
    hold_idx <- (N - H + 1):N
    train_initial_end <- N - H
    if (forecast_origin_mode == "expanding") {
      for (j in seq_along(hold_idx)) {
        tr_idx <- 1:(train_initial_end + (j - 1))
        if (trendmethod == "loess") {
          fit1 <- loess(Value ~ timenum, data = df[tr_idx, ],
                        span = span_opt, degree = 2,
                        family = if (robust) "symmetric" else "gaussian")
          tr1 <- as.numeric(predict(fit1, newdata = df[hold_idx[j], "timenum", drop = FALSE]))
        } else {
          if (use_gamm) {
            rand_effect <- if (!is.null(random_effect)) random_effect else {
              if (!is.null(group_var)) setNames(list(~1), group_var) else NULL
            }
            correlation <- NULL
            if (cor_struct != "none") {
              if (cor_struct == "ar1") {
                form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
                correlation <- nlme::corAR1(form = form)
              } else if (cor_struct == "arma") {
                form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
                correlation <- nlme::corARMA(p = arma_p, q = arma_q, form = form)
              }
            }
            gfit <- mgcv::gamm(Value ~ s(timenum, bs = "cs"),
                               data = df[tr_idx, ],
                               family = gaussian(), method = "REML",
                               random = rand_effect, correlation = correlation)
            tr1 <- try(as.numeric(stats::predict(gfit$lme,
                                                 newdata = df[hold_idx[j], , drop = FALSE],
                                                 level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                       silent = TRUE)
            if (inherits(tr1, "try-error") || any(!is.finite(tr1))) {
              tr1 <- as.numeric(predict(gfit$gam, newdata = df[hold_idx[j], , drop = FALSE]))
            }
          } else {
            gfit <- gam(Value ~ s(timenum, bs = "cs"),
                        data = df[tr_idx, ], method = "REML",
                        family = if (robust) mgcv::scat() else gaussian())
            tr1 <- as.numeric(predict(gfit, newdata = df[hold_idx[j], , drop = FALSE]))
          }
        }
        fcomp <- 0
        if (usefourier && selectedK > 0 && !is.null(full_fourier_design)) {
          resid_train <- df$Value[tr_idx] - (if (trendmethod == "loess")
            as.numeric(predict(fit1, newdata = df[tr_idx, "timenum", drop = FALSE])) else
            if (use_gamm)
              tryCatch(as.numeric(stats::predict(gfit$lme,
                                                 newdata = df[tr_idx, , drop = FALSE],
                                                 level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                       error = function(e) as.numeric(predict(gfit$gam, newdata = df[tr_idx, , drop = FALSE])))
            else
              as.numeric(predict(gfit, newdata = df[tr_idx, , drop = FALSE])))
          if (forecast_lock_K) {
            terms <- paste(unlist(lapply(seq_len(selectedK), function(k) c(paste0("sin", k), paste0("cos", k)))), collapse = " + ")
            fml <- as.formula(paste0("resid_train ~ ", terms))
            fit2 <- lm(fml, data = cbind(df[tr_idx, ], full_fourier_design[tr_idx, , drop = FALSE]))
            fcomp <- as.numeric(predict(fit2, newdata = cbind(df[hold_idx[j], , drop = FALSE],
                                                              full_fourier_design[hold_idx[j], , drop = FALSE])))
          } else {
            bestK <- 0; bestCrit <- Inf; bestFit <- NULL
            for (K in 1:fourierK_max) {
              terms <- paste(unlist(lapply(seq_len(K), function(k) c(paste0("sin", k), paste0("cos", k)))), collapse = " + ")
              fml <- as.formula(paste0("resid_train ~ ", terms))
              fitK <- try(lm(fml, data = cbind(df[tr_idx, ], full_fourier_design[tr_idx, , drop = FALSE])), silent = TRUE)
              if (inherits(fitK, "try-error")) next
              n <- length(tr_idx); kparams <- length(coef(fitK)); rss <- sum(residuals(fitK)^2)
              crit <- if (fourier_selection_criterion == "BIC") n*log(rss/n) + kparams*log(n) else {
                aic <- n*log(rss/n) + 2*kparams
                if (n - kparams - 1 > 0) aic + (2*kparams*(kparams+1))/(n - kparams - 1) else aic
              }
              if (crit < bestCrit) { bestCrit <- crit; bestK <- K; bestFit <- fitK }
            }
            if (!is.null(bestFit) && bestK > 0) {
              fcomp <- as.numeric(predict(bestFit, newdata = cbind(df[hold_idx[j], , drop = FALSE],
                                                                   full_fourier_design[hold_idx[j], , drop = FALSE])))
            } else fcomp <- 0
          }
        }
        yhat[j] <- tr1 + fcomp
      }
    } else {
      if (is.null(train_window)) train_window <- train_initial_end
      for (j in seq_along(hold_idx)) {
        endj <- hold_idx[j] - 1
        startj <- max(1L, endj - train_window + 1L)
        tr_idx <- startj:endj
        if (trendmethod == "loess") {
          fit1 <- loess(Value ~ timenum, data = df[tr_idx, ],
                        span = span_opt, degree = 2,
                        family = if (robust) "symmetric" else "gaussian")
          tr1 <- as.numeric(predict(fit1, newdata = df[hold_idx[j], "timenum", drop = FALSE]))
        } else {
          if (use_gamm) {
            rand_effect <- if (!is.null(random_effect)) random_effect else {
              if (!is.null(group_var)) setNames(list(~1), group_var) else NULL
            }
            correlation <- NULL
            if (cor_struct != "none") {
              if (cor_struct == "ar1") {
                form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
                correlation <- nlme::corAR1(form = form)
              } else if (cor_struct == "arma") {
                form <- if (!is.null(group_var)) as.formula(paste0("~ timenum | ", group_var)) else ~ timenum
                correlation <- nlme::corARMA(p = arma_p, q = arma_q, form = form)
              }
            }
            gfit <- mgcv::gamm(Value ~ s(timenum, bs = "cs"),
                               data = df[tr_idx, ],
                               family = gaussian(), method = "REML",
                               random = rand_effect, correlation = correlation)
            tr1 <- try(as.numeric(stats::predict(gfit$lme,
                                                 newdata = df[hold_idx[j], , drop = FALSE],
                                                 level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                       silent = TRUE)
            if (inherits(tr1, "try-error") || any(!is.finite(tr1))) {
              tr1 <- as.numeric(predict(gfit$gam, newdata = df[hold_idx[j], , drop = FALSE]))
            }
          } else {
            gfit <- gam(Value ~ s(timenum, bs = "cs"),
                        data = df[tr_idx, ], method = "REML",
                        family = if (robust) mgcv::scat() else gaussian())
            tr1 <- as.numeric(predict(gfit, newdata = df[hold_idx[j], , drop = FALSE]))
          }
        }
        fcomp <- 0
        if (usefourier && selectedK > 0 && !is.null(full_fourier_design)) {
          resid_train <- df$Value[tr_idx] - (if (trendmethod == "loess")
            as.numeric(predict(fit1, newdata = df[tr_idx, "timenum", drop = FALSE])) else
            if (use_gamm)
              tryCatch(as.numeric(stats::predict(gfit$lme,
                                                 newdata = df[tr_idx, , drop = FALSE],
                                                 level = if (!is.null(group_var) || !is.null(random_effect)) 1 else 0)),
                       error = function(e) as.numeric(predict(gfit$gam, newdata = df[tr_idx, , drop = FALSE])))
            else
              as.numeric(predict(gfit, newdata = df[tr_idx, , drop = FALSE])))
          terms <- paste(unlist(lapply(seq_len(selectedK), function(k) c(paste0("sin", k), paste0("cos", k)))), collapse = " + ")
          fml <- as.formula(paste0("resid_train ~ ", terms))
          fit2 <- lm(fml, data = cbind(df[tr_idx, ], full_fourier_design[tr_idx, , drop = FALSE]))
          fcomp <- as.numeric(predict(fit2, newdata = cbind(df[hold_idx[j], , drop = FALSE],
                                                            full_fourier_design[hold_idx[j], , drop = FALSE])))
        }
        yhat[j] <- tr1 + fcomp
      }
    }
    err <- yhat - yobs
    rmse <- sqrt(mean(err^2, na.rm = TRUE))
    mae  <- mean(abs(err), na.rm = TRUE)
    mape <- mean(abs(err / yobs)[is.finite(yobs) & yobs != 0], na.rm = TRUE)
    smape <- mean(2 * abs(err) / (abs(yobs) + abs(yhat)), na.rm = TRUE)
    forecast <- list(
      holdout_index = hold_idx,
      y_obs = yobs,
      y_hat = yhat,
      errors = err,
      metrics = list(RMSE = rmse, MAE = mae, MAPE = mape, sMAPE = smape),
      mode = forecast_origin_mode,
      H = H
    )
  }

  # ---- Summary table ----
  summarytable <- data.frame(
    Metric = c("Sampling regularity",
               "Trend method",
               "Fourier harmonics (K)",
               "Fourier selection criterion",
               "CI method",
               "Trend presence (ADF+KPSS)",
               "Dominant cycle (days)",
               "Detected cycle periods (days)",
               "Change-points (count)",
               "Change-point dates",
               "Residual autocorrelation",
               "Residual normality (Shapiro Francia)",
               "Forecast hold-out (H)",
               "Forecast mode",
               "Forecast RMSE",
               "Forecast MAE",
               "Forecast MAPE",
               "Forecast sMAPE"),
    Value = c(ifelse(irregular, "Irregular", "Regular"),
              trendlabel,
              if (usefourier) selectedK else "None",
              if (usefourier) fourier_selection_criterion else "N/A",
              cimethod,
              trendpresence,
              ifelse(is.na(dominantperioddays), "N/A", round(dominantperioddays, 2)),
              ifelse(length(significantperiods), paste(significantperiods, collapse = ", "), "None"),
              length(cpdates),
              ifelse(length(cpdates), paste(format(cpdates), collapse = "; "), " "),
              if (ljungboxtest$p.value < 0.05) "Significant" else if (acfsig) "Weak" else "None",
              if (normalitytest$p.value > 0.05) "Pass" else "Fail",
              if (!is.null(forecast)) forecast$H else 0,
              if (!is.null(forecast)) forecast$mode else "N/A",
              if (!is.null(forecast)) round(forecast$metrics$RMSE, 4) else "N/A",
              if (!is.null(forecast)) round(forecast$metrics$MAE, 4)  else "N/A",
              if (!is.null(forecast)) round(forecast$metrics$MAPE, 4) else "N/A",
              if (!is.null(forecast)) round(forecast$metrics$sMAPE, 4) else "N/A"
    ),
    stringsAsFactors = FALSE
  )

  # ---- Export plot & DOCX ----
  if (exportplot) {
    ggsave("analysisplot.png", combinedplot, width = 8, height = 12, dpi = 300)
    if (!is.null(fourierplot)) {
      ggsave("fourierplot.png", fourierplot, width = 8, height = 4, dpi = 300)
    }
  }
  if (exportdocx) {
    testresults <- data.frame(
      Test    = c("Shapiro Francia", "Ljung Box", "ADF", "KPSS"),
      PValue  = signif(c(normalitytest$p.value,
                         ljungboxtest$p.value,
                         adftest$p.value,
                         kpsstest$p.value), 4)
    )
    cycleinfo <- if (length(significantperiods)) {
      data.frame(SignificantPeriodsDays = significantperiods)
    } else data.frame(SignificantPeriodsDays = "None")

    metatable <- data.frame(
      Field = c("Project ID", "Cohort ID", "Assay version", "Analyst", "Run date", "Notes",
                "GAMM used", "Group var", "Correlation"),
      Value = c(project_id %||% "",
                cohort_id %||% " ",
                assay_version %||% " ",
                analyst %||% " ",
                as.character(run_date %||% Sys.Date()),
                notes %||% " ",
                if (use_gamm) "Yes" else "No",
                group_var %||% " ",
                if (use_gamm) cor_struct else " "),
      stringsAsFactors = FALSE
    )

    param_list <- list(
      normalize = normalize,
      trendmethod = trendmethod,
      robust = robust,
      use_gamm = use_gamm,
      group_var = group_var %||% NA,
      cor_struct = if (use_gamm) cor_struct else NA,
      arma_p = if (use_gamm && cor_struct == "arma") arma_p else NA,
      arma_q = if (use_gamm && cor_struct == "arma") arma_q else NA,
      loess_span_mode = loess_span_mode,
      loess_span_fixed = loess_span_fixed %||% NA,
      loess_span_grid = paste(round(loess_span_grid, 2), collapse = ","),
      loess_cv_k = loess_cv_k,
      usefourier = usefourier,
      selectedK = if (usefourier) selectedK else NA,
      fourierK_max = fourierK_max,
      fourier_selection_criterion = if (usefourier) fourier_selection_criterion else NA,
      auto_seasonality = auto_seasonality,
      seasonalfrequency = seasonalfrequency,
      perioddays = round(perioddays, 3),
      specspans = paste(specspans, collapse = ","),
      cimethod = cimethod,
      nboot = nboot,
      blocksize = blocksize %||% NA,
      blocklength_mode = blocklength_mode,
      forecast_holdout_h = forecast_holdout_h,
      forecast_origin_mode = forecast_origin_mode,
      train_window = train_window %||% NA,
      forecast_lock_K = forecast_lock_K
    )
    parameters_table <- data.frame(
      Parameter = names(param_list),
      Value = unlist(param_list),
      stringsAsFactors = FALSE
    )

    doc <- read_docx()
    if (!is.null(logopath) && file.exists(logopath)) {
      doc <- doc %>% body_add_img(src = logopath, width = logowidth, height = logoheight)
    }
    doc <- doc %>%
      body_add_par("Signal Analysis Report", style = "heading 1") %>%
      body_add_par(interpretationsummary, style = "Normal") %>%
      body_add_par("Summary", style = "heading 2") %>%
      body_add_flextable(flextable(summarytable)) %>%
      body_add_par("Key Tests", style = "heading 2") %>%
      body_add_flextable(flextable(testresults)) %>%
      body_add_par("Detected Cycles", style = "heading 2") %>%
      body_add_flextable(flextable(cycleinfo)) %>%
      body_add_par("Metadata", style = "heading 2") %>%
      body_add_flextable(flextable(metatable))

    if (!is.null(forecast)) {
      fct <- data.frame(
        Index = forecast$holdout_index,
        Date = df$Date[forecast$holdout_index],
        Observed = forecast$y_obs,
        Predicted = forecast$y_hat,
        Error = forecast$errors
      )
      doc <- doc %>%
        body_add_par("Forecast Evaluation (Rolling Origin)", style = "heading 2") %>%
        body_add_par(sprintf("Mode: %s | H=%d | RMSE=%.4f | MAE=%.4f | MAPE=%.4f | sMAPE=%.4f",
                             forecast$mode, forecast$H,
                             forecast$metrics$RMSE, forecast$metrics$MAE,
                             forecast$metrics$MAPE, forecast$metrics$sMAPE), style = "Normal") %>%
        body_add_flextable(flextable(fct))
    }

    if (include_parameters_appendix) {
      doc <- doc %>%
        body_add_par("Appendix: Parameters", style = "heading 2") %>%
        body_add_flextable(flextable(parameters_table))
    }

    if (exportplot && file.exists("analysisplot.png")) {
      doc <- doc %>% body_add_par("Figures", style = "heading 2") %>%
        body_add_img(src = "analysisplot.png", width = 6, height = 8)
    }
    if (exportplot && file.exists("fourierplot.png")) {
      doc <- doc %>% body_add_img(src = "fourierplot.png", width = 6, height = 4)
    }
    print(doc, target = outputpath)
    message("Report exported to: ", outputpath)
  }

  # ---- Return ----
  list(
    trendplot = ptrend,
    periodogramplot = pperiod,
    acfplot = pacf,
    fourierplot = fourierplot,
    spectrum = spectrumdf,
    significantspectrum = significantspectrum,
    significantperiodsdays = significantperiods,
    trendmethod = trendmethod,
    use_gamm = use_gamm,
    cimethod = cimethod,
    irregulardates = irregular,
    trendci95 = data.frame(lower = cilower, upper = ciupper),
    residuals = residuals,
    normalitytest = normalitytest,
    ljungboxtest = ljungboxtest,
    acfvalues = acfvalues,
    ADF = adftest,
    KPSS = kpsstest,
    interpretation = interpretationsummary,
    summarytable = summarytable,
    changepointsindex = cppoints,
    changepointsdates = cpdates,
    dominantperioddays = dominantperioddays,
    dominantfrequency = dominantfrequency,
    fouriercomponent = if (usefourier) df$FourierComponent else NULL,
    auto_seasonality = auto_seasonality,
    selectedK = if (usefourier) selectedK else NULL,
    fourier_selection_criterion = if (usefourier) fourier_selection_criterion else NULL,
    forecast = forecast
  )
}
#End of code
