# Changelog

## cycleTrendR 0.2.0

## cycleTrendR 0.1.0

- Initial release of **cycleTrendR**.
- Added the main function
  [`adaptive_cycle_trend_analysis()`](https://PietroPiu-labstats.github.io/cycleTrendR/reference/adaptive_cycle_trend_analysis.md)
  supporting:
  - LOESS, GAM, and GAMM trend estimation
  - Automatic Fourier harmonic selection (AICc/BIC)
  - Lomb–Scargle periodogram for irregular sampling
  - Bootstrap confidence intervals (IID and MBB)
  - Change-point detection
  - Rolling-origin forecasting
- Added publication-quality ggplot2 visualizations for:
  - Trend + CI
  - Periodogram
  - Residual ACF
  - Diagnostics and summary tables
- Added a comprehensive vignette: *cycleTrendR-overview*.
- Added README with installation instructions and examples.
- Ensured full CRAN compliance: 0 errors, 0 warnings, 0 notes.
