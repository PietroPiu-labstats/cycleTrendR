Time‑series data rarely behave nicely. They’re noisy, irregularly sampled, autocorrelated, and often contain hidden cycles. Traditional tools such as STL or simple GAMs struggle when the data violates assumptions.

cycleTrendR is an R package that unifies:

LOESS, GAM, and GAMM trend estimation

Automatic Fourier selection

Lomb–Scargle periodograms for irregular sampling

Bootstrap confidence intervals

Change‑point detection

Rolling-origin forecasting

All in one function: adaptive_cycle_trend_analysis().

cycleTrendR adapts to the data instead of forcing the data to fit a rigid model.

A quick example

set.seed(1526)

  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 200)
  signal <- sin(2*pi*1:200/30) + rnorm(200, 0, 0.2)

library(cycleTrendR)

res <- adaptive_cycle_trend_analysis(
  signal = signal,
  dates  = dates,
  usefourier = TRUE,
  trendmethod = "gam"
)

plot(res)

The output includes:

Trend line

95% bootstrap CI

Outliers

Change-points

Periodogram

Residual ACF


What makes cycleTrendR different?

Works with irregular time series

Automatically detects cycles

Handles autocorrelation via GAMM

Provides robust uncertainty estimates

Produces publication-ready plots


When should you use it?

Biological assay drift correction

Environmental monitoring

Epidemiological seasonality

Financial cycle detection

Sensor data with autocorrelation
