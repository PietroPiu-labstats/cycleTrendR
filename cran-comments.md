## Test environments
* Local Windows 11, R 4.5.2 (ucrt)
* GitHub Actions (Windows, macOS, Ubuntu) — R release, devel, oldrel
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Additional comments
* This is the first submission of **cycleTrendR** to CRAN.
* The package provides adaptive trend estimation, cycle detection, spectral analysis, bootstrap confidence intervals, change‑point detection, and rolling‑origin forecasting for irregularly sampled time series.
* All examples have been reduced in computational cost to comply with CRAN policies.
* The package does not write to the user’s home directory or external files unless explicitly requested by the user (e.g., DOCX export).
* The package passes all checks on multiple platforms, including GitHub Actions and win-builder.
