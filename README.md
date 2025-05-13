# npq: Interval score optimised pairs of nonparametric regression quantiles

R-Codefiles for "Distribution-free prediction intervals from interval score optimised pairs of nonparametric regression quantiles" from Harry Haupt, Joachim Schnurbus, and Ida Bauer.

| Codefile | Description |
|------------------------------------|------------------------------------|
| `NPQ_Functions.R` | Functions to perform interval forecasting with mixed kernel quantile regression (is sourced in other files). |
| `NPQ_Control.R` | File for performing the bandwidth training and creating test-forecasts for $IS$ and $IS_c$ and for measuring the performance (is sourced for example computation). |
| `NPQ_ControlMC.R` | like `NPQ_Control.R`, but tailored to Monte Carlo simulation (is sourced in Monte Carlo simulation). |
| `NPQ_ComputeISS.R` | Function to perform $ISS$ estimation (is sourced in example computation). |
| `NPQ_ComputeISSforMC.R` | Function to perform $ISS$ estimation (is sourced in Monte Carlo simulation). |
| `NPQ_ComputeExamples_90.R` | Function for computing industrial production example (90% coverage). |
| `NPQ_ComputeExamples_80.R` | Function for computing industrial production example (80% coverage). |
| `NPQ_ComputeMC_CC1.R` | Function for computing Monte Carlo simulation. |
