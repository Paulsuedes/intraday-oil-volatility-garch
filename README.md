# Intraday Oil Volatility Forecasting with (G)ARCH (R)

Forecasting crude oil volatility using GARCH-family models in R with out-of-sample evaluation.
The project compares model classes (sGARCH / EGARCH / GJR-GARCH), distributions (Normal vs Student-t),
different volatility proxies (daily r² / Parkinson / realized variance from 5-minute returns),
and performance in a crisis subsample.

## Key components
- **Data & preprocessing:** daily WTI prices and 5-minute intraday prices → log returns
- **Volatility proxies:** daily r², Parkinson (high/low), realized variance (5-min)
- **Models:** sGARCH(1,1), EGARCH(1,1), GJR-GARCH(1,1)
- **Distributions:** Normal vs Student-t
- **Forecasting:** rolling and expanding windows, 1-step ahead variance forecasts
- **Evaluation:** MSE, QLIKE, Diebold–Mariano tests (QLIKE)
- **Stress period:** crisis subsample (|return| > 3%) + simple 1% VaR

## Repository structure (planned)
- `scripts/` preprocessing, modeling, forecasting, evaluation
- `figures/` key plots
- `data/` not included (licensing)

## Note on data
The original dataset is not included due to licensing.
The code is written to work with any intraday price series in the documented format.

## Tech stack
R / RStudio: `dplyr`, `lubridate`, `rugarch`, `forecast`, `ggplot2`, `tseries`, `moments`
