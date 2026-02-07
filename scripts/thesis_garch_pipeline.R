############################################################
## Packages
############################################################

library(dplyr)
library(lubridate)
library(ggplot2)
library(moments)
library(zoo)
library(rugarch)
library(forecast)

############################################################
## 1. Read data (local files expected in /data)
## Note: Original dataset not included due to licensing.
############################################################

daily_path <- file.path("data", "WTI_1_day_data.csv")
hf_path    <- file.path("data", "WTI_5m_data.csv")

if (!file.exists(daily_path) || !file.exists(hf_path)) {
  stop(
    "Data files not found in /data.\n",
    "Note: The original dataset is not included in this repository due to licensing.\n",
    "See data/README.md for the expected file names and format."
  )
}

WTI.1.day.data <- read.csv(daily_path, header = FALSE, sep = ";", stringsAsFactors = FALSE)
WTI.5.m.data   <- read.csv(hf_path,    header = FALSE, sep = ";", stringsAsFactors = FALSE)

############################################################
## 2. Set column names (OHLCV)
############################################################

colnames(WTI.5.m.data)   <- c("date_chr", "time_chr", "open", "high", "low", "close", "volume")
colnames(WTI.1.day.data) <- c("date_chr", "time_chr", "open", "high", "low", "close", "volume")

############################################################
## 3. Prepare daily data (daily log-returns)
############################################################

wti_daily <- WTI.1.day.data %>%
  mutate(
    date  = dmy(date_chr),
    close = as.numeric(close)
  ) %>%
  filter(!is.na(close), close > 0) %>%
  arrange(date) %>%
  transmute(
    date  = date,
    close = close,
    r_log = 100 * (log(close) - log(lag(close)))   # daily log returns in %
  ) %>%
  filter(!is.na(r_log))    # drop first NA return

############################################################
## 4. Prepare 5-minute data (intraday log-returns)
############################################################

wti_5m = WTI.5.m.data %>%
  mutate(
    date     = dmy(date_chr),
    time     = hm(time_chr),
    datetime = as.POSIXct(date + time, tz = "UTC"),
    close    = as.numeric(close)
  ) %>%
  arrange(datetime) %>%
  group_by(date) %>%
  mutate(
    r_5m = 100 * (log(close) - log(lag(close)))     # 5-minute log returns in %
  ) %>%
  ungroup()


############################################################
## 5. Daily realized variance from 5-minute returns
############################################################

rv_5m_daily = wti_5m %>%
  group_by(date) %>%
  summarise(
    n_5m = sum(!is.na(r_5m)),               # number of 5-min returns per day
    RV_5m = sum(r_5m^2, na.rm = TRUE)       # realized variance (in %^2)
  ) %>%
  filter(n_5m >= 50) %>%                    # keep only days with at least 50 obs
  ungroup()


############################################################
## 6. Merge: daily returns + RV_5m (final data for GARCH)
############################################################

data_for_garch = wti_daily %>%
  inner_join(rv_5m_daily, by = "date") %>%
  arrange(date)

# quick check
head(data_for_garch)
summary(data_for_garch$RV_5m)


############################################################
## 7. Distribution of daily returns
############################################################

# 7.1 Basic statistics
ret = data_for_garch$r_log

summary(ret)
skewness(ret)
kurtosis(ret)   # > 3 => fat tails

# 7.2 Parameters for normal and t distribution
mu     = mean(ret)
sigma  = sd(ret)
k_hat  = kurtosis(ret)

# simple df estimate from kurtosis (only meaningful if k_hat > 3)
df_t = (4 * k_hat - 6) / (k_hat - 3)
df_t


############################################################
## 8. Histogram with normal and t density (with legend)
############################################################

ggplot(data_for_garch, aes(x = r_log)) +
  geom_histogram(aes(y = ..density..),
                 bins = 80,
                 fill = "grey80",
                 color = "grey40") +
  stat_function(
    fun = function(x) dnorm(x, mean = mu, sd = sigma),
    aes(color = "Normal distribution"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  stat_function(
    fun = function(x) dt((x - mu) / sigma, df = df_t) / sigma,
    aes(color = "t distribution"),
    linewidth = 1
  ) +
  scale_color_manual(
    name = "Densities",
    values = c("Normal distribution" = "blue", "t distribution" = "red")
  ) +
  labs(
    title = "Distribution of daily returns with normal and t density",
    x = "Daily log return in %",
    y = "Density"
  ) +
  theme_minimal()


############################################################
## 9. QQ plots: normal vs t
############################################################

# 9.1 QQ plot against normal distribution
qqnorm(ret, main = "QQ-plot returns vs normal distribution")
qqline(ret, col = 2)

# 9.2 QQ plot against t distribution
n = length(ret)
p = ppoints(n)

# theoretical quantiles of t distribution
q_t   = qt(p, df = df_t)
q_emp = sort((ret - mu) / sigma)  # standardized returns

plot(q_t, q_emp,
     xlab = "Theoretical t quantiles",
     ylab = "Empirical quantiles of standardized returns",
     main = "QQ-plot returns vs t distribution")
abline(0, 1, col = 2)


############################################################
## 10. Time series plots: price, returns, variance (from 2008)
############################################################

# 10.1 Daily data from 2008 onwards
wti_daily_2008 = wti_daily %>%
  filter(date >= as.Date("2008-01-01"))

# 10.2 Price evolution
ggplot(wti_daily_2008, aes(x = date, y = close)) +
  geom_line() +
  labs(
    title = "WTI price evolution (from 2008)",
    x = "Date",
    y = "Closing price"
  ) +
  theme_minimal()

# 10.3 Return series
ggplot(wti_daily_2008, aes(x = date, y = r_log)) +
  geom_line() +
  labs(
    title = "Daily log returns (from 2008)",
    x = "Date",
    y = "r_log in %"
  ) +
  theme_minimal()

# 10.4 Realized variance (log scale to see volatility clustering)
data_for_garch = data_for_garch %>%
  mutate(log_RV_5m = log(RV_5m))

ggplot(data_for_garch, aes(x = date, y = log_RV_5m)) +
  geom_line() +
  labs(
    title = "Log realized variance (5-minute based)",
    x = "Date",
    y = "log(RV_5m)"
  ) +
  theme_minimal()


############################################################
## Identify crisis, calm and transition periods
############################################################

# Base data frame with additional variables
df = data_for_garch %>%
  arrange(date) %>%
  mutate(
    abs_ret  = abs(r_log),
    big_move = abs_ret > 3    # |ret| > 3% = extreme move
  )


############################################################
## 11. Crisis periods: cluster days with |ret| > 3%
############################################################

max_gap_days = 30

crisis_days = df %>%
  filter(big_move) %>%
  arrange(date) %>%
  mutate(
    prev_date   = lag(date),
    gap_days    = as.numeric(date - prev_date),
    new_cluster = ifelse(is.na(gap_days) | gap_days > max_gap_days, 1L, 0L),
    cluster_id  = cumsum(new_cluster)
  )

crisis_clusters = crisis_days %>%
  group_by(cluster_id) %>%
  summarise(
    start_date   = min(date),
    end_date     = max(date),
    n_extreme    = n(),                                   # number of |ret|>3% days
    span_days    = as.numeric(max(date) - min(date)) + 1, # calendar span
    mean_abs_ret = mean(abs_ret),
    mean_RV_5m   = mean(RV_5m)
  ) %>%
  arrange(start_date)

crisis_clusters


############################################################
## 11a. Crisis period for modelling: Global Financial Crisis
############################################################

# pick the cluster whose start is in 2008/2009 (financial crisis)
financial_crisis = crisis_clusters %>%
  filter(start_date >= as.Date("2008-09-01"),
         start_date <= as.Date("2009-12-31")) %>%
  slice(1)

financial_crisis

crisis_start = financial_crisis$start_date
crisis_end   = financial_crisis$end_date

crisis_start
crisis_end

crisis_period = df %>%
  filter(date >= crisis_start, date <= crisis_end)


############################################################
## 12. Calm period (long low-vol window)
############################################################

calm_span = 300          # ~ 1 year of trading days
allowed_share_big = 0.03 # at most 3% of days with |ret| > 3%

# rolling mean of RV_5m and share of big moves in each window
roll_mean_RV = zoo::rollapply(
  df$RV_5m,
  width = calm_span,
  FUN   = mean,
  align = "left",
  fill  = NA
)

roll_share_big = zoo::rollapply(
  as.numeric(df$big_move),
  width = calm_span,
  FUN   = mean,
  align = "left",
  fill  = NA
)

calm_windows = df %>%
  mutate(
    roll_mean_RV   = roll_mean_RV,
    roll_share_big = roll_share_big,
    window_start   = date,
    window_end     = lead(date, calm_span - 1)
  ) %>%
  filter(
    !is.na(roll_mean_RV),
    roll_share_big <= allowed_share_big
  )

if (nrow(calm_windows) == 0) {
  stop("No calm windows found. Try reducing 'calm_span' or increasing 'allowed_share_big'.")
}

main_calm_window = calm_windows %>%
  arrange(roll_mean_RV) %>%
  slice(1)

calm_start = main_calm_window$window_start
calm_end   = main_calm_window$window_end

calm_start
calm_end

calm_period = df %>%
  filter(date >= calm_start, date <= calm_end)


############################################################
## 13. Transition period: calm to COVID crisis
############################################################

pre_days  = 200
post_days = 200

# find the cluster corresponding to COVID crisis (around Feb 2020)
covid_cluster = crisis_clusters %>%
  filter(start_date >= as.Date("2020-02-01"),
         start_date <= as.Date("2020-03-31")) %>%
  slice(1)

covid_cluster

transition_anchor = covid_cluster$start_date

transition_start = max(transition_anchor - pre_days, min(df$date))
transition_end   = min(transition_anchor + post_days, max(df$date))

transition_period = df %>%
  filter(date >= transition_start, date <= transition_end)

transition_start
transition_end


############################################################
## 14. Quick check of ranges and regime summary
############################################################

range(crisis_period$date)      # Financial crisis
range(calm_period$date)        # Calm
range(transition_period$date)  # Calm -> COVID crisis

ggplot(crisis_period, aes(x = date, y = log(RV_5m))) +
  geom_line() +
  labs(
    title = "Crisis period (Global Financial Crisis)",
    x = "Date",
    y = "log(RV_5m)"
  ) +
  theme_minimal()

ggplot(calm_period, aes(x = date, y = log(RV_5m))) +
  geom_line() +
  labs(
    title = "Calm period (long low-vol window)",
    x = "Date",
    y = "log(RV_5m)"
  ) +
  theme_minimal()

ggplot(transition_period, aes(x = date, y = log(RV_5m))) +
  geom_line() +
  labs(
    title = "Transition period (calm to COVID crisis)",
    x = "Date",
    y = "log(RV_5m)"
  ) +
  theme_minimal()

df %>%
  mutate(regime = case_when(
    date >= crisis_start     & date <= crisis_end     ~ "crisis",
    date >= calm_start       & date <= calm_end       ~ "calm",
    date >= transition_start & date <= transition_end ~ "transition",
    TRUE ~ "other"
  )) %>%
  filter(regime %in% c("crisis", "calm", "transition")) %>%
  group_by(regime) %>%
  summarise(
    mean_RV     = mean(RV_5m),
    median_RV   = median(RV_5m),
    mean_abs_ret = mean(abs(r_log))
  )


############################################################
## 15. GARCH helpers: common utilities
############################################################

library(rugarch)
library(forecast)

## 15.1 Prepare data: sort, select needed columns, drop NAs
prepare_garch_data = function(data) {
  data %>%
    arrange(date) %>%
    select(date, r_log, RV_5m) %>%
    filter(!is.na(r_log), !is.na(RV_5m))
}

## 15.2 Build rugarch specification for given model and distribution
build_garch_spec = function(model_type, dist_type) {
  ugarchspec(
    variance.model = list(
      model      = model_type,
      garchOrder = c(1, 1)
    ),
    mean.model = list(
      armaOrder    = c(0, 0),
      include.mean = TRUE
    ),
    distribution.model = dist_type
  )
}

############################################################
## 16. In-sample GARCH fit (single fit on full sample)
############################################################

garch_insample_helper = function(data,
                                 model_type = c("sGARCH", "eGARCH", "gjrGARCH"),
                                 dist_type  = c("norm", "std"),
                                 make_plots = TRUE) {
  
  model_type = match.arg(model_type)
  dist_type  = match.arg(dist_type)
  
  df = prepare_garch_data(data)
  
  # 1) Specification and fit
  spec = build_garch_spec(model_type, dist_type)
  fit  = ugarchfit(spec = spec, data = df$r_log, solver = "hybrid")
  
  # 2) In-sample variance and loss vs RV_5m
  sigma_vec  = as.numeric(sigma(fit))
  sigma2_vec = sigma_vec^2
  
  results = df %>%
    mutate(
      sigma  = sigma_vec,
      sigma2 = sigma2_vec,
      MSE    = (RV_5m - sigma2)^2,
      QLIKE  = log(sigma2) + RV_5m / sigma2
    )
  
  # 2a) Information criteria (infocriteria() liefert Matrix)
  ic     = infocriteria(fit)
  loglik = likelihood(fit)
  
  loss_summary = data.frame(
    n_obs       = nrow(results),
    MSE_mean    = mean(results$MSE,   na.rm = TRUE),
    QLIKE_mean  = mean(results$QLIKE, na.rm = TRUE),
    RV_mean     = mean(results$RV_5m, na.rm = TRUE),
    sigma2_mean = mean(results$sigma2, na.rm = TRUE),
    logLik      = as.numeric(loglik),
    AIC         = as.numeric(ic[1]),  # "Akaike"
    BIC         = as.numeric(ic[2])   # "Bayes"
  )
  
  # 3) Plot (optional)
  if (make_plots) {
    title_txt = paste0(
      "In-sample variance vs realized variance (",
      model_type, ", dist=", dist_type, ")"
    )
    
    p = ggplot(results, aes(x = date)) +
      geom_line(aes(y = RV_5m,  colour = "Realized variance (RV_5m)")) +
      geom_line(aes(y = sigma2, colour = "Fitted variance (GARCH)")) +
      labs(
        title = title_txt,
        x = "Date",
        y = "Variance"
      ) +
      scale_colour_manual(
        name   = "",
        values = c("Realized variance (RV_5m)" = "black",
                   "Fitted variance (GARCH)"   = "blue")
      ) +
      theme_minimal()
    
    print(p)
  }
  
  return(list(
    model_type   = model_type,
    dist_type    = dist_type,
    loss_summary = loss_summary,
    results      = results,
    fit_object   = fit
  ))
}

############################################################
## 17. Out-of-sample GARCH (rolling / expanding, 80/20 split)
##     → eigene Schleife, KEIN ugarchroll (um den Fehler zu vermeiden)
############################################################

garch_oos_helper = function(data,
                            model_type   = c("sGARCH", "eGARCH", "gjrGARCH"),
                            dist_type    = c("norm", "std"),
                            window_type  = c("rolling", "expanding"),
                            window_size  = 1000,      # nur für rolling
                            train_frac   = 0.8,
                            make_plots   = TRUE) {
  
  model_type  = match.arg(model_type)
  dist_type   = match.arg(dist_type)
  window_type = match.arg(window_type)
  
  df = prepare_garch_data(data)
  n  = nrow(df)
  
  if (n < 100) warning("Very few observations in this sample (n < 100).")
  
  # 17.1 Train/Test-Split
  n_train = floor(train_frac * n)
  n_test  = n - n_train
  
  if (n_test < 20) warning("Test sample is very small (n_test < 20).")
  if (n_test <= 0) stop("Not enough observations for out-of-sample split.")
  
  if (window_type == "rolling") {
    window_size = min(window_size, n_train)  # maximal Trainingslänge
  }
  
  spec   = build_garch_spec(model_type, dist_type)
  ret    = df$r_log
  
  dates_test = df$date[(n_train + 1):n]
  RV_test    = df$RV_5m[(n_train + 1):n]
  
  sigma_fc   = rep(NA_real_, n_test)
  
  # 17.2 Schleife über alle OOS-Punkte
  for (i in seq_len(n_test)) {
    
    if (window_type == "expanding") {
      idx_start = 1
      idx_end   = n_train + i - 1
    } else { # rolling window
      idx_end   = n_train + i - 1
      idx_start = max(1, idx_end - window_size + 1)
    }
    
    y_train = ret[idx_start:idx_end]
    
    fit_i = ugarchfit(spec = spec, data = y_train, solver = "hybrid")
    fc_i  = ugarchforecast(fit_i, n.ahead = 1)
    
    sigma_fc[i] = as.numeric(sigma(fc_i)[1])
  }
  
  # 17.3 Ergebnisse in ein Data Frame
  results = data.frame(
    date   = dates_test,
    sigma  = sigma_fc,
    sigma2 = sigma_fc^2
  )
  
  results = df %>%
    filter(date %in% dates_test) %>%
    select(date, r_log, RV_5m) %>%
    left_join(results, by = "date") %>%
    mutate(
      error = RV_5m - sigma2,
      MSE   = error^2,
      QLIKE = log(sigma2) + RV_5m / sigma2
    )
  
  loss_summary = results %>%
    summarise(
      n_obs       = n(),
      MSE_mean    = mean(MSE,   na.rm = TRUE),
      QLIKE_mean  = mean(QLIKE, na.rm = TRUE),
      RV_mean     = mean(RV_5m, na.rm = TRUE),
      sigma2_mean = mean(sigma2, na.rm = TRUE)
    )
  
  # 17.4 Plot (optional)
  if (make_plots) {
    plot_title = paste0(
      "OOS forecast vs realized variance (",
      model_type, ", ", window_type, ", dist=", dist_type, ")"
    )
    
    p = ggplot(results, aes(x = date)) +
      geom_line(aes(y = RV_5m,  colour = "Realized variance (RV_5m)")) +
      geom_line(aes(y = sigma2, colour = "Forecast variance (GARCH)")) +
      labs(
        title = plot_title,
        x = "Date",
        y = "Variance"
      ) +
      scale_colour_manual(
        name   = "",
        values = c("Realized variance (RV_5m)" = "black",
                   "Forecast variance (GARCH)" = "red")
      ) +
      theme_minimal()
    
    print(p)
  }
  
  return(list(
    model_type   = model_type,
    dist_type    = dist_type,
    window_type  = window_type,
    window_size  = if (window_type == "rolling") window_size else NA,
    train_frac   = train_frac,
    loss_summary = loss_summary,
    results      = results
  ))
}

############################################################
## 18. Run GARCH models: in-sample and out-of-sample
############################################################

# 18.1 Regimes
regimes = list(
  crisis     = crisis_period,
  calm       = calm_period,
  transition = transition_period
)

models  = c("sGARCH", "eGARCH", "gjrGARCH")
dists   = c("norm", "std")
windows = c("rolling", "expanding")

############################################################
## 18.2 In-sample results
############################################################

insample_results = list()
insample_summary = tibble()

for (reg_name in names(regimes)) {
  dat = regimes[[reg_name]]
  
  for (m in models) {
    for (d in dists) {
      
      cat("In-sample:", reg_name, m, d, "\n")
      
      res = garch_insample_helper(
        data       = dat,
        model_type = m,
        dist_type  = d,
        make_plots = FALSE
      )
      
      insample_results[[paste(reg_name, m, d, sep = "_")]] = res
      
      row = cbind(
        regime = reg_name,
        model  = m,
        dist   = d,
        res$loss_summary
      )
      
      insample_summary = dplyr::bind_rows(insample_summary, row)
    }
  }
}

insample_summary  # jetzt mit AIC/BIC gefüllt

############################################################
## 18.3 OOS results für alle Regimes / Modelle / Dists / Windows
############################################################

oos_results = list()
oos_summary = tibble()

for (reg_name in names(regimes)) {
  dat = regimes[[reg_name]]
  
  base_n    = nrow(prepare_garch_data(dat))
  ws_default = min(1000, floor(0.8 * base_n))
  
  for (m in models) {
    for (d in dists) {
      for (w in windows) {
        
        cat("OOS:", reg_name, m, d, w, "\n")
        
        res = garch_oos_helper(
          data        = dat,
          model_type  = m,
          dist_type   = d,
          window_type = w,
          window_size = ws_default,
          train_frac  = 0.8,
          make_plots  = FALSE
        )
        
        oos_results[[paste(reg_name, m, d, w, sep = "_")]] = res
        
        row = cbind(
          regime = reg_name,
          model  = m,
          dist   = d,
          window = w,
          res$loss_summary
        )
        
        oos_summary = dplyr::bind_rows(oos_summary, row)
      }
    }
  }
}

oos_summary

############################################################
## 19. Diebold–Mariano Test (Beispiel-Funktion)
############################################################

dm_test_garch = function(res1, res2) {
  df1 = res1$results
  df2 = res2$results
  
  df_join = df1 %>%
    select(date, RV_5m, sigma2_1 = sigma2, error_1 = error) %>%
    inner_join(
      df2 %>% select(date, RV_5m2 = RV_5m, sigma2_2 = sigma2, error_2 = error),
      by = "date"
    )
  
  if (!all(df_join$RV_5m == df_join$RV_5m2)) {
    stop("RV_5m does not match between the two models.")
  }
  
  e1 = df_join$error_1
  e2 = df_join$error_2
  
  forecast::dm.test(e1, e2, h = 1, power = 2)
}

# 7.2 Pairwise DM comparisons within a regime
# Vergleicht alle Modelle eines Regimes mit gleicher Verteilung & gleichem Window

compare_dm_regime = function(regime_name,
                             dist_filter   = "std",
                             window_filter = "rolling") {
  
  # Alle Keys in oos_results, die zu diesem Regime gehören
  keys_regime = names(oos_results)[grepl(paste0("^", regime_name, "_"), names(oos_results))]
  if (length(keys_regime) == 0) {
    stop("No OOS results found for regime: ", regime_name)
  }
  
  # Keys in Komponenten zerlegen: regime_model_dist_window
  info = do.call(rbind, lapply(keys_regime, function(k) {
    parts = strsplit(k, "_")[[1]]
    data.frame(
      key    = k,
      regime = parts[1],
      model  = parts[2],
      dist   = parts[3],
      window = parts[4],
      stringsAsFactors = FALSE
    )
  }))
  
  info = as_tibble(info) %>%
    dplyr::filter(
      dist   == dist_filter,
      window == window_filter
    )
  
  if (nrow(info) < 2) {
    warning("Not enough models in regime=", regime_name,
            " for dist=", dist_filter,
            " and window=", window_filter,
            " (need at least 2).")
    return(tibble())
  }
  
  dm_list = list()
  
  # Alle Modellpaare i<j im gewählten Regime / dist / window
  cnt = 0
  for (i in 1:(nrow(info) - 1)) {
    for (j in (i + 1):nrow(info)) {
      res1 = oos_results[[ info$key[i] ]]
      res2 = oos_results[[ info$key[j] ]]
      
      dm   = dm_test_garch(res1, res2)
      
      cnt  = cnt + 1
      dm_list[[cnt]] = tibble(
        regime = regime_name,
        dist   = dist_filter,
        window = window_filter,
        model_1 = info$model[i],
        model_2 = info$model[j],
        DM_stat = as.numeric(dm$statistic),
        p_value = dm$p.value
      )
    }
  }
  
  bind_rows(dm_list)
}




# 7.3 Beispiele: alle Modelle pro Regime (std, rolling)

dm_crisis_std_rolling     = compare_dm_regime("crisis",     dist_filter = "std", window_filter = "rolling")
dm_calm_std_rolling       = compare_dm_regime("calm",       dist_filter = "std", window_filter = "rolling")
dm_transition_std_rolling = compare_dm_regime("transition", dist_filter = "std", window_filter = "rolling")
dm_full_std_rolling       = compare_dm_regime("full",       dist_filter = "std", window_filter = "rolling")

dm_crisis_std_rolling
dm_calm_std_rolling
dm_transition_std_rolling
dm_full_std_rolling








