# tidy_data.R

suppressPackageStartupMessages(library(dplyr))

agg.noniter <- tbl_df(read.csv('runs/agg_maxima.csv', header = TRUE)) %>%
  select(everything(), -c, -m, -r, -b) %>%
  rename(t_ml = t) %>%
  mutate(tolerance = NA, n_eval = NA, method = "normal")

agg.iter <- tbl_df(read.csv('runs/agg_tol.csv', header = TRUE)) %>%
  select(everything(), -ml_brent, -n_eval_brent) %>%
  rename(t_ml = ml_t, t_hat = ml_est) %>%
  mutate(method = "iterative")

agg.brent <- tbl_df(read.csv('runs/agg_tol.csv', header = TRUE)) %>%
  select(everything(), -ml_est, -n_eval) %>%
  rename(t_ml = ml_t, t_hat = ml_brent, n_eval = n_eval_brent) %>%
  mutate(method = "brent")

# from lcfit_compare.cc:
# Limits are from Bio++: minimum branch length which may be considered is 1e-6

lcfit <- bind_rows(agg.noniter, agg.iter, agg.brent) %>%
  rename(model = model_name, rate.dist = rdist_name,
         rate = branch_length_rate) %>%
  mutate(t_hat = ifelse(t_hat < 1e-6, 1e-6, t_hat)) %>%
  mutate(rate = sprintf("Exp(mu=%.2f)", 1/rate), 
         err = t_ml - t_hat, abs.err = abs(err), rel.err = abs.err / t_ml,
         within.tol = (abs.err <= tolerance),
         success = !is.na(t_hat))
