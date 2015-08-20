source("bin/lcfit.R")

rejection_sampler <- function(m, lambda, N) {
  ml_t <- lcfit_bsm_ml_t(m)
  ml_ll <- lcfit_bsm_log_like(ml_t, m)

  if (is.finite(ml_t)) {
    f <- function(t) {
      exp(lcfit_bsm_log_like(t, m) + log(dexp(t, lambda)) - ml_ll - log(dexp(ml_t, lambda)))
    }
    g <- function(t) {
      exp(log(dexp(t, lambda)) - log(dexp(ml_t, lambda)))
    }
  } else {
    f <- function(t) { exp(lcfit_bsm_log_like(t, m) - ml_ll) }
    g <- function(t) { 1.0 }
  }

  s <- vector(length = N)
  i <- 0
  trials <- 0

  pb <- txtProgressBar(min = 0, max = N, style = 3)

  while (i <= N) {
    t <- rexp(1, lambda)
    u <- runif(1)

    trials <- trials + 1
    if (u * g(t) < f(t)) {
      i <- i + 1
      s[i] <- t
      setTxtProgressBar(pb, i)
    }
  }

  close(pb)
  print(sprintf("acceptance ratio: %g", i / trials))
  return(s)
}

#####

source ("bin/log_sampling.R")
library(dplyr)

print(system.time(rej.samples.s <- rejection_sampler(m.s, lambda, 1000)))
rej.data.s <- lcfit_sample_exp_prior_compare(m.s, lambda, rej.samples.s)
p.rej.s <- p.s %+% rej.data.s + ggtitle("rejection sampling, small model")

print(system.time(rej.samples.m <- rejection_sampler(m.m, lambda, 1000)))
rej.data.m <- lcfit_sample_exp_prior_compare(m.m, lambda, rej.samples.m)
p.rej.m <- p.s %+% rej.data.m + ggtitle("rejection sampling, medium model")

# m.w <- list(c = 2000, m = 500, r = 2, b = 0.4)
# print(system.time(rej.samples.w <- rejection_sampler(m.w, lambda, 1000)))
# rej.data.w <- lcfit_sample_exp_prior_compare(m.w, lambda, rej.samples.w)
# p.rej.w <- p.s %+% rej.data.w + ggtitle("rejection sampling, weird model")

combined.s <- tbl_df(rbind(data.frame(t = samples.s, method = "exact"),
                           data.frame(t = rej.samples.s, method = "rejection")))

p.all.s <- ggplot(combined.s, aes(t, fill = method)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5)

combined.m <- tbl_df(rbind(data.frame(t = samples.m, method = "exact"),
                           data.frame(t = rej.samples.m, method = "rejection")))

p.all.m <- p.all.s %+% combined.m
