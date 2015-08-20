library(ggplot2)
source("bin/lcfit.R")

rejection_sampler <- function(m, lambda, N) {
  ml_t <- lcfit_bsm_ml_t(m)
  ml_ll <- lcfit_bsm_log_like(ml_t, m)

  f <- function(t) { exp(lcfit_bsm_log_like(t, m) - ml_ll) }

  i <- 0
  trials <- 0

  pb <- txtProgressBar(min = 0, max = N, style = 3)

#   # iterative implementation, for reference
#   s <- vector(length = N)
#   while (i <= N) {
#     t <- rexp(1, lambda)
#     u <- runif(1)
#
#     trials <- trials + 1
#     if (u < f(t)) {
#       i <- i + 1
#       s[i] <- t
#       setTxtProgressBar(pb, i)
#     }
#   }

  # vectorized implementation of the above
  s <- NULL
  batch_size <- min(N, 10000)
  while (i <= N) {
    t <- rexp(batch_size, lambda)
    u <- runif(batch_size)

    trials <- trials + batch_size
    t <- t[u < f(t)]
    s <- c(s, t)

    i <- i + length(t)
    setTxtProgressBar(pb, min(N, i))
  }

  s <- s[1:N]

  close(pb)
  print(sprintf("acceptance ratio: %g", i / trials))
  return(s)
}

#####

source ("bin/log_sampling.R")

# small model
print(system.time(rej.samples.s <- rejection_sampler(m.s, lambda, 1000)))
rej.data.s <- lcfit_sample_exp_prior_compare(m.s, lambda, rej.samples.s)
p.rej.s <- p.s %+% rej.data.s + ggtitle("rejection sampling, small model")

# medium model
print(system.time(rej.samples.m <- rejection_sampler(m.m, lambda, 1000)))
rej.data.m <- lcfit_sample_exp_prior_compare(m.m, lambda, rej.samples.m)
p.rej.m <- p.s %+% rej.data.m + ggtitle("rejection sampling, medium model")

# # weird model
# print(system.time(rej.samples.w <- rejection_sampler(m.w, lambda, 1000)))
# rej.data.w <- lcfit_sample_exp_prior_compare(m.w, lambda, rej.samples.w)
# p.rej.w <- p.s %+% rej.data.w + ggtitle("rejection sampling, weird model")

# small model, combined
combined.s <- rbind(data.frame(t = samples.s, method = "exact"),
                    data.frame(t = rej.samples.s, method = "rejection"))
hist.s <- hist(combined.s$t, breaks = 100, plot = FALSE)
p.all.s <- ggplot(combined.s, aes(t, fill = method)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5,
                 breaks = hist.s$breaks)

# medium model, combined
combined.m <- rbind(data.frame(t = samples.m, method = "exact"),
                    data.frame(t = rej.samples.m, method = "rejection"))
hist.m <- hist(combined.m$t, breaks = 100, plot = FALSE)
p.all.m <- ggplot(combined.m, aes(t, fill = method)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5,
                 breaks = hist.m$breaks)
