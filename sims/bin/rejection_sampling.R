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

test_sampler <- function(sampler, m, lambda, N, label = "unknown") {
  message(paste("testing", label), ", model = (", m$c, ", ", m$m, ", ", m$r, ", ", m$b, ")")
  print(system.time(samples <- sampler(m, lambda, N)))
  data <- lcfit_sample_exp_prior_compare(m, lambda, samples)
  p <- ggplot(data, aes(x = t)) + geom_line(aes(y = expected)) +
    geom_bar(aes(y = observed), stat = "identity", alpha = 0.4) +
    ylab("probability") +
    xlab("branch length") +
    ggtitle(bquote(atop(.(label), list(c==.(m$c), m==.(m$m), r==.(m$r), b==.(m$b), lambda==.(lambda)))))

  return(list(model = m, lambda = lambda, samples = samples,
              data = data, label = label, plot = p))
}

compare_results <- function(results1, results2)
{
  m1 <- results1$model
  m2 <- results2$model
  stopifnot(m1$c == m2$c, m1$m == m2$m, m1$r == m2$r, m1$b == m2$b,
            results1$lambda == results2$lambda)

  data <- rbind(data.frame(t = results1$samples, method = results1$label),
                data.frame(t = results2$samples, method = results2$label))
  hobj <- hist(data$t, breaks = 100, plot = FALSE)

  p <- ggplot(data, aes(t, fill = method)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5,
                   breaks = hobj$breaks) +
    ylab("probability") +
    xlab("branch length") +
    ggtitle(bquote(list(c==.(m1$c), m==.(m1$m), r==.(m1$r), b==.(m1$b),
                        lambda==.(results1$lambda))))

  return(list(data = data, plot = p))
}

#####

source ("bin/log_sampling.R")

set.seed(0)

lambda <- 0.1
N <- 1000

# small model
m.s <- list(c = 5, m = 8, r = 1, b = 0.5)

inv.s <- test_sampler(lcfit_sample_exp_prior, m.s, lambda, N,
                      "inversion sampler")
rejR.s <- test_sampler(rejection_sampler, m.s, lambda, N,
                       "rejection sampler (R)")
rejC.s <- test_sampler(lcfit::lcfit_bsm_sample, m.s, lambda, N,
                       "rejection sampler (C++)")

inv_vs_rejR.s <- compare_results(inv.s, rejR.s)
inv_vs_rejC.s <- compare_results(inv.s, rejC.s)
rejR_vs_rejC.s <- compare_results(rejR.s, rejC.s)

# medium model
m.m <- list(c = 1100, m = 800, r = 2, b = 0.5)

inv.m <- test_sampler(lcfit_sample_exp_prior, m.m, lambda, N,
                      "inversion sampler")
rejR.m <- test_sampler(rejection_sampler, m.m, lambda, N,
                       "rejection sampler (R)")
rejC.m <- test_sampler(lcfit::lcfit_bsm_sample, m.m, lambda, N,
                       "rejection sampler (C++)")

inv_vs_rejR.m <- compare_results(inv.m, rejR.m)
inv_vs_rejC.m <- compare_results(inv.m, rejC.m)
rejR_vs_rejC.m <- compare_results(rejR.m, rejC.m)

# long model
m.l <- list(c = 2000, m = 500, r = 2, b = 0.5)
inv.l <- test_sampler(lcfit_sample_exp_prior, m.l, lambda, N,
                      "inversion sampler")
rejR.l <- test_sampler(rejection_sampler, m.l, lambda, N,
                       "rejection sampler (R)")
rejC.l <- test_sampler(lcfit::lcfit_bsm_sample, m.l, lambda, N,
                       "rejection sampler (C++)")

inv_vs_rejR.l <- compare_results(inv.l, rejR.l)
inv_vs_rejC.l <- compare_results(inv.l, rejC.l)
rejR_vs_rejC.l <- compare_results(rejR.l, rejC.l)
