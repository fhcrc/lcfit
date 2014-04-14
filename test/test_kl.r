#!/usr/bin/env Rscript

# Generate test data

library(entropy)

logsumexp <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}

kl_log <- function(log_freqs1, log_freqs2) {
  log_freqs1 <- log_freqs1 - logsumexp(log_freqs1)
  log_freqs2 <- log_freqs2 - logsumexp(log_freqs2)
  LR <- log_freqs1 - log_freqs2
  sum(exp(log_freqs1) * LR) / log(2)
}

pv <- function(v) paste('[', paste(v, collapse=', '), ']', sep='')

set.seed(1)
for(i in 1:4) {
  x <- log(runif(4) / 1000)
  y <- log(runif(4) / 1000)
  cat(sprintf("        x = %s\n        y = %s\n        expected = %f  # %f\n        self._test(x, y, expected)\n\n",
              pv(x),
              pv(y),
              kl_log(x, y),
              KL.plugin(exp(x) / sum(exp(x)),
                        exp(y) / sum(exp(y)), unit='log2')))
}
