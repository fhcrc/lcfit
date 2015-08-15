source("bin/lcfit.R")

# Compute log(sum(exp(x))).
# http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
#
# Note that sorting x in ascending order may increase precision of the result.
logsumexp <- function(x) {
  i <- which.max(x)
  log1p(sum(exp(x[-i] - x[i]))) + x[i]
}

# Compute the logarithm of the sampling helper function K(x) given a model and
# an exponential prior.
#
# K(x) = \sum_{i = 0}^{c} {
#          \sum_{j = 0}^{m} {
#            {{c}\choose{i}} {{m}\choose{j}} (i + j + \lambda / r)^{-1} (-1)^j x^{i + j + \lambda / r}
#          }
#        }
#
# Note that the inner sum is an alternating series. To improve numerical
# stability, the equation is rearranged so that sums can be delayed as long as
# possible by storing even terms and odd terms in separate vectors. The equation
# as implemented is
#
# K(x) = x^{\lambda / r} \left[
#          \sum_{i = 0}^{c} {
#            {{c}\choose{i}} \sum_{\substack{j = 0 \\ j \text{ even}}}^{m} {
#              {{m}\choose{j}} (i + j + \lambda / r)^{-1} x^{i + j}
#            }
#          } -
#          \sum_{i = 0}^{c} {
#            {{c}\choose{i}} \sum_{\substack{j = 1 \\ j \text{ odd}}}^{m} {
#              {{m}\choose{j}} (i + j + \lambda / r)^{-1} x^{i + j}
#            }
#          }
#        \right]
#
# For brevity, the variables below are abbreviated. For example,
#
#   lve, lvo: "log v even" and "log v odd", the log of the terms of the outer sums
#         lu: "log u", the log of the terms of the inner sums
# lsue, lsuo: "log sum u even" and "log sum u odd", the log of the sums of the even and odd inner terms
#
# etc.
lcfit_k_exp_prior_ln <- function(m, lambda, x) {
  # Allocate vectors for the outer sums (even and odd terms).
  lve <- vector(length = m$c+1)
  lvo <- vector(length = m$c+1)

  for (i in seq(0, m$c)) {
    j <- seq(0, m$m)

    # Compute all the terms of the inner sums. The even and odd terms are
    # separated below.
    lu <- lchoose(m$m, j) - log(i + j + lambda / m$r) + (i + j) * log(x)

    stopifnot(is.finite(lu))

    # Separate and sum the even and odd terms.
    lsue <- logsumexp(lu[c(TRUE, FALSE)]) # 0, 2, 4, ...
    lsuo <- logsumexp(lu[c(FALSE, TRUE)]) # 1, 3, 5, ...

    stopifnot(is.finite(lsue), is.finite(lsuo))

    # Multiply the even and odd sums by ${{c}\choose{i}}$ and store.
    lve[i+1] <- lchoose(m$c, i) + lsue
    lvo[i+1] <- lchoose(m$c, i) + lsuo
  }

  stopifnot(is.finite(lve), is.finite(lvo))

  # Sum the even and odd terms.
  lsve <- logsumexp(lve)
  lsvo <- logsumexp(lvo)

  stopifnot(is.finite(lsve), is.finite(lsvo))

  #print(sprintf("%.50f > %.50f: %s", lsve, lsvo, lsve > lsvo))
  stopifnot(lsve > lsvo)

  # Take the difference of the sums of the even and odd terms and finish up by multiplying by $x^{\lambda / r}$.
  lsv <- lsve + log1p(-exp(lsvo - lsve)) + (lambda / m$r) * log(x)

  stopifnot(is.finite(lsv))

  return(lsv)
}

# Find the value x such that K(x) ~= u given a model and an exponential prior.
lcfit_inv_exp_prior <- function(m, lambda, u, tolerance = 1e-6) {
  a <- 0.0
  b <- 1.0
  x <- 0.5
  al <- exp(lcfit_k_exp_prior_ln(m, lambda, x));

  while (abs(al - u) > tolerance && b - a > x * .Machine$double.eps) {
    if (al > u) {
      b <- x
    } else {
      a <- x
    }

    x <- (a + b) / 2.0
    al <- exp(lcfit_k_exp_prior_ln(m, lambda, x))
  }

  return(x)
}

# Find the value x such that log(K(x)) ~= log(u) given a model and an exponential prior.
lcfit_inv_exp_prior_ln <- function(m, lambda, lu, tolerance = 1e-6) {
  a <- 0.0
  b <- 1.0
  x <- 0.5
  lal <- lcfit_k_exp_prior_ln(m, lambda, x);

  while (abs(exp(lal) - exp(lu)) > tolerance && b - a > x * .Machine$double.eps) {
    if (lal > lu) {
      b <- x
    } else {
      a <- x
    }

    x <- (a + b) / 2.0
    lal <- lcfit_k_exp_prior_ln(m, lambda, x)
  }

  return(x)
}

# Generate N samples from the posterior on branch lengths given a model and an
# exponential prior.
lcfit_sample_exp_prior <- function(m, lambda, N) {
  k0 <- exp(lcfit_k_exp_prior_ln(m, lambda, exp(-m$r * m$b)))
  vinv <- Vectorize(lcfit_inv_exp_prior, vectorize.args = "u")
  u <- runif(N)

  x <- vinv(m, lambda, (1 - u) * k0)
  t <- -1.0 / m$r * log(x) - m$b

  return(t)
}

lcfit_sample_exp_prior_ln <- function(m, lambda, N) {
  lk0 <- lcfit_k_exp_prior_ln(m, lambda, exp(-m$r * m$b))
  vinv <- Vectorize(lcfit_inv_exp_prior_ln, vectorize.args = "lu")
  u <- runif(N)

  x <- vinv(m, lambda, log1p(-u) + lk0)
  t <- -1.0 / m$r * log(x) - m$b

  return(t)
}

# Generate N samples from the posterior on branch lengths given a model and an
# exponential prior and plot a comparison with the expected distribution.
lcfit_sample_exp_prior_compare <- function(m, lambda, N) {
  # Generate samples.
  s <- lcfit_sample_exp_prior(m, lambda, N)

  # Compute sample histogram.
  nbins <- 100
  hobj <- hist(s, breaks = nbins, plot = FALSE)
  t <- hobj$mids

  # Compute expected likelihoods and (unnormalized) densities.
  # x <- exp(-m$r*(m$b + t))
  # l <- (1 + x)^m$c * (1 - x)^m$m * exp(-lambda * t)
  l <- exp(lcfit_bsm_log_like(t, m) - lambda * t)
  d <- l * m$r / exp(lambda * m$b)

  # Normalize densities by approximating the function and dividing by its
  # integral.
  dfun <- approxfun(t, d)
  cons <- integrate(dfun, min(t), max(t))$value
  d <- d / cons

  data <- data.frame(t = t, expected = d, observed = hobj$density)

  ggplot(data, aes(x = t)) + geom_line(aes(y = expected)) +
    geom_bar(aes(y = observed), stat = "identity", alpha = 0.4) +
    ylab('probability') +
    xlab('branch length')
}

lambda <- 0.1

m.s <- list(c = 5, m = 8, r = 1, b = 0.5)
print(system.time(p.s <- lcfit_sample_exp_prior_compare(m.s, lambda, 1000)))

#m.m <- list(c = 100, m = 1, r = 1, b = 0.5)
#print(lcfit_sample_exp_prior_compare(m.m, lambda, 1000))
