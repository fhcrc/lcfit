# Compute log(sum(exp(x))).
# http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
#
# Note that sorting x in descending order may increase precision of the result.
logsumexp <- function(x) {
  #sort(x, decreasing = TRUE)
  #i <- 1
  i <- which.max(x)
  log1p(sum(exp(x[-i] - x[i]))) + x[i]
}

# Compute the logarithm of the sampling helper function K(x) given a model and
# an exponential prior by summing the series.
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
.k_exp_prior_series_ln <- function(m, lambda, x) {
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

# Integrand of the integral representation of Appell's F1 function with y = -x.
# http://functions.wolfram.com/HypergeometricFunctions/AppellF1/07/ShowAll.html
.xxx_appell_integrand <- function(t, a, b1, b2, c, x) {
  (1 - t)^(-1 - a + c) * t^(-1 + a) / ((1 - t*x)^b1 * (1 + t*x)^b2)
}

# Integrand of the integral representation of Appell's F1 function with y = -x.
# http://functions.wolfram.com/HypergeometricFunctions/AppellF1/07/ShowAll.html
.appell_integrand <- function(t, a, b1, b2, c, x) {
  y <- (-1 - a + c)*log(1 - t) + (-1 + a)*log(t) - b1*log(1 - t*x) - b2*log(1 + t*x)
  exp(y)
}

# Compute the logarithm of the sampling helper function K(x) given a model and
# an exponential prior by approximating the integral of Appell's F1 function.
# http://functions.wolfram.com/HypergeometricFunctions/AppellF1/07/ShowAll.html
.k_exp_prior_appell_ln <- function(m, lambda, x) {
  a <- lambda / m$r
  b1 <- -m$m
  b2 <- -m$c
  c <- (m$r + lambda) / m$r

  lv <- lgamma(c) - (lgamma(a) + lgamma(-a + c)) +
    log(integrate(.appell_integrand, 0, 1, a, b1, b2, c, x)$value)

  ly <- log(m$r) + (lambda / m$r) * log(x) - log(lambda) + lv
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
  lk0 <- lcfit_k_exp_prior_ln(m, lambda, exp(-m$r * m$b))
  vinv <- Vectorize(lcfit_inv_exp_prior_ln, vectorize.args = "lu")
  u <- runif(N)

  x <- vinv(m, lambda, log1p(-u) + lk0)
  t <- -1.0 / m$r * log(x) - m$b

  return(t)
}

# Compare branch lengths sampled from a model and an exponential prior with the
# expected distribution.
lcfit_sample_exp_prior_compare <- function(m, lambda, s, breaks = 100) {
  # Compute sample histogram.
  hobj <- hist(s, breaks = breaks, plot = FALSE)
  t <- hobj$mids

  # Compute expected likelihoods and (unnormalized) densities.
  x <- exp(-m$r*(m$b + t))
  ll <- m$c*log(1 + x) + m$m*log(1 - x) - lambda*t
  ld <- ll + log(m$r) - lambda*m$b

  # Normalize densities by approximating the function and dividing by its
  # integral.
  ldfun <- approxfun(t, ld)
  dfun <- function(t) { exp(ldfun(t)) }
  cons <- integrate(dfun, min(t), max(t))$value
  d <- exp(ld - log(cons))

  data <- data.frame(t = t, expected = d, observed = hobj$density)
}

#####

set.seed(0)

# Choose the implementation of K(x) to use.
lcfit_k_exp_prior_ln <- .k_exp_prior_appell_ln

lambda <- 0.1

m.s <- list(c = 5, m = 8, r = 1, b = 0.5)
print(system.time(samples.s <- lcfit_sample_exp_prior(m.s, lambda, 1000)))
data.s <- lcfit_sample_exp_prior_compare(m.s, lambda, samples.s)

p.s <- ggplot(data.s, aes(x = t)) + geom_line(aes(y = expected)) +
  geom_bar(aes(y = observed), stat = "identity", alpha = 0.4) +
  ylab('probability') +
  xlab('branch length')

m.m <- list(c = 1100, m = 800, r = 2, b = 0.5)
print(system.time(samples.m <- lcfit_sample_exp_prior(m.m, lambda, 1000)))
data.m <- lcfit_sample_exp_prior_compare(m.m, lambda, samples.m)

p.m <- p.s %+% data.m
