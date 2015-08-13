logsumexp <- function(x) {
  i <- which.max(x)
  log1p(sum(exp(x[-i] - x[i]))) + x[i]
}

LK <- function(m, lambda, x) {
  ls <- -.Machine$double.xmax

  lve <- vector(length = m$c+1)
  lvo <- vector(length = m$c+1)

  for (i in seq(0, m$c)) {
    j <- seq(0, m$m)
    lu <- lchoose(m$m, j) - log(i + j + lambda / m$r) + (i + j) * log(x)

    stopifnot(is.finite(lu))

    lsue <- logsumexp(lu[c(TRUE, FALSE)]) # 0, 2, 4, ...
    lsuo <- logsumexp(lu[c(FALSE, TRUE)]) # 1, 3, 5, ...

    stopifnot(is.finite(lsue), is.finite(lsuo))

    lve[i+1] <- lchoose(m$c, i) + lsue
    lvo[i+1] <- lchoose(m$c, i) + lsuo
  }

  stopifnot(is.finite(lve), is.finite(lvo))

  lsve <- logsumexp(lve)
  lsvo <- logsumexp(lvo)

  stopifnot(is.finite(lsve), is.finite(lsvo))

  print(sprintf("%.50f > %.50f: %s", lsve, lsvo, lsve > lsvo))
  stopifnot(lsve > lsvo)

  lsv <- lsve + log1p(-exp(lsvo - lsve)) + (lambda / m$r) * log(x)

  stopifnot(is.finite(lsv))

  return(lsv)
}

LK.wrap <- function(m, lambda, t) {
  x <- exp(-m$r * (m$b + t))
  LK(m, lambda, x)
}
