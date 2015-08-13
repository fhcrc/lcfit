K <- function(m, lambda, x) {
  ve <- vector(length = m$c+1)
  vo <- vector(length = m$c+1)

  for (i in seq(0, m$c)) {
    j <- seq(0, m$m)
    u <- choose(m$m, j) / (i + j + lambda / m$r) * x^(i + j)

    stopifnot(is.finite(u))

    sue <- sum(u[c(TRUE, FALSE)]) # 0, 2, 4, ...
    suo <- sum(u[c(FALSE, TRUE)]) # 1, 3, 5, ...

    stopifnot(is.finite(sue), is.finite(suo))

    ve[i+1] <- choose(m$c, i) * sue
    vo[i+1] <- choose(m$c, i) * suo
  }

  stopifnot(is.finite(ve), is.finite(vo))

  sve <- sum(ve)
  svo <- sum(vo)

  print(sprintf("%f > %f: %s", sve, svo, sve > svo))
  stopifnot(sve > svo)

  sv <- (sve - svo) * x^(lambda / m$r)

  stopifnot(is.finite(sv))

  return(sv)
}

K.wrap <- function(m, lambda, t) {
  x <- exp(-m$r * (m$b + t))
  K(m, lambda, x)
}
