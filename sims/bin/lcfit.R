lcfit_bsm_ml_t <- function(m) {
  if (m$c > m$m) {
    ml_t <- ((log((m$c - m$m) / (m$c + m$m))) / (-m$r)) - m$b
  } else {
    ml_t <- Inf
  }

  if (is.finite(ml_t) && ml_t < 0) {
    ml_t <- 0.0
  }
  return(ml_t)
}

lcfit_bsm_log_like <- function(t, m) {
  expterm <- exp(-m$r * (t + m$b))
  ll <- m$c * log((1 + expterm) / 2) + m$m * log((1 - expterm) / 2)
  return(ll)
}

lcfit_bsm_regime <- function(m) {
  if (m$c == m$m) { return(-1) }
  if (m$c < m$m) { return(4) }

  lhs <- exp(m$b * m$r)
  rhs <- (sqrt(m$c) + sqrt(m$m))^2 / (m$c - m$m)

  if (lhs > rhs) { return(3) }
  if (m$b > 0.0) { return(2) }

  return(1)
}
