library(ggplot2)
library(gridExtra)

source("lcfit.R")

model.regime <- function(m) {
  if (m$c == m$m) { return(-1) }
  if (m$c < m$m) { return(4) }

  lhs <- exp(m$b * m$r)
  rhs <- (sqrt(m$c) + sqrt(m$m))^2 / (m$c - m$m)

  if (lhs > rhs) { return(3) }
  if (m$b > 0.0) { return(2) }

  return(1)
}

t <- seq(1e-6, 2, length=1000)

# regime 1: c > m, b = 0, exp(br) <= (sqrt(c) + sqrt(m))^2 / (c - m)
m1 <- list(c = 10, m = 1, r = 1, b = 0)
ll1 <- lcfit_bsm_log_like(t, m1)

# regime 2: c > m, b > 0, exp(br) <= (sqrt(c) + sqrt(m))^2 / (c - m)
m2 <- list(c = 10, m = 1, r = 1, b = 0.1)
ll2 <- lcfit_bsm_log_like(t, m2)

# regime 3: c > m, b > 0, exp(br) > (sqrt(c) + sqrt(m))^2 / (c - m)
m3 <- list(c = 10, m = 1, r = 1, b = 1)
ll3 <- lcfit_bsm_log_like(t, m3)

# regime 4: c < m
m4 <- list(c = 1, m = 10, r = 1, b = 0.1)
ll4 <- lcfit_bsm_log_like(t, m4)

regimes <- data.frame(t=t, ll1=ll1, ll2=ll2, ll3=ll3, ll4=ll4)

my_theme <- theme_bw() + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                               axis.text.y = element_blank(), axis.title.y = element_blank(),
                               axis.ticks = element_blank())

p1 <- ggplot(regimes, aes(x=t, y=ll1)) + geom_line() + my_theme
p2 <- ggplot(regimes, aes(x=t, y=ll2)) + geom_line() + my_theme
p3 <- ggplot(regimes, aes(x=t, y=ll3)) + geom_line() + my_theme
p4 <- ggplot(regimes, aes(x=t, y=ll4)) + geom_line() + my_theme

svg("regimes.svg")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()
