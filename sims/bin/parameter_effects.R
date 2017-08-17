#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

theme_set(theme_bw(16))
theme_update(legend.position="bottom", legend.key=element_blank(),
             plot.title=element_text(rel(1.0), family = "Helvetica"))

bsm_loglike <- function(t, c, m, r, b) c*log((1+exp(-r*(t+b)))/2)+m*log((1-exp(-r*(t+b)))/2)

varying_title <- function(param)
  ggtitle(bquote(paste(italic(.(param)), .(" varying"))))
linetype <- function(param)
  scale_linetype_discrete(name=bquote(italic(.(param))))
xlabel <- xlab("branch length")
ylabel <- ylab("log likelihood")

x <- seq(0, 1, length.out=200)

r_values <- c(0.5, 1.0, 2.0, 4.0)
b_values <- c(0.25, 0.5, 1.0, 1.5)
m_values <- c(100, 250, 500, 1000)
c_values <- c(500, 1000, 2000, 5000)

b_varying <- expand.grid(t=x, b=b_values)

p_b <- ggplot(b_varying, aes(x=t, y=bsm_loglike(t, 1500, 1000, 4.0, b), linetype=factor(b))) +
  geom_line() +
  linetype('b') +
  xlabel +
  ylabel +
  annotate(geom="text", label="c=1500,m=1000,r=4.0", x=0.8, y=-1685) +
  varying_title('b')

# Plot varying R values
r_varying <- expand.grid(t=x, r=r_values)
p_r <- ggplot(r_varying, aes(x=t, y=bsm_loglike(t, 1500, 1000, r, .5), linetype=factor(r))) +
  geom_line() +
  linetype('r') +
  xlabel +
  ylabel +
  annotate(geom="text", label="c=1500,m=1000,b=0.5", x=0.8, y=-2350) +
  varying_title('r')

# Plot varying m values
m_varying <- expand.grid(t=x, m=m_values)
p_m <- ggplot(m_varying, aes(x=t, y=bsm_loglike(t, 1500, m, 2.0, .5), linetype=factor(m))) +
  geom_line() +
  linetype('m') +
  xlabel +
  ylabel +
  annotate(geom="text", label="c=1500,r=2.0,b=0.5", x=0.8, y=-1500) +
  varying_title('m')


c_varying <- expand.grid(t=x, c=c_values)
p_c <- ggplot(c_varying, aes(x=t, y=bsm_loglike(t, c, 1000, 2.0, .5), linetype=factor(c))) +
  geom_line() +
  linetype('c') +
  xlabel +
  ylabel +
  annotate(geom="text", label="m=1000,r=2.0,b=0.5", x=0.8, y=-2500) +
  varying_title('c')

pdf("parameter_effects.pdf", width=9.25, height=10)
grid.arrange(p_c, p_m, p_r, p_b, nrow=2)
dev.off()
