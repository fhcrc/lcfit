#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(scales)


# setup theme
theme_set(theme_bw(20))
theme_update(legend.key=element_blank())

# From http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                     mean = mean   (xx[,col], na.rm=na.rm),
                     sd   = sd     (xx[,col], na.rm=na.rm)
                     )
                 },
                 measurevar,
                 na.rm
                 )

  # Rename the "mean" column
  datac <- rename(datac, c("mean"=measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

maxima <- read.csv('lcfit_maxima.csv', as.is=TRUE)

m <- melt(maxima, id.vars=c('node_id', 'model_name', 'ml_tree', 'tolerance'), measure.vars=c('brent_n', 'lcfit_n'))
m <- transform(m, variable=ifelse(variable=='brent_n', 'Brent', 'lcfit'))


# Summarize
m_summary <- summarySE(m, "value", c("model_name", "tolerance", "variable"))

# Plot
xscale <- scale_x_continuous(trans=log10_trans(),
                            breaks=trans_breaks("log10", function(x) 10^x),
                            labels=trans_format("log10", math_format(10^.x)))

p <- ggplot(m_summary, aes(x=tolerance, y=value, ymin=value-se, ymax=value+se,
                           color=variable, shape=variable)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(width=0.15) +
  xscale +
  scale_color_discrete(name="Fit method") +
  scale_shape_discrete(name="Fit method") +
  theme(legend.position="bottom") +
  xlab("Tolerance") +
  ylab("Number of LL calls") +
  ylim(0, max(m_summary$value) + 2)
ggsave('ml_steps_by_tolerance.svg')
ggsave('ml_steps_by_tolerance.png')


p <- ggplot(m, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(x='Fit method', y='# of LL calls') +
theme(legend.position='none') +
  facet_wrap(~model_name) +
  ggsave('ml_lcfit_nsteps.svg')

p <- ggplot(maxima, aes(x=lcfit_t, y=brent_t, color=model_name)) +
  geom_abline(slope=1, intercept=0) +
  geom_point(size=1.2) +
  ggtitle("ML branch length estimates") +
  facet_wrap(~model_name) +
  theme(legend.position='none')

ggsave('ml_lcfit_comparison.svg')
