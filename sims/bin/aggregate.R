#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
library(plyr)
library(RJSONIO)
suppressPackageStartupMessages(library(distrEx))
suppressPackageStartupMessages(library(entropy))

theme_set(theme_bw(16) + theme(strip.background=element_blank(),
                               legend.key = element_blank(),
                               legend.position='bottom'))

args <- commandArgs(TRUE)
stopifnot(length(args) >= 1)
control_paths <- args

controls <- llply(control_paths, fromJSON)

ess <- function(ll) { ll <- ll - max(ll); exp(-log(sum(exp(2 * ll))) + 2 * log(sum(exp(ll)))) }

logsumexp <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}

kl_log <- function(log_freqs1, log_freqs2) {
  log_freqs1 <- log_freqs1 - logsumexp(log_freqs1)
  log_freqs2 <- log_freqs2 - logsumexp(log_freqs2)
  LR <- log_freqs1 - log_freqs2
  sum(exp(log_freqs1) * LR) / log(2)
}

compare_bls <- function(d) {
  branch_length <- d[['branch_length']]
  bpp_ll <- d[['bpp_ll']]
  fit_ll <- d[['fit_ll']]

  bpp_l <- exp(bpp_ll - max(bpp_ll))
  bpp_l <- bpp_l / sum(bpp_l)
  fit_l <- exp(fit_ll - max(fit_ll))
  fit_l <- fit_l / sum(fit_l)

  bpp_d <- DiscreteDistribution(seq_along(branch_length), prob = bpp_l)
  fit_d <- DiscreteDistribution(seq_along(branch_length), prob = fit_l)
  data.frame(hellinger=HellingerDist(bpp_d, fit_d),
             kl=KL.plugin(bpp_l, fit_l, unit = 'log2'),
             my_kl=kl_log(bpp_ll, fit_ll),
             bpp_rel_ess=ess(bpp_ll) / length(bpp_ll),
             fit_rel_ess=ess(fit_ll) / length(bpp_ll))
}

error <- ldply(controls, function(p) {
  bls_file <- p$lcfit[1]
  maxima_file <- p$lcfit[2]
  fit_file <- p$lcfit[3]

  bls <- read.csv(bls_file, as.is=TRUE)
  maxima <- read.csv(maxima_file, as.is=TRUE)
  error <- ddply(bls, .(node_id), function(d) compare_bls(d))
  error <- join(error, maxima, by=c('node_id'), type='left', match='first')
  transform(error,
            model=p$model_name,
            rate=p$rdist_name,
            initial=p$initial,
            n_leaves=p$n_leaves,
            seed=p$seed,
            tree=p$source_tree,
            branch_length_rate=p$branch_length_rate,
            rel_err=abs(t - t_hat) / t,
            abs_err=abs(t - t_hat))
}, .progress = 'text')

p1 <- ggplot(error, aes(x=model, y=hellinger, fill=rate)) +
  geom_boxplot() +
  facet_grid(seed ~ branch_length_rate) +
  ggtitle("Hellinger distance between lcfit and empirical likelihoods
at sampled branch lengths") +
  theme(legend.position='bottom')

p2 <- ggplot(error, aes(x=model, y=my_kl, fill=rate)) +
  geom_boxplot() +
  facet_grid(seed ~ branch_length_rate) +
  ggtitle("KL divergence between lcfit and empirical likelihoods
at sampled branch lengths") +
  theme(legend.position='bottom') +
  ylab("KL divergence (bits)")

p3 <- ggplot(subset(error, branch_length_rate == 10), aes(x=t, y=my_kl, color=model)) +
  facet_grid(rate ~ model) +
  geom_point() +
  ylab('KL') +
  xlab('ML branch length')
ggsave('ml_branch_length_vs_kl.png', p3, width=14, height=9, dpi=72)

local({
  error <- na.omit(error)
  error_noshort <- subset(error, t >= 0.01)
  print(quantile(error$rel_err * 100))
  print(quantile(error_noshort$rel_err * 100))
  p1 <- ggplot(error_noshort, aes(x=t, y=rel_err, color=factor(branch_length_rate))) +
    facet_grid(rate ~ model) +
    geom_hline(yintercept=1.0, linetype='dashed') +
    geom_point() +
    ylab(expression(abs(frac(hat(t) - t, t)))) +
    scale_y_log10() +
    xlab('ML branch length') +
    scale_color_discrete(name = 'Branch length rate') +
    theme(legend.position = 'bottom') +
    ggtitle(expression(paste("Relative error: ", t >= 0.01)))
  p2 <- ggplot(error, aes(x=t, y=t_hat, color=factor(branch_length_rate))) +
    facet_grid(rate ~ model) +
    geom_abline(color='grey', linetype='dashed') +
    geom_point() +
    coord_equal() +
    ylab(expression(hat(t))) +
    xlab(expression(t)) +
    scale_color_discrete(name = 'Branch length rate') +
    theme(legend.position = 'bottom') +
    ggtitle("Comparison of fit branch lengths")

  svg('ml_branch_length_vs_rel_err.svg', width=14, height=12)
  tryCatch(grid.arrange(p1, p2, nrow=2), finally=dev.off())
})

p5 <- ggplot(error, aes(x=t, y=t_hat, color=factor(branch_length_rate))) +
  facet_grid(rate ~ model) +
  geom_point() +
  geom_abline(slope=1) +
  ylab(expression(hat(t))) +
  scale_y_log10() +
  scale_x_log10() +
  xlab('ML branch length') +
  theme(legend.position = 'bottom')
ggsave('t_vs_that.svg', p5, width=14, height=9, dpi=72)

ggsave('hellinger.png', p1, height=14, width=7, dpi=72)
ggsave('kl.png', p2, height=14, width=7, dpi=72)
pdf('kl_hellinger.pdf', width=21, height=14, useDingbats=FALSE)
tryCatch(grid.arrange(p1, p2, p3, ncol=3), finally=dev.off())
write.csv(error, 'agg.csv', row.names=FALSE)
