#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(RJSONIO)
suppressPackageStartupMessages(library(distrEx))
suppressPackageStartupMessages(library(entropy))

theme_set(theme_bw(16) + theme(strip.background=element_blank()))

args <- commandArgs(TRUE)
stopifnot(length(args) >= 2)
control_paths <- args[-1]

controls <- llply(control_paths, fromJSON)

hellinger_bls <- function(d) {
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
             raw_rss=sum((bpp_ll - fit_ll)^2))
}

hellinger <- ldply(controls, function(p) {
  bls_file <- p$lcfit[1]
  maxima_file <- p$lcfit[2]
  fit_file <- p$lcfit[3]

  bls <- read.csv(bls_file, as.is=TRUE)
  hellinger <- ddply(bls, .(node_id), function(d) hellinger_bls(d))
  transform(hellinger,
            model=p$model_name,
            rate=p$rdist_name,
            initial=p$initial,
            n_leaves=p$n_leaves,
            seed=p$seed,
            tree=p$source_tree,
            branch_length_rate=p$branch_length_rate)
})

p1 <- ggplot(hellinger, aes(x=model, y=hellinger, fill=rate)) +
  geom_boxplot() +
  facet_grid(seed ~ branch_length_rate) +
  ggtitle("Hellinger distance between lcfit and empirical likelihoods
at sampled branch lengths") +
  theme(legend.position='bottom')
p2 <- ggplot(hellinger, aes(x=model, y=raw_rss, fill=rate)) +
  geom_boxplot() +
  facet_grid(seed ~ branch_length_rate) +
  ggtitle("RSS between lcfit and empirical log-likelihoods
at sampled branch lengths") +
  theme(legend.position='bottom')

p3 <- ggplot(hellinger, aes(x=model, y=kl, fill=rate)) +
  geom_boxplot() +
  facet_grid(seed ~ branch_length_rate) +
  ggtitle("KL divergence between lcfit and empirical likelihoods
at sampled branch lengths") +
  theme(legend.position='bottom') +
  ylab("KL divergence (bits)")

ggsave('hellinger.png', p1, height=14, width=7, dpi=72)
ggsave('rss.png', p2, height=14, width=7, dpi=72)
ggsave('kl.png', p3, height=14, width=7, dpi=72)
