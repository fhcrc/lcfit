library(dplyr)
library(tidyr)
library(ggplot2)

#
# plot setup
#

linetypes <- c(empirical = "dashed",
               gamma = "solid",
               weibull = "solid",
               lcfit  = "solid",
               lcfit2 = "solid")

colors <- c(empirical = "black",
            gamma = "green",
            weibull = "red",
            lcfit = "purple",
            lcfit2 = "orange")

labels <- c(empirical = "empirical",
            gamma = "gamma",
            weibull = "weibull",
            lcfit = "lcfit",
            lcfit2 = "lcfit2")

plot_curves <- FALSE

#
# load data
#

# don't load the data again if the variables are already defined
if (!exists("lnl")) { lnl <- tbl_df(read.csv("lnl.agg.csv")) }
if (!exists("lcfit")) { lcfit <- tbl_df(read.csv("lcfit.agg.csv")) }

#
# tidy data
# 

lnl_t <- lnl %>% 
  gather("distribution", "lnl", c(empirical, lcfit))

#
# compute measures
#

kl_divergence <- function(p, q) {
  sum(p * log2(p / q))
}

hellinger_distance <- function(p, q) {
  (1 / sqrt(2)) * sqrt(sum( (sqrt(p) - sqrt(q))^2 ))
}

# note that the log-likelihoods are already normalized
freqs <- lnl %>%
  group_by(source_tree, model_name, rdist_name, branch_length_rate, node_id) %>%
  select(-c(n_sites, n_leaves, seed)) %>%
  mutate(empirical = exp(empirical),
         lcfit = exp(lcfit)) %>%
  mutate(empirical = empirical / sum(empirical),
         lcfit = lcfit / sum(lcfit))

measures <- freqs %>%
  summarize(kl = kl_divergence(empirical, lcfit),
            hellinger = hellinger_distance(empirical, lcfit))

#
# normalized log-likelihood curves
#

if (plot_curves) {
  # unite run parameters into a single key and drop redundant columns.
  # currently n_sites is always 1000, n_leaves is always 10, and the seed is 
  # included in source_tree.
  lnl_tu <- lnl_t %>%
    unite(key, source_tree, model_name, rdist_name, branch_length_rate) %>%
    select(-c(n_sites, n_leaves, seed))

  lnl_plots <- lnl_tu %>%
    group_by(key, node_id) %>%
    do(lnl_plot = ggplot(., aes(x = t)) +
         geom_line(aes(y = lnl, color = distribution, linetype = distribution)) +
         ylab("normalized log-likelihood") +
         xlab("branch length") +
         scale_linetype_manual(values = linetypes, labels = labels) +
         scale_color_manual(values = colors, labels = labels) +
         theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
         ggtitle(sprintf("%s %s", first(.$key), first(.$node_id))))

  pdf("curves.pdf")
  print(lnl_plots$lnl_plot)
  dev.off()
}

#
# distribution measures
#

# adjust branch length rate labels to show mu (scale) instead of lambda (rate)
measures$branch_length_rate <- factor(measures$branch_length_rate,
                                      labels = c(bquote(mu==0.1),
                                                 bquote(mu==0.01)))

p.measure <- ggplot(measures, aes(x = model_name, fill = model_name)) +
  xlab("model") +
  facet_grid(. ~ branch_length_rate, labeller = label_parsed) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")

# compute the upper limits of the boxplot whiskers for each rate distribution
upper_limits <- measures %>%
  group_by(rdist_name, branch_length_rate, model_name) %>%
  summarize(kl = boxplot.stats(kl)$stats[5],
            hellinger = boxplot.stats(hellinger)$stats[5]) %>%
  group_by(rdist_name) %>%
  summarize(top_kl = max(kl),
            top_hellinger = max(hellinger))

measure_plots <- measures %>%
  group_by(rdist_name) %>%
  inner_join(upper_limits) %>%
  do(kl = p.measure %+% . +
       geom_boxplot(aes(y = kl)) +
       coord_cartesian(ylim = 1.05 * c(0.0, .$top_kl[1])) +
       ylab("KL divergence (bits)"),
     hellinger = p.measure %+% . +
       geom_boxplot(aes(y = hellinger)) +
       coord_cartesian(ylim = 1.05 * c(0.0, .$top_hellinger[1])) +
       ylab("Hellinger distance"))

# summarize the outliers above the thresholds
upper_limits.t <- upper_limits %>%
  rename(kl = top_kl, hellinger = top_hellinger) %>%
  gather(measure, limit, kl, hellinger)

outliers <- measures %>%
  gather(measure, value, kl, hellinger) %>%
  group_by(rdist_name, branch_length_rate, measure) %>%
  inner_join(upper_limits.t) %>%
  filter(value > limit) %>%
  summarize(count = n(),
            lower_threshold = first(limit),
            mean = mean(value),
            median = median(value),
            max = max(value)) %>%
  ungroup() %>%
  arrange(measure, rdist_name, branch_length_rate)

# asymptotic error
# GOTCHA: there's one data point for which the asymptotic error is 
# about -3500, so we filter that one out for plotting
p.err <- p.measure %+% filter(lcfit, abs(err_max_t) < 3000) +
  geom_boxplot(aes(y = err_max_t)) +
  ylab("asymptote error (nats)")

#pdf("kl.pdf")
#print(measure_plots$kl)
#dev.off()

pdf("kl-gamma.pdf")
print(measure_plots$kl[[1]])
dev.off()

pdf("kl-uniform.pdf")
print(measure_plots$kl[[2]])
dev.off()

pdf("hellinger-gamma.pdf")
print(measure_plots$hellinger[[1]])
dev.off()

pdf("hellinger-uniform.pdf")
print(measure_plots$hellinger[[2]])
dev.off()

pdf("asymptote_err.pdf")
print(p.err)
dev.off()

# write KL outliers table to file
outliers.kl <- filter(outliers, measure == "kl") %>%
  select(-measure)
colnames(outliers.kl) <- c("Rate Distribution", "Branch Length Mean", "Count", "Plot Threshold", "Mean", "Median", "Maximum")
print.xtable(xtable(outliers.kl, label = "TABLEkloutliers", digits = 4),
             file = "kl-outliers.texi")
