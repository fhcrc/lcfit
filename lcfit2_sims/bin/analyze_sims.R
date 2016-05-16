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

#
# load data
#

lnl <- tbl_df(read.csv("lnl.agg.csv"))
lcfit2 <- tbl_df(read.csv("lcfit2.agg.csv"))

#
# tidy data
# 

lnl_t <- lnl %>% 
  gather("distribution", "lnl", c(empirical, lcfit2))

# unite run parameters into a single key and drop redundant columns.
# currently n_sites is always 1000, n_leaves is always 10, and the seed is 
# included in source_tree.
lnl_tu <- lnl_t %>%
  unite(key, source_tree, model_name, rdist_name, branch_length_rate) %>%
  select(-c(n_sites, n_leaves, seed))

#
# plots
#

# normalized log-likelihood plots
lnl_gtu <- lnl_tu %>%
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
print(lnl_gtu$lnl_plot)
dev.off()
