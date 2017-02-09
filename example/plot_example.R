library(dplyr)
library(tidyr)
library(ggplot2)

# this is a pared-down version of sims/bin/analyze_sims.R

#
# plot setup
#

linetypes <- c(empirical = "dashed",
               lcfit  = "solid")

colors <- c(empirical = "black",
            lcfit = "purple")

labels <- c(empirical = "empirical",
            lcfit = "lcfit")

#
# load data
#

lnl <- tbl_df(read.csv("lnl.csv"))

#
# tidy data
# 

lnl_t <- lnl %>% 
  gather("distribution", "lnl", c(empirical, lcfit))

#
# normalized log-likelihood curves
#

lnl_plots <- lnl_t %>%
  group_by(node_id) %>%
  do(lnl_plot = ggplot(., aes(x = t)) +
       geom_line(aes(y = lnl, color = distribution, linetype = distribution)) +
       ylab("normalized log-likelihood") +
       xlab("branch length") +
       scale_linetype_manual(values = linetypes, labels = labels) +
       scale_color_manual(values = colors, labels = labels) +
       theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
       ggtitle(sprintf("%s", first(.$node_id))))
  
pdf("curves.pdf")
print(lnl_plots$lnl_plot)
dev.off()
