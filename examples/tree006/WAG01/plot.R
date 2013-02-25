#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

theme_fullframe <- function (base_size = 12){
  structure(list(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0, "lines"),
        axis.ticks.margin = unit(0, "lines"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin = unit(0, "lines"),
        plot.background = element_blank(),
        plot.margin = unit(0*c(-1.5, -1.5, -1.5, -1.5), "lines")
        ), class = "theme")
}

theme_set(theme_bw(19))

bls <- read.csv('lcfit_bls.csv', as.is=TRUE)

node_plot <- function(node) {
  p <- ggplot(subset(bls, node_id == node),
              aes(x=branch_length, y=bpp_ll)) +
    geom_line() +
    xlab("Branch Length") +
    ylab("LnL")
}

svg('node3_node6.svg', width=14, height=7)
grid.arrange(node_plot(3), node_plot(6), nrow=1)
dev.off()

# All nodes
svg('all_nodes.svg', width=6, height=6)
nodes <- as.list(unique(bls$node_id)[1:16])
node_plots <- lapply(nodes, function(node) {
  ggplot(subset(bls, node_id == node),
              aes(x=branch_length, y=bpp_ll)) +
    geom_line() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.ticks=element_blank())
})
node_plots$nrow <- 4
do.call('grid.arrange', node_plots)
dev.off()
