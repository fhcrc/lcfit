        aggfile <- "../../test.csv"
 rm(list = ls())     # clear objects  
 graphics.off()      # close graphics windows   

library(ggplot2)
library(gridExtra)

#create dummy data
test= data.frame(
  group = c("Group 1", "Group 1", "Group 1","Group 2", "Group 2", "Group 2"), 
  x = c(1 ,2,3,1,2,3 ),
  y = c(33,25,27,36,23,25),
  n=c(71,55,65,58,65,58),
  ypos=c(18,18,18,17,17,17)

  )


p1 <- qplot(x=x, y=y, data=test, colour=group) +
  ylab("Mean change from baseline") + 
  theme(plot.margin = unit(c(1,3,8,1), "lines")) +
  geom_line()+
  scale_x_continuous("Visits", breaks=seq(-1,3) ) +
  theme(legend.position="bottom",
       legend.title=element_blank()) 

# Create the textGrobs 
for (ii in 1:nrow(test))
{
  #display numbers at each visit
  p1=p1+ annotation_custom(grob = textGrob(test$n[ii]),  
                           xmin = test$x[ii], 
                           xmax = test$x[ii], 
                           ymin = test$ypos[ii], 
                           ymax = test$ypos[ii])
  }




  # Code to override clipping
  tmp <- ggplot_build(p1)
  gt <- ggplot_gtable(ggplot_build(p1))
  gt$layout$clip[gt$layout$name=="panel"] <- "off"
  grid.draw(gt)

library(plyr)

mytitle(ggplot_build(p1), "left")

mytitle <- function(data, text) {
  plot <- data$plot
  panel <- data$panel
  data <- data$data
  theme <- ggplot2:::plot_theme(plot)

  build_grob <- function(layer, layer_data) {
    if (nrow(layer_data) == 0) return()

    dlply(layer_data, "PANEL", function(df) {
      panel_i <- match(df$PANEL[1], panel$layout$PANEL)
      layer$make_grob(df, scales = panel$ranges[[panel_i]], cs = plot$coord)
    }, .drop = FALSE)
  }

    # List by layer, list by panel
  geom_grobs <- Map(build_grob, plot$layer, data)

  plot_table <- ggplot2:::facet_render(plot$facet, panel, plot$coordinates,
                             theme, geom_grobs)

  title <- ggplot2:::element_render(theme, "plot.title", text)
  title_height <- grobHeight(title) +
    if (is.null(plot$labels$title)) unit(0, "lines") else unit(0.5, "lines")

  pans <- plot_table$layout[grepl("^panel", plot_table$layout$name), ,
    drop = FALSE]

  plot_table <- gtable_add_rows(plot_table, title_height, pos = 0)
  plot_table <- gtable_add_grob(plot_table, title, name = "title",
    t = 1, b = 1, l = min(pans$l), r = max(pans$r), clip = "off")

}
