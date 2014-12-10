library(ggplot2)

pp <- function (n,r=4) {
 x <- seq(-r*pi, r*pi, len=n)
 df <- expand.grid(x=x, y=x)
 df$r <- sqrt(df$x^2 + df$y^2)
 df$z <- cos(df$r^2)*exp(-df$r/6)
 df
}
        p <- ggplot(pp(20)[sample(20*20, size=200),], aes(x=x, y=y))
        p <- p + geom_tile(aes(fill=z))
        p <- p +  theme(panel.grid = element_blank())
        p <- p +  theme(strip.background = element_blank())
        print(p)

