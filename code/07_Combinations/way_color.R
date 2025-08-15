library(akima)
library(grDevices)

x <- c(0, 1)
y <- c(0, 1)

r <- matrix(c(0, 1, 0, 0), nrow = length(x), ncol = length(y))
g <- matrix(c(0, 0, 0, 1), nrow = length(x), ncol = length(y))
b <- matrix(c(0, 0, 0, 0), nrow = length(x), ncol = length(y))

colors <- c("#000000", "#c51b7d", "#000000", "#4d9221")
colors <- col2rgb(colors)
r <- matrix(colors[1, ]/255, nrow = length(x), ncol = length(y))
g <- matrix(colors[2, ]/255, nrow = length(x), ncol = length(y))
b <- matrix(colors[3, ]/255, nrow = length(x), ncol = length(y))


# z <- matrix(c(
#   0, 0, 0,
#   1, 0, 0,
#   0, 1, 0,
#   1, 1, 1
# ))

plotdata <- expand.grid(x = seq(0, 1, length.out = 100), y = seq(0, 1, length.out = 100))

plotdata$r <- bilinear(x, y, r, x0=plotdata$x, y0=plotdata$y)$z
plotdata$g <- bilinear(x, y, g, x0=plotdata$x, y0=plotdata$y)$z
plotdata$b <- bilinear(x, y, b, x0=plotdata$x, y0=plotdata$y)$z
plotdata$rgb <- rgb(plotdata$r, plotdata$b, plotdata$g)

ggplot(plotdata, aes(x, y)) + geom_tile(aes(fill = rgb)) + scale_fill_identity()

fourway_colors <- c("#000000", "#c51b7d", "#4d9221", "#DDDDDD")
threeway_colors <- c("#000000", "#c51b7d", "#4d9221")

map_colors <- function(x0, y0, reverse = FALSE) {
  r <- matrix(c(0, 1, 0, 0), nrow = length(x), ncol = length(y))
  g <- matrix(c(0, 0, 0, 1), nrow = length(x), ncol = length(y))
  b <- matrix(c(0, 0, 0, 0), nrow = length(x), ncol = length(y))
  
  colors <- c(threeway_colors[1], threeway_colors[2], threeway_colors[1], threeway_colors[3])
  if (reverse){
    colors <- c(colors[1], colors[2], colors[4], colors[3])
  }
  colors <- col2rgb(colors)
  r <- matrix(colors[1, ]/255, nrow = length(x), ncol = length(y))
  g <- matrix(colors[3, ]/255, nrow = length(x), ncol = length(y))
  b <- matrix(colors[2, ]/255, nrow = length(x), ncol = length(y))
  
  plotdata <- data.frame(x0=x0, y0=y0)
  
  plotdata$r <- bilinear(x, y, r, x0=plotdata$x, y0=plotdata$y)$z
  plotdata$g <- bilinear(x, y, g, x0=plotdata$x, y0=plotdata$y)$z
  plotdata$b <- bilinear(x, y, b, x0=plotdata$x, y0=plotdata$y)$z
  plotdata$rgb <- rgb(plotdata$r, plotdata$b, plotdata$g)
  
  plotdata$rgb
}

plotdata <- expand.grid(x = seq(0, 1, length.out = 100), y = seq(0, 1, length.out = 100))
plotdata$rgb <- map_colors(plotdata$x, plotdata$y)

ggplot(plotdata, aes(x, y)) + geom_tile(aes(fill = rgb)) + scale_fill_identity()


