## Install required packages
pkgs <- c(
  "gamlss.dist", "ggplot2", "hexSticker",
  "showtext", "dplyr"
)

for(p in pkgs) {
  if(!p %in% installed.packages()) {
    install.packages(p)
  }
}

## Load required packages
library("gamlss.dist")
library("ggplot2")
library("hexSticker")
library("showtext")
library("dplyr")
library("gamlss2")

## Add a Google font to use in the logo text
font_add_google(name = "Source Code Pro", family = "sourcecode")

## Enable automatic use of showtext for better font rendering
showtext_auto()

## Generate some data
#set.seed(1)
#x <- seq(-3, 3, length = 2000)
#y <- sin(x) + rnorm(length(x), sd = exp(-1.5 + cos(x)))
#plot(x, y)

## Load data
data("mcycle", package = "MASS")
x <- mcycle$times
y <- mcycle$accel

## Fit model
b <- gamlss2(y ~ s(x, k = 20) | s(x, k = 20))

## Predict quantiles
qmat <- quantile(b, prob = seq(0.01, 0.99, length = 501))

## Prepare color palette
k <- (ncol(qmat) - 1) / 2
col <- hcl.colors(floor(1.05 * k), "Blue-Yellow", rev = TRUE, alpha = 0.05)[1:k]

## Prepare basic plot
p <- ggplot() + theme_void()

## Plot ribbons manually from outer to inner
center <- (ncol(qmat) + 1) / 2
for(j in rev(0:(k - 1))) {
  lower_idx <- center - j
  upper_idx <- center + j
  ribbon_df <- data.frame(
    x = x,
    ymin = qmat[, lower_idx],
    ymax = qmat[, upper_idx]
  )
  
  p <- p + geom_ribbon(
    data = ribbon_df,
    aes(x = x, ymin = ymin, ymax = ymax),
    fill = col[j + 1],
    inherit.aes = FALSE
  )
}

center_line_df <- data.frame(x = seq(min(x), max(x), length = 300))
par <- predict(b, newdata = center_line_df)
center_line_df$y <- family(b)$mean(par)
p <- p + geom_line(
  data = center_line_df,
  aes(x = x, y = y),
  color = "white",
  linewidth = 1, ## Adjust thickness as you like (maybe 1.2-1.5 looks best)
  inherit.aes = FALSE
)

## Create the hex sticker with customized layout and fonts
sargs <- list(
  subplot = p,
  package = "gamlss2",
  p_size = 60,
  p_fontface = "bold",
  p_color = "white",
  s_width = 2.2,
  s_height = 1,
  s_x = 1,
  s_y = 1.36,
  p_x = 1,
  p_y = 0.8,
  h_color = hcl(180, 60, 75),
  h_fill = hcl(250, 60, 35),
  h_size = 1.5,
  l_alpha = 0.4,
  dpi = 600,
  spotlight = TRUE
)

## logo version with spotlight
set.seed(2)
s_logo <- do.call("sticker", sargs)

## sticker version with spotlight and URL
sargs <- c(sargs, list(
  url = "https://gamlss-dev.github.io/gamlss2/",
  u_size = 8.52,
  u_color = "white",
  u_x = 0.99,
  u_y = 0.08,
  filename = "gamlss2_sticker.png"
))
set.seed(2)
s_sticker <- do.call("sticker", sargs)
