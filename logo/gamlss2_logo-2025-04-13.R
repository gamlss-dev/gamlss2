## Install required packages.
pkgs <- c(
  "gamlss.dist", "ggplot2", "hexSticker",
  "showtext", "dplyr", "colorspace"
)

for(p in pkgs) {
  if(!p %in% installed.packages()) {
    install.packages(p)
  }
}

## Load required packages.
library("gamlss.dist")
library("ggplot2")
library("hexSticker")
library("showtext")
library("dplyr")
library("colorspace")

## Add a Google font to use in the logo text.
font_add_google(name = "Source Code Pro", family = "sourcecode")

## Enable automatic use of showtext for better font rendering.
showtext_auto()

## Close any open graphics devices.
graphics.off()

## Generate the y-axis values.
y <- seq(0, 10.5, length = 300)

## Calculate density curves using the BCPE distribution.
dy  <- dBCPE(y, mu = 2, sigma = 3, nu = 3.5, tau = 3)
dy0 <- dBCPE(y, mu = 3.5, sigma = 2.5, nu = 3, tau = 4)
dy1 <- rev(dBCPE(y, mu = 2, sigma = 2.5, nu = 3, tau = 4))  ## Not used in plot.

## Combine data into a tibble.
df <- tibble(
  y = y,
  dy = dy,
  dy0 = dy0,
  dy1 = dy1
)

## Convert data to long format for plotting two curves.
df_long <- df %>%
  select(y, dy0, dy) %>%
  tidyr::pivot_longer(cols = -y, names_to = "curve", values_to = "density")

## Create long format for only the white-overlaid curve.
df_long2 <- df %>%
  select(y, dy0) %>%
  tidyr::pivot_longer(cols = -y, names_to = "curve", values_to = "density")

## Compute quantiles from the BCPE distribution.
q <- list()
for(i in seq(0.0001, 0.99999, length = 9)) {
  q[[paste0(i * 100, "%")]] <- qBCPE(i, mu = 5, sigma = 0.7, nu = 3, tau = 3)
}

## Helper function to generate polygons under the BCPE curve.
gen_pol <- function(x, mu = 3.5, sigma = 2.5, nu = 3, tau = 4) {
  xx <- seq(0, x, length = 300)
  tibble(
    x = c(xx, rev(xx)),
    y = c(rep(-0.013  , 300), rev(dBCPE(xx, mu, sigma, nu, tau)))
  )
}

## Generate a color palette for the quantile polygons.
col <- rainbow_hcl(length(q))

## Create layered polygons from outer to inner quantiles.
polys <- lapply(rev(2:length(q)), function(i) {
  gen_pol(q[[i]]) %>%
    mutate(fill = col[i], alpha = 1)
})

## Add a base shadow-like polygon behind all others.
polys <- c(
  list(gen_pol(10, mu = 2, sigma = 3, nu = 3.5, tau = 3) %>%
         mutate(fill = rgb(0.1, 0.1, 0.1), alpha = 0.4)),
  polys
)

## Combine all polygons into one data frame.
poly_df <- dplyr::bind_rows(polys, .id = "group")

## Create the plot using ggplot2.
p <- ggplot() +
  geom_polygon(data = poly_df,
               aes(x = x, y = y, group = group, fill = fill, alpha = alpha),
               color = NA) +
  geom_line(data = df_long2,
            aes(x = y, y = density, color = curve),
            linewidth = 1.3) +
  scale_color_manual(values = c(dy = "#1B1B1B", dy0 = "white")) +
  scale_fill_identity() +
  scale_alpha_identity() +
  theme_void() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(min(poly_df$y), max(df$dy) * 1.1))

## Create the hex sticker with customized layout and fonts.
s <- sticker(
  subplot = p,
  package = "gamlss2",
  p_size = 60,
  p_fontface = "bold",
  s_width = 1.85,
  s_height = 0.9,
  s_x = 1,
  s_y = 1.4,
  p_x = 1,
  p_y = 0.75,
  h_color = "#1B1B1B",
  h_size = 1.5,
  spotlight = TRUE,
  l_alpha = 0.4,
  url = "https://gamlss-dev.github.io/gamlss2/",
  u_size = 8.52,
  u_color = "white",
  u_x = 0.99,
  u_y = 0.08,
  dpi = 600
)

## Display the sticker in the plotting window.
plot(s)

