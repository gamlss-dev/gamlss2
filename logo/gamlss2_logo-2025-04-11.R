if(!"gamlss.dist" %in% installed.packages()) {
  install.packages("gamlss.dist")
}
if(!"berryFunctions" %in% installed.packages()) {
  install.packages("berryFunctions")
}
if(!"showtext" %in% installed.packages()) {
  install.packages("showtext")
}

library("gamlss.dist")
library("berryFunctions")
library("showtext")

# Load a Google Font
##font_add_google(name = "Special Elite", family = "typewriter")
##font_add_google(name = "IBM Plex Sans", family = "logo")
##font_add_google(name = "Orbitron", family = "orbitron")
##font_add_google(name = "Exo 2", family = "exo2")
##font_add_google(name = "Courier Prime", family = "texttt")
font_add_google(name = "Source Code Pro", family = "sourcecode")

# Activate showtext
showtext_auto()

graphics.off()

png <- FALSE

if(png) {
  png("gamlss2.png", units = "in",
    width = 5, height = 5, res = 500,
    bg = "transparent")
} else {
  x11(width = 5, height = 5)
}

par(mar = c(0, 0, 0, 0))#, xaxs = "i", yaxs = "i")

y <- seq(0, 10, length = 300)
dy <- dBCPE(y, mu = 2, sigma = 3, nu = 3.5, tau = 3)
dy0 <- dBCPE(y, mu = 3.5, sigma = 2.5, nu = 3, tau = 4)
dy1 <- rev(dBCPE(y, mu = 2, sigma = 2.5, nu = 3, tau = 4))

q1 <- qBCPE(0.1, mu = 5, sigma = 0.7, nu = 3, tau = 3)
q2 <- qBCPE(0.2, mu = 5, sigma = 0.7, nu = 3, tau = 3)
q3 <- qBCPE(0.8, mu = 5, sigma = 0.7, nu = 3, tau = 3)
q4 <- qBCPE(0.9, mu = 5, sigma = 0.7, nu = 3, tau = 3)

gen_pol <- function(x, mu = 3.5, sigma = 2.5, nu = 3, tau = 4) {
  x <- seq(0, x, length = 300)
  rbind(
    cbind(x, -0.02),
    cbind(rev(x), rev(dBCPE(x, mu = mu, sigma = sigma, nu = nu, tau = tau)))
  )
}

xlim <- range(y)
ylim <- range(dy)

ylim <- ylim + c(-0.64, 0.15) * abs(diff(ylim))

plot(NA, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")

rx <- abs(diff(xlim))*0.015
ry <- abs(diff(ylim))*0.015

roundedRect(xlim[1] + rx, ylim[1] + ry, xlim[2] - rx, ylim[2] - ry,
  rounding = 0.2, lwd = 20, col = "#CBAF7A", border = "#1B1B1B")

i <- 1:length(y)
i <- i[10:(length(i) - 10)]

lines(dy[i] ~ y[i], lwd = 20, col = "#1B1B1B")

polygon(gen_pol(10, mu = 2, sigma = 3, nu = 3.5, tau = 3),
  col = rgb(0.1, 0.1, 0.1, alpha = 0.2), border = NA)

polygon(gen_pol(10), col = "#888C8D", border = NA)
polygon(gen_pol(q4), col = "#5A8CA5", border = NA)
polygon(gen_pol(q3), col = "#D08C2A", border = NA)
polygon(gen_pol(q2), col = "#5A8CA5", border = NA)
polygon(gen_pol(q1), col = "#888C8D", border = NA)

lines(dy0[i] ~ y[i], lwd = 20, col = "#A33D2D")
#lines(dy1[i] ~ y[i], lwd = 20, col = rgb(0.1, 0.1, 0.1, alpha = 0.2))

## lines(rep(0, length(y))[i] ~ y[i], lwd = 20, col = "#A33D2D")

roundedRect(xlim[1] + rx, ylim[1] + ry, xlim[2] - rx, ylim[2] - ry,
  rounding = 0.2, lwd = 20, border = "#1B1B1B")

text(5, -0.13, "gamlss2", font = 2, cex = 6 , col = "#1B1B1B")

if(png)
  dev.off()

