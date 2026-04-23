#library("gamlss2")

#d <- bamlss::GAMart()
#d$x4 <- runif(nrow(d))

#b <- gamlss2(num ~ s(x1) + s(x2) + s(x3) + s(x4), data = d)

#ll1 <- family(b)$pdf(d$num, predict(b), log = TRUE)

#nd <- d
#nd$x2 <- nd$x2 + 0.00001

#ll2 <- family(b)$pdf(d$num, predict(b, newdata = nd), log = TRUE)

#dx2 <- (ll1 - ll2) / 0.00001

#x11(); plot(nd$x2, abs(dx2), main = sum(abs(dx2)))
