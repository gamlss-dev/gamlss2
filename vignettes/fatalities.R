## https://www.euromomo.eu/index.html
tdir <- tempfile()
dir.create(tdir)
file_url <- "https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/demo_r_mwk3_t.tsv.gz"
download.file(file_url, file.path(tdir, "fatalities.tsv.gz"))

## Background information:
## https://www.statcan.gc.ca/eng/statistical-programs/document/3233_D5_V1

d <- read.csv(gzfile(file.path(tdir, "fatalities.tsv.gz")), sep = "\t")

countries <- d$unit.geo.time
countries <- gsub("NR,", "", countries)

md <- list()
ld <- list()
for(j in unique(countries)) {
  dd <- unlist(d[countries == j, ][-1])
  md[[j]] <- as.numeric(dd)
  names(md[[j]]) <- names(dd)
  ld[[j]] <- length(na.omit(md[[j]]))
}

ld <- unlist(ld)
barplot(ld)

fd <- list()
for(j in "AT") {
  td <- subset(d, countries == j)[, -1]
  td <- unlist(td)
  n <- names(td)
  td <- gsub("p", "", td)
  td <- gsub(" ", "", td)
  td <- as.numeric(td)
  names(td) <- n
  td <- na.omit(td)
  n <- names(td)
  td <- as.numeric(td)
  n <- gsub("X", "", n)
  n <- strsplit(n, "W")
  year <- sapply(n, function(x) as.integer(x[1]))
  week <- sapply(n, function(x) as.integer(x[2]))

  td <- data.frame("num" = td, "year" = year, "week" = week)
  td <- td[order(td$week), ]
  td <- td[order(td$year), ]
  rownames(td) <- NULL

  fd[[j]] <- td
}

fatalities <- data.frame("num" = fd[[1]]$num, "year" = fd[[1]]$year, "week" = fd[[1]]$week)

plot(num ~ week, data = fatalities)
d2 <- subset(fatalities, year >= 2020)
points(d2$week, d2$num, pch=16,col=2)

library("gamlss2")

f <- list(num~s(week,bs="cc"),~s(week,bs="cc"),~s(week,bs="cc"),~s(week,bs="cc"))

b <- gamlss2(f, data = fatalities, family = BCPE)

nd <- data.frame("week" = 1:53)

par <- predict(b, newdata = nd, type = "parameter")

nd$q95 <- family(b)$q(0.95, par)
nd$q50 <- family(b)$q(0.5, par)
nd$q05 <- family(b)$q(0.05, par)

plot(num ~ week, data = fatalities, pch = 16,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.4))
with(nd, matplot(week, cbind(q95, q50, q05),
  type = "l", lty = c(2, 1, 2), lwd = 3,
  col = 4, add = TRUE))
with(d2, points(week, num, pch = 16, col = 2))

