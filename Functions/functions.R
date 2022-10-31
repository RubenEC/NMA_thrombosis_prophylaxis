#residual deviance plot for detecting outliers (slight modification of function to include names)
mtc.devplot.annotated <- function (x, ...) 
{
  require(tidyverse)
  stopifnot(class(x) == "mtc.deviance")
  if (is.null(x[["dev.re"]])) {
    tpl <- x[["dev.ab"]]
    study <- matrix(rep(1:nrow(tpl), times = ncol(tpl)), 
                    nrow = nrow(tpl), ncol = ncol(tpl))
    study <- t(study)[t(!is.na(tpl))]
    devbar <- t(x[["dev.ab"]])[t(!is.na(tpl))]
    title <- "Per-arm residual deviance"
    xlab <- "Arm"
  }
  else {
    nd <- c(x[["nd.ab"]], x[["nd.re"]])
    devbar <- c(apply(x[["dev.ab"]], 1, sum, na.rm = TRUE), 
                x[["dev.re"]])/nd
    study <- 1:length(devbar)
    title <- "Per-study mean per-datapoint residual deviance"
    xlab <- "Study"
  }
  plot(devbar, ylim = c(0, max(devbar, na.rm = TRUE)), ylab = "Residual deviance", 
       xlab = xlab, main = title, pch = c(1, 22)[(study%%2) + 1], ...)
  for (i in 1:length(devbar)) {
    lines(c(i, i), c(0, devbar[i]))
  }
  names<- subset(as.data.frame(tpl),V1>2|V2>2|V3>2)
  names <-names %>% mutate(
    V1 = ifelse(V1>2, rownames(names),NA),
    V2 = ifelse(V2>2, rownames(names),NA),
    V3 = ifelse(V3>2, rownames(names),NA)) %>%
    pivot_longer(cols = c(V1, V2, V3)) %>%
    filter(!is.na(value))%>%
    pull(value)
  values <- which(devbar >2, arr.ind=TRUE)
  text(x=values, y=devbar[devbar>2], labels = names, cex = 0.5, pos = 3)
}

#leverage plot for detecting outliers (with slight modification of function to include names)
mtc.levplot.annotated <- function (x, ...) 
{
  stopifnot(class(x) == "mtc.deviance")
  fit.ab <- apply(x[["fit.ab"]], 1, sum, na.rm = TRUE)
  dev.ab <- apply(x[["dev.ab"]], 1, sum, na.rm = TRUE)
  lev.ab <- dev.ab - fit.ab
  fit.re <- x[["fit.re"]]
  dev.re <- x[["dev.re"]]
  lev.re <- dev.re - fit.re
  nd <- c(x[["nd.ab"]], x[["nd.re"]])
  w <- sqrt(c(dev.ab, dev.re)/nd)
  lev <- c(lev.ab, lev.re)/nd
  plot(w, lev, xlim = c(0, max(c(w, 2.5))), ylim = c(0, max(c(lev,4))), 
       xlab = "Square root of residual deviance", 
       ylab = "Leverage", main = "Leverage versus residual deviance", 
       ...)
  text(w[w>1.3], lev[w>1.3], labels = names(w[w>1.3]), cex = 0.3, pos = 3)
  mtext("Per-study mean per-datapoint contribution")
  x <- seq(from = 0, to = 3, by = 0.05)
  for (c in 1:4) {
    lines(x, c - x^2)
  }
}

#sucra
sucra = function(x, lower.is.better = FALSE) {
  
  rank.probability = x
  
  
  # Convert rank.probability to matrix
  mat = as.matrix(rank.probability)
  
  # Loop over treatments, for each treatment: calculate SUCRA
  a = ncol(mat)
  j = nrow(mat)
  names = rownames(mat)
  
  sucra = numeric()
  for (x in 1:j) {
    sucra[x] = sum(cumsum(mat[x, 1:(a - 1)]))/(a - 1)
  }
  
  # If condition for lower.is.better
  if (lower.is.better == TRUE) {
    sucra = numeric()
    for (x in 1:j) {
      sucra[x] = 1 - sum(cumsum(mat[x, 1:(a - 1)]))/(a - 1)
    }
  }
  
  # Make data.frame
  res = data.frame(Treatment = names, SUCRA = sucra)
  
  # Order
  res = res[order(-res$SUCRA), ]
  rownames(res) = 1:j
  
  rownames(res) = res$Treatment
  res$Treatment = NULL
  
  class(res) = "sucra"
  
  invisible(res)
  
  res
  
}