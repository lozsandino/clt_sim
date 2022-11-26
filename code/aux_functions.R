dSeq <- function(x) {
  res <- seq(min(x), max(x), by = (max(x) - min(x))/1000)
  
  return(res)
}

getMDD <- function(sd, nX, isNorm = TRUE) {
  # Estimates the minimum detectable difference.
  # sd: standard deviation of the sample difference
  # n: size of sample X (nY = nX)
  
  sig <- 0.05 # statistical significance
  pow <- 0.8 # statistical power
  
  if (isNorm) {
    res <- (qnorm(1 - sig/2) + qnorm((1 - pow), lower.tail = FALSE)) * sd * sqrt(2/(nX - 1))
  } else {
    res <- (qt(1 - sig/2, df = nX - 1) + qt((1 - pow), df = nX - 1, lower.tail = FALSE)) * sd * sqrt(2/(nX - 1))
  }
  
  return(res)
}

randomAvgDiff <- function(seed, x, y) {
  # Estimates the average difference of re-sampled x and y.
  
  set.seed(seed)
  
  resX <- mean(sample(x, replace = TRUE))
  resY <- mean(sample(y, replace = TRUE))
  
  res <- resX - resY
  
  return(res)
}

bootsAvgDiff <- function(x, y, iter) {
  # Creates a bootstrap sample of the average difference between x and y.
  
  res <- sapply(1:iter, FUN = randomAvgDiff, x = x, y = y)
  
  return(res)
}

ggTheme <- function() {
  theme(
    plot.title = element_text(hjust = 0.5), 
    text = element_text(size = 14, family = "serif", color = "#3b3b3b"),
    plot.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0"),
    legend.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0")
  )
}


ggHist <- function(x, meanX, densX = NULL, nameX = NULL) {
  
  # Set up
  n <- length(x)
  k <- ceiling(log2(n)) + 1 # Sturges
  avgX <- mean(x)
  
  if (is.null(nameX)) {
    nameX <- substitute(x)
  }
  
  # Plot
  p <- ggplot() +
    geom_histogram(
      aes(x = x, y = after_stat(density)),
      bins = k,
      alpha = 0.8
    ) +
    geom_vline(
      aes(xintercept = avgX, color = "avg"),
      linetype = "dashed", 
      linewidth = 1.2
    ) +
    geom_vline(
      aes(xintercept = meanX, color = "mean"),
      linetype = "dashed", 
      alpha = 0.6,
      linewidth = 1.2
    ) +
    xlab("") +
    scale_color_manual(
      name = "",
      values = c(avg = "#702963", mean = "#CD7F32"), 
      labels = c("sample mean", "true mean")
    ) +
    ggtitle(paste("Histogram of", nameX)) +
    ggTheme()
  
  if (!is.null(densX)) {
    p <- p + 
      geom_line(
        aes(x = dSeq(x), y = densX),
        alpha = 0.8,
        linewidth = 1
      )
  }
  
  return(p)
}

ggConv <- function(x, realParam, var = FALSE) {
  # Plots the convergence of the estimations to the real param. If var = FALSE,
  # then it plots the mean.
  
  # Set up
  n <- length(x)
  
  if (n > 1000) {
    # To decrease the dimension, while maintaining the objective of the plot.
    index <- sample(30:n, 1000, replace = FALSE)
    index <- index[order(index)]
  } else {
    index <- 30:n
  }
  
  est <- c()
  
  
  for (k in index) {
    # randomly select k observations and estimate the average.
    if (var) {
      est <- c(est, var(sample(x, k)))
      type <- "Variance of Difference"
    } else {
      est <- c(est, mean(sample(x, k)))
      type <- "Mean Difference"
    }
    
  }
  
  # Plot
  p <- ggplot() +
    geom_line(
      aes(
        x = 2 * index, # to present the total size of the sample.
        y = est, 
        color = "estim"
      )
    ) +
    geom_hline(
      aes(yintercept = realParam, color = "real"),
      alpha = 0.6,
      linewidth = 1.2
    ) +
    xlab(expression("n")) +
    ylab("value") +
    scale_color_manual(
      name = "",
      # values = c(estim = "#0000FF", real = "#FF0000"), 
      values = c(estim = "#702963", real = "#CD7F32"), 
      labels = c("sample", "true")
    ) +
    ggtitle(paste("Convergence of", type)) +
    ggTheme()
  
  return(p)
}

ggNorm <- function(x) {
  p <- ggplot() +
    geom_qq(
      aes(sample = x),
      alpha = 0.5
    ) +
    stat_qq_line(
      aes(sample = x),
      linewidth = 0.9,
      color = "#000000"
    ) +
    xlab("theoretical") +
    ylab("empirical") +
    ggtitle("Normal QQ Plot") +
    ggTheme()
  
  return(p)
}

ggMDD <- function(x, y, boots) {
  # Plots the densities representing the minimum detectable difference.
  
  # Set up
  sig <- 0.05
  pow <- 0.8
  
  nX <- length(x)
  sdX <- sd(x)
  sdY <- sd(y)
  sd <- sqrt(sdX^2 + sdY^2)
  mdd <- getMDD(sd, nX)
  hvals <- dSeq(c(dSeq(boots) - mdd - 2 * sd * sqrt(2/(nX - 1)), dSeq(boots) + mdd + 2 * sd * sqrt(2/(nX - 1))))
  
  # Densities
  densNull <- dnorm(hvals, sd = sd * sqrt(2/(nX - 1)))
  densAlt <- dnorm(hvals, mean = mdd, sd = sd * sqrt(2/(nX - 1)))
  criticalValue <- qnorm(1 - sig/2, sd = sd * sqrt(2/(nX - 1)))
  
  p <- ggplot() +
    geom_line(
      aes(
        x = hvals,
        y = densAlt,
        color = "alt"
      ),
      linewidth = 1.2
    ) +
    geom_area(
      aes(
        x = ifelse(hvals >= criticalValue, hvals, NA),
        y = densAlt
      ),
      fill = "#702963",
      alpha = 0.3
    ) +
    geom_vline(
      aes(
        xintercept = mdd,
        color = "alt"
      ),
      linetype = "dashed", 
      linewidth = 1.2
    ) +
    geom_line(
      aes(
        x = hvals,
        y = densNull,
        color = "nul"
      ),
      linewidth = 1.2
    ) + 
    geom_area(
      aes(
        x = ifelse(hvals >= criticalValue, hvals, NA),
        y = densNull
      ),
      fill = "#CD7F32",
      alpha = 0.3
    ) +
    geom_area(
      aes(
        x = ifelse(hvals <= -criticalValue, hvals, NA),
        y = densNull
      ),
      fill = "#CD7F32",
      alpha = 0.3
    ) +
    geom_vline(
      aes(
        xintercept = 0,
        color = "nul"
      ),
      linetype = "dashed", 
      linewidth = 1.2
    ) +
    xlab("") +
    xlim(min(hvals), max(hvals)) +
    ylab("density") +
    scale_color_manual(
      name = "",
      values = c(alt = "#702963", nul = "#CD7F32"), 
      labels = c(
        "mdd",
        expression(paste(mu, "= 0"))
      )
    ) +
    ggTheme()
  
  return(p)
}