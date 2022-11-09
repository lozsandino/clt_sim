dSeq <- function(x) {
  res <- seq(min(x), max(x), by = (max(x) - min(x))/1000)
  
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
    theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 14, family = "serif", color = "#3b3b3b"),
      plot.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0"),
      legend.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0")
    )
  
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
  n <- min(1000, length(x))
  est <- c()
  
  for (k in 30:n) {
    # randomly select k observations and estimate the average.
    if (var) {
      est <- c(est, var(sample(x, k)))
      type <- "Variance"
    } else {
      est <- c(est, mean(sample(x, k)))
      type <- "Mean Difference"
    }
    
  }
  
  # Plot
  p <- ggplot() +
    geom_line(
      aes(x = 30:n, y = est, color = "estim")
    ) +
    geom_hline(
      aes(yintercept = realParam, color = "real"),
      alpha = 0.6,
      linewidth = 1.2
    ) +
    xlab("n") +
    ylab("mean difference") +
    scale_color_manual(
      name = "",
      # values = c(estim = "#0000FF", real = "#FF0000"), 
      values = c(estim = "#702963", real = "#CD7F32"), 
      labels = c("sample", "true")
    ) +
    ggtitle(paste("Convergence of", type)) +
    theme(
      plot.title = element_text(hjust = 0.5), 
      text = element_text(size = 14, family = "serif", color = "#3b3b3b"),
      plot.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0"),
      legend.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0")
    )
  
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
    theme(
      plot.title = element_text(hjust = 0.5), 
      text = element_text(size = 14, family = "serif", color = "#3b3b3b"),
      plot.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0"),
      legend.background = element_rect(fill = "#e0e0e0", color = "#e0e0e0")
    )
  
  return(p)
}
