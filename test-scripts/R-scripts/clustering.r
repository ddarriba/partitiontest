# Sample from a Chinese Restaurant Process
#
# num_elements: number of samples to draw
# alpha: alpha.  If large then there will tend to be more categories.
#
rchinese <- function(num_elements, alpha) {
  stopifnot(length(alpha) == 1)
  samples = c()
  next_category = 1
  for (i in 1:num_elements) {
    denom = i - 1 + alpha
    p_new = alpha / denom
    rvalue = runif(1)
    ci =
      if (rvalue < p_new) {
        next_category = next_category + 1;
        next_category - 1
      }
      else {
	samples[floor(rvalue * i)]
      }
    samples = c(samples, ci)
  }
  samples
}

rchinese2 <- function(num_elements, alpha) {
  stopifnot(length(alpha) == 1)
  samples = c()
  new_category = 1
  for (i in 1:num_elements) {
    denom = i - 1 + alpha
    p_new = alpha / denom
    ci =
      if (runif(1) < p_new) {
        new_category = new_category + 1;
        new_category - 1
      }
      else {
        counts = unname(table(factor(samples, levels = 1:(i-1))))
        pmf = counts / sum(counts)
        rdiscrete(1, pmf)
      }
    samples = c(samples, ci)
  }
  samples
}

#rchinese(10,0.5)
sample_cluster <- function(max_boxes, num_elements) {
  boxes <- sample(1:max_boxes,num_elements,replace=T)
  return(boxes)
}

cluster <- function(n_boxes, n_elements) {
  if (n_boxes > n_elements) {
    cat("Error in cluster(n_boxes,n_elements\n  n_elements must be >= n_boxes\n")
    return(NA)
  } else {
    boxes <- 1:n_boxes
    remain <- n_elements - n_boxes
    boxes[(n_boxes+1):n_elements] <- sample_cluster(n_boxes, remain)
    return(boxes)
  }
}
