# Sample from a Chinese Restaurant Process
#
# size: number of samples to draw
# a: alpha.  If large then there will tend to be more categories.
#
rchinese <- function(size, a) {
  stopifnot(length(a) == 1)
  samples = c()
  new.category = 1
  for (i in 1:size) {
    denom = i - 1 + a
    p.new = a / denom
    ci =
      if (runif(1) < p.new) {
        new.category = new.category + 1;
        new.category - 1
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
cluster <- function(num_elements, max_boxes) {
  boxes <- sample(1:num_elements,max_boxes,replace=T)
  return(boxes)
}
