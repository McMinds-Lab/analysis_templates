monteCarloP <- function(x, pn='p') {
  if(pn == 'n') {
    res <- (1 + sum(x >= 0, na.rm=TRUE)) / (1 + length(x))
  } else if(pn == 'p') {
    res <- (1 + sum(x <= 0, na.rm=TRUE)) / (1 + length(x))
  }
} # https://arxiv.org/pdf/1603.05766.pdf
