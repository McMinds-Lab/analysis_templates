logit <- function(x) { log(x / (1 - x)) }
inv_logit <- function(x) { 1 / (1 + exp(-x)) }

