## these functions are convenient to have in ~/.Rprofile for interactive use. I also put logit/inv_logit from 'link_functions' in  my .Rprofile

# head, but subset columns as well as rows
headc <- function(x,nr=5,nc=nr) { head(x,nr)[,1:nc] }

# i always forget the conversion between fastq Q-scores and error rate
q2p <- \(q) { 10^(-0.1 * q) }
p2q <- \(p) { -10 * log10(p) }
