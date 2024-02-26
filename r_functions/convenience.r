## these functions are convenient to have in ~/.Rprofile for interactive use. I also put logit/inv_logit from 'link_functions' in  my .Rprofile

# head, but subset columns as well as rows
headc <- function(x,nr=5,nc=nr) { head(x,nr)[,1:nc] }
