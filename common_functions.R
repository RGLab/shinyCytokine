library(testthat)
## functions common to one or more of 'preprocess.R', 'ui.R', 'server.R'

## get combinations up to an order n
## combine is a function that tells us how
## we combine a particular set of combos
## note that order of elements is 'preserved'
combinations <- function(x, combine=c, ...) {
  combine <- match.fun(combine)
  output <- vector("list", length(x))
  for( i in 1:length(x) ) {
    tmp <- combn(x, i, simplify=FALSE)
    output[[i]] <- lapply(tmp, combine, ...)
  }
  return( unlist(output, recursive=FALSE) )
}

expect_identical( 
  unlist( combinations(1:6)[1:6] ),
  1:6
)
