simpi <- function(n=1000){
  z <- 0
  res <- .C("simpi",
            PACKAGE=.package.Name,
            as.integer(n),
            as.integer(z))[[2]]
  estimate <- res/n*4
  estimate
}
