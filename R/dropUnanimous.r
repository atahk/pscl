## drop unanimous cols of a (rollcall) matrix
dropUnanimous <- function(obj, lop = 0) {
  UseMethod("dropUnanimous")
}

dropUnanimous.rollcall <- function(obj, lop = 0) {
  if (!inherits(obj, "rollcall")) {
    stop("dropUnanimous.rollcall only defined for objects of class rollcall")
  }
  dropRollCall(obj,
    dropList = list(lop = lop)
  )
}

dropUnanimous.matrix <- function(obj, lop = 0) {
  if (!is.matrix(obj)) {
    stop("dropUnanimous.matrix only defined for objects of class matrix")
  }

  if (lop > 1 | lop < 0 | is.na(lop) | !is.numeric(lop) | length(lop) != 1) {
    stop("bad value for lop, must be a single proportion")
  }

  goodObj <- !is.na(obj)
  if (!all(as.vector(obj[goodObj]) %in% c(0, 1, NA))) {
    stop("rollcall matrix contains codes other than 0, 1, and NA.")
  }

  m <- apply(obj, 2, minMargin)
  drop <- m <= lop
  out <- obj[, !drop]
  out
}
