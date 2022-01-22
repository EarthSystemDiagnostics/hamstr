#' Extract Stan Parameter Index
#'
#' @param x vector of stan parameter names 
#'
#' @return
#' @keywords internal
get_par_idx <- function(x) {
  
  a <- strsplit(x, "\\[")
  a <- sapply(a, function(x) x[2])
  b <- unlist(strsplit(a, "\\]"))
  
  as.numeric( b )
}