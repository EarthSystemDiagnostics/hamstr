

# #debugonce(GetIndices)
# 
# debugonce(GetBrksHalfOffset)
# 
#brks <- GetBrksHalfOffset(K_fine = c(99), K_factor = 2)
#brks
#brks <- rev(brks)

#GetIndices(brks = brks)
#GetIndices(nK = c(2,3,5))


# #debugonce(GetIndices)
# nK <- c(1, 3,5)
# GetIndices(nK)



#hierarchical_depths2(fit_HP2$data)

#' #' Get Nearest Primes to Geometric Series
#' #'
#' #' @param m Terminal number
#' #' @param f Approximate ratio of series
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' GetNextMultiPrimes(100, 3)
#' GetNextMultiPrimes <- function(m, f = 2){
#'   x <- 1
#'   np <- 1
#'   
#'   while(np < m){
#'     np <- as.numeric(gmp::nextprime(f*tail(x, 1)-1))
#'     x <- append(x, np)
#'   }
#'   
#'   x <- x[x <= m]
#'   
#'   if (max(x) < 0.5*m){
#'    x <- append(x, m)
#'   }else{
#'     x <- append(head(x, -1), m)
#'   }
#'   
#'   return(x)
#' }





