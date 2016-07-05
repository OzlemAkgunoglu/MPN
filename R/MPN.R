#' Most Probable Number (MPN) calculator
#'
#' @param g the number of positive (or growth) tubes in the jth dilution
#' @param m the amount of the original sample put in the tubes at each dilution.
#' @param m the number of tubes at each dilution
#' @return The most probable  \code{x} and \code{y}.
#' @examples
#' MPN(g=c(5,3,0),m=c(0.1,0.01,0.001),t=10)

MPN <- function(g,m,t) {
  upper <- sum(g)/sum((t-g)*m)
  min_fun <- function(lambda){
    abs(sum(g*m/(1-exp(-lambda*m)))-sum(t*m))}
  solution <- optimize(f=min_fun,
        lower=0,
        upper=upper) #setting this too high can lead it to stick to these high values
  return(solution$minimum)
}


#Add upper and lower bound calculations

# MPN <- function(g,m,t) {
#   min_fun <- function(lambda){
#     abs(sum(g*m/(1-exp(-lambda*m)))-sum(t*m))}
#   solution <- optim(fn=min_fun,
#       par=c(lambda=1),
#       method="Brent",
#       lower=0,
#       upper=100)
# return(solution)
# }



