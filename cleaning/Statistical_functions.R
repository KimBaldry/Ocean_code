# RSS, total variance and remove outlier functions


# function to calculate the residual sum of squares (RSS)
RSS = function(y,fy){sum((y-fy)^2,na.rm = T)}

# function to calculate total variance
tot_var = function(y){sum((y-mean(y,na.rm = T))^2, na.rm = T)}


# flag fluorescence profiles with high variability
var_test = function(fluor){
  if(var(fluor,na.rm = T) > 0.2 * diff(range(fluor,na.rm = T))){FLAG = "High"}else{FLAG = "LOW"}
  return(FLAG)
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}