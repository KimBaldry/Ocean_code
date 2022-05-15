library(caTools)
## remove outlying data points
## remove outlying data points
rm_out <- function(y){
  y <- as.numeric(y)
  ny <- length(y)
  
  if(ny <= 5 & ny > 3){m <- mean(y, na.rm = T)
  s <- sd(y, na.rm = T)
  }
  
  if(ny > 5){
    m = runmean(y, k = 5, endrule = "constant")
    s = runsd(y, k = 5, endrule = "constant")
      }
  if(ny > 3){
    out_idx = which(y < (m - 3*s) | y > (m + 3*s))
    y[out_idx] = NA}
  y
}


mean_5m <- function(depths,fx,depth_out,window = 5){
  ddx = which(depths < depth_out + window/2 & depths > depth_out - window/2)
  if(length(ddx) == 0){m = NA}else{
    sub_x = fx[ddx]
  if(length(which(is.finite(sub_x))) < 5){m = mean(sub_x,na.rm = T)}else{
    m <- mean(sub_x, na.rm = T)
    s <- sd(sub_x, na.rm = T)
    sub_x[which(sub_x < (m - 3*s) | sub_x > (m + 3*s))] = NA
    m = mean(sub_x,na.rm = T)
  if(is.nan(m)){m = NA}}}
  m
}

