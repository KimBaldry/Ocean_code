### Broken stick model
# Author: Kimberlee Baldry
#
#
# This simple broken stick model
# this only computes for required max_depth
# the model choses break point by adding breakpoints in locations that minimises the RSS of the fit,
# only adding a new breakpoint if the additional break point reduces the RSS from the previous fit and
# reduces the coefficient of determination by more than 0.01 (or thresh)
# this means that noise wont cause additional breakpoints to be fit. If the profile is of low resolution, all points should become break-points.
# The user should specify a minimum resolution based on their problem, scale of their study feature. 


# function to find next break points
next_bp = function(depth, fluor,bps){
  # calculate RSS for each depth if a new fit was made.
  RSS_vec = lapply(1:length(fluor), FUN = function(x){
    bps_new = rbind(bps,list(depth[x],fluor[x]))
    bps_new = bps_new[order(depth)]
    RSS(fluor,approx(bps_new$depth,bps_new$fluor, xout = depth)$y)
  })
  new_bp = list()
  # find the point that results in the minimum RSS
  idx = which.min(unlist(RSS_vec))
  # return depth, fluorescence and the new RSS of the fit
  new_bp$depth = depth[idx]
  new_bp$fluor = fluor[idx]
  new_bp$RSS = unlist(RSS_vec)[idx]
  new_bp
}


# Windowing
### get first two break points ###
# median of first and last 5m - if only one point sets to this one point. if two takes the mean
windowing = function(depth, fluor, thresh = 0.01){
  fluor = fluor[which(is.finite(fluor))]
  depth = depth[which(is.finite(fluor))]
  # get the median of the first and last 5m. This will be the first and last point of the model. Note, if only two points takes the mean.
  min_depth = min(depth)
  max_depth = max(depth)
  idx1 = which(depth > max_depth-5)
  idx2 = which(depth < min_depth + 5)
  if(length(idx1 > 2)){bp1 = c(median(depth[idx1]),median(fluor[idx1]))}else{bp1 = c(mean(depth[idx1]),mean(fluor[idx1]))}
  if(length(idx2 > 2)){bp2 = c(median(depth[idx2]),median(fluor[idx2]))}else{bp2 = c(mean(depth[idx2]),mean(fluor[idx2]))}
  # set up break point data frame to record results
  bps = data.table(rbind(bp2,bp1))
  colnames(bps) = c("depth","fluor")
  # rm(bp1, bp2)
  # chop depth and fluor to within first and last break points
  idx = which(depth > bps$depth[1]& depth < bps$depth[2])
  depth = depth[idx]
  fluor = fluor[idx]
  
  ### fit 3rd break point using next_bp()
  # calculate sotal sum of squares
  TSS = tot_var(fluor)
  # calculate first RSS
  RSS_old = RSS(fluor, approx(bps$depth,bps$fluor, xout = depth)$y)
  # coefficient of determination or r2 
  # Note - can have less than 0 as the first fit can be worse than the mean
  CD_old = (TSS-RSS_old)/TSS
  
  # third break point
  new_bp = next_bp(depth, fluor,bps)
  RSS_new = new_bp$RSS
  CD_new = (TSS-RSS_new)/TSS
  
  ### Fit all break points
  # Check break point meets requirements. If it does continue to fit break points until it doesnt
  # if RSS is bigger - stop - or if there is not a significant increase in coefficient of determination (5% inc. of var explained)
  while( (RSS_old > RSS_new) & ( (CD_new - CD_old) > thresh) ){
  RSS_old = RSS_new
  CD_old = CD_new
  bps = rbind(bps,list(new_bp$depth,new_bp$fluor))
  bps = bps[order(depth)]
  # next bp
  new_bp = next_bp(depth, fluor,bps)
  RSS_new = new_bp$RSS
  CD_new = (TSS-RSS_new)/TSS
  
  }
  # now we have identified the break points for the model
  
  ### return results
  result = list()
  bps = rbind(bps,list(new_bp$depth,new_bp$fluor))
  bps = bps[order(depth)]
  # number of break points
  result$n = nrow(bps)
  ## this extends to min and max depths - makes broken stick look nice
  #mod1 = lm(fluor ~ depth , data = bps[1:2])
  #mod2 = lm(fluor ~ depth , data = bps[(nrow(bps)-1):nrow(bps)])
  #bps = rbind(bps,list(min_depth,mod1$coefficients[1] + mod1$coefficients[2]*min_depth))
  #bps = rbind(bps,list(max_depth,mod2$coefficients[1] + mod2$coefficients[2]*max_depth))
  bps = bps[order(depth)]
  # break points data frame
  result$break_points = bps
  # RSS of final fit
  result$RSS = RSS_new
  # Coefficient of determination of final fit
  result$CD = CD_new
  # number of sticks/windows
  result$n_wind = nrow(bps) - 1
  # Lower list of break point depths
  result$L = bps$depth[1:(nrow(bps)-1)]
  # Upper list of breakpoint depths
  result$U = bps$depth[2:nrow(bps)]
  return(result)
}
