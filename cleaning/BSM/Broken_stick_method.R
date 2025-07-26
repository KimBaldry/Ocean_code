### Broken stick model
#
#
#


# function to find next break points
next_bp = function(depth, fluor,bps, res){
  # calculate RSS for each depth if a new fit was made
  RSS_vec = lapply(1:length(fluor), FUN = function(x){
    bps_new = rbind(bps,list(depth[x],fluor[x]))
    bps_new = bps_new[order(bps_new[,1])]
    RSS(fluor,approx(bps_new$depth,bps_new$fluor, xout = depth)$y)
  })
  new_bp = list()
  idx_rss = which.min(unlist(RSS_vec))
  idx_min = which(depth >= depth[idx_rss] - res/2)[1]
  idx_max = rev(which(depth <= depth[idx_rss] + res/2))[1]
  new_bp$depth = depth[idx_rss]
  if(length(c(idx_min:idx_max)) <= 2){new_bp$fluor = mean(fluor[idx_min:idx_max],na.rm = T)}else{new_bp$fluor = median(fluor[idx_min:idx_max],na.rm = T)}
  nb = rbind(bps,list(new_bp$depth,new_bp$fluor))
  nb = nb[order(nb[,1]),]
  nb = list("depth" = new_bp$depth, "fluor" = new_bp$fluor, "RSS" = RSS(fluor,approx(nb$depth,nb$fluor, xout = depth)$y), "idx_rss" = idx_rss)

  nb
}

# A broken stick model for the PPZ
# this only computes for required max_depth
# the model choses break point by adding breakpoints in locations that minimises the RSS of the fit,
# only adding a new breakpoint if the additional break point reduces the RSS from the previous fit and
# reduces the coefficient of determination by more than 0.01 (or thresh)
# this means that noise wont cause additional breakpoints to be fit. If the profile is of low resolution, all points should become break-points.
# The user should specify a minimum resolution based on their problem, scale of their study feature. 

# Windowing
  ### get first two break points ###
  # median of first and last 5m - if only one point sets to this one point. if two takes the mean
windowing = function(depth, fluor, thresh = 0.01, resolution_thresh = 5){
  fluor = fluor[which(is.finite(fluor))]
  depth = depth[which(is.finite(fluor))]
  min_depth = min(depth)
  max_depth = max(depth)
  idx1 = which(depth > max_depth-5)
  idx2 = which(depth < min_depth + 5)
  # if(length(idx1) > 2){bp1 = c(median(depth[idx1]),median(fluor[idx1]))}else{bp1 = c(mean(depth[idx1]),mean(fluor[idx1]))}
  # if(length(idx2) > 2){bp2 = c(median(depth[idx2]),median(fluor[idx2]))}else{bp2 = c(mean(depth[idx2]),mean(fluor[idx2]))}
  bp1 = c(mean(depth[idx1], na.rm = T),mean(fluor[idx1], na.rm = T))
  bp2 = c(mean(depth[idx2], na.rm = T),mean(fluor[idx2],na.rm = T))
  # set up break point data frame
  bps = data.table(rbind(bp2,bp1))
  bps = bps[order(bps[,1]),]
  colnames(bps) = c("depth","fluor")
  # rm(bp1, bp2)
  # chop depth and fluor
  # idx = which(depth > bps$depth[1]& depth < bps$depth[2])
  # depth = depth[idx]
  # fluor = fluor[idx]
  # calculate sotal sum of squares
  TSS = tot_var(fluor)
  # calculate first RSS
  RSS_old = RSS(fluor, approx(bps$depth,bps$fluor, xout = depth)$y)
  # coefficient of determination or r2 
  # Note - can have less than 0 as the first fit can be worse than the mean
  CD_old = (TSS-RSS_old)/TSS
  
  # third break point
  new_bp = next_bp(depth, fluor,bps, resolution_thresh)
  # ### if close to ends, re-fit ends
  # if(new_bp$idx_rss %in% idx1){bps = bps[-nrow(bps),]}
  # if(new_bp$idx_rss %in% idx2){bps = bps[-1,]}
  # if(length(idx1) > 2){bp1 = c(median(depth[idx1]),median(fluor[idx1]))}else{bp1 = c(mean(depth[idx1]),mean(fluor[idx1]))}
  # if(length(idx2) > 2){bp2 = c(median(depth[idx2]),median(fluor[idx2]))}else{bp2 = c(mean(depth[idx2]),mean(fluor[idx2]))}
  RSS_new = new_bp$RSS
  CD_new = (TSS-RSS_new)/TSS
  
  # # refit ends
  # RSS_old = RSS_new
  # CD_old = CD_new
  # bps = rbind(bps,list(new_bp$depth,new_bp$fluor))
  # bps = bps[order(bps[,1])]
  # 
  # # fourth bp
  # new_bp = next_bp(depth, fluor,bps, resolution_thresh)
  # RSS_new = new_bp$RSS
  # CD_new = (TSS-RSS_new)/TSS
  
  
  # if RSS is bigger - stop - or if there is not a significant increase in coefficient of determination (5% inc. of var explained)
  if((RSS_old < RSS_new)){print("RSS_old < RSS_new")}
  while( (RSS_old >= RSS_new) & ( (CD_new - CD_old) > thresh) ){
  RSS_old = RSS_new
  CD_old = CD_new
  bps = rbind(bps,list(new_bp$depth,new_bp$fluor))
  bps = bps[order(bps[,1]),]
  # next bp
  new_bp = next_bp(depth, fluor,bps, resolution_thresh)
  # if(new_bp$idx_rss %in% idx1){bps = bps[-nrow(bps),]}
  # if(new_bp$idx_rss %in% idx2){bps = bps[-1,]}
  RSS_new = new_bp$RSS
  CD_new = (TSS-RSS_new)/TSS
  
  }
  
  ### now we have identified the break points for the model
  # lets output this information
  
  result = list()
  bps = rbind(bps,list(new_bp$depth,new_bp$fluor))
  bps = bps[order(bps[,1]),]
  result$n = nrow(bps)
  ## this extends to min and max depths - makes broken stick look nice
  #mod1 = lm(fluor ~ depth , data = bps[1:2])
  #mod2 = lm(fluor ~ depth , data = bps[(nrow(bps)-1):nrow(bps)])
  #bps = rbind(bps,list(min_depth,mod1$coefficients[1] + mod1$coefficients[2]*min_depth))
  #bps = rbind(bps,list(max_depth,mod2$coefficients[1] + mod2$coefficients[2]*max_depth))
  bps = bps[order(bps[,1]),]
  result$break_points = bps
  result$RSS = RSS_new
  result$CD = CD_new
  result$n_wind = nrow(bps) - 1
  result$L = bps$depth[1:(nrow(bps)-1)]
  result$U = bps$depth[2:nrow(bps)]
  return(result)
}


# fit_bp_segmodel_2 = function(depth, fluor,thresh = 0.01){
#   #CHL50 = CHL_50(depth[is.finite(fluor)],fluor[is.finite(fluor)])
#   # EMLD = Eco_MLD(depth[is.finite(fluor)],prof_data$CTDFLUOR[is.finite(fluor)])
#   # if(EMLD$QI > 0.3){
#   #   depth = depth[which(depth < (EMLD$EMLD + 50))]
#   #   fluor = fluor[which(depth < (EMLD$EMLD + 50))]
#   # }
#  
#   sp = 1
#   
#   fit_lm = lm(fluor~depth)
#   mod = segmented.lm(fit_lm, seg.Z = ~depth, control = seg.control(alpha = 0.0001, n.boot = 50, h = 0.5), npsi = sp)
#   BIC_new = BIC(mod)
#   TSS = tot_var(fluor)
#   # calculate first RSS
#   RSS_new = RSS(fluor, mod$fitted.values)
#   # coefficient of determination or r2 
#   # Note - can have less than 0 as the first fit can be worse than the mean
#   CD_new = (TSS-RSS_new)/TSS
#   BIC_old = BIC_new
#   CD_old = 0
#   RSS_old = RSS_new
#   br = F
#   
#   while( (RSS_old >= RSS_new) & ( (CD_new - CD_old) > thresh) & (BIC_old >= BIC_new) & !br){
#     
#     sp= sp+1
#     tryCatch({mod = segmented.lm(fit_lm, seg.Z = ~depth, control = seg.control(alpha = 0.0001, n.boot = 50, h = 0.5), npsi = sp)}, warning=function(w){ 
#       br = T}, error = function(e){br = T})
#     BIC_old = BIC_new
#     CD_old = CD_new
#     RSS_old = RSS_new
#     
#     BIC_new = BIC(mod)
#     RSS_new = RSS(fluor, mod$fitted.values)
#     CD_new = (TSS-RSS_new)/TSS
#     
#   }
#   final_sp = sp - 1
#   #Two chances for better fit
#   br = F
#   chance = 1
#   while(chance <=2){
#     sp = sp+1
#     tryCatch({mod1 = segmented.lm(fit_lm, seg.Z = ~depth, control = seg.control(alpha = 0.0001, n.boot = 50, h = 0.5), npsi = sp)}, warning=function(w){
#       br = T}, error = function(e){br = T})
#     BIC_new = BIC(mod1)
#     RSS_new = RSS(fluor, mod1$fitted.values)
#     CD_new = (TSS-RSS_new)/TSS
#     
#     if((RSS_old >= RSS_new) & ( (CD_new - CD_old) > thresh) & (BIC_old >= BIC_new) & !br){
#       while( (RSS_old >= RSS_new) & ( (CD_new - CD_old) > thresh) & (BIC_old >= BIC_new) & !br){
#         sp= sp+1
#         tryCatch({mod = segmented.lm(fit_lm, seg.Z = ~depth, control = seg.control(alpha = 0.0001, n.boot = 50, h = 0.5), npsi = sp)}, warning=function(w){
#           br = T}, error = function(e){br = T})
#         
#         BIC_old = BIC_new
#         CD_old = CD_new
#         RSS_old = RSS_new
#         
#         BIC_new = BIC(mod)
#         RSS_new = RSS(fluor, mod$fitted.values)
#         CD_new = (TSS-RSS_new)/TSS
#       }
#       final_sp = sp-1
#       
#     }else{chance = chance + 1
#     br = F}
#   }
#   
#   
#   mod = segmented.lm(fit_lm, seg.Z = ~depth, control = seg.control(alpha = 0.0001, n.boot = 50, h = 0.5), npsi = final_sp)
#   return(list("mod" = mod, "BIC" = BIC_new, "RSS" = RSS_new, "depth.in" = depth, "fluor.in" = fluor))    
# }



alpha = 4
# local fitting
# calculate the regression for each window
# need to decide if we want an overlap or not - I dont think so as there ae sharp increases and we dont want to contaminate segments "too much"
# choose an alpha that corresponds to d = 6 sd

## adjust so only influenced by neighboring segments
# triDiag(1,1,1,W)
PBSM_fit = function(result, depth, fluor){
  fluor = fluor[which(is.finite(fluor))]
  depth = depth[which(is.finite(fluor))]
  
  U = result$U # upper break-point
  L = result$L # lower break-point
  dt = data.table(d = depth, f = fluor)
  fits = list() # linear fits in each stick
  mps = list() # mid-points
  sds = list()
  W = length(U)
  for(w in 1:W){
    if(w == 1){
      L[w] = min(dt$d)
      fits[[w]] = lm(f ~ d, data = dt[d <= U[w]])
    mps[[w]] = ( L[w] + U[w] ) /2
    sds [[w]] = (  U[w] - L[w] ) /alpha
    }
    if(w == W){
      U[w] = max(dt$d)
      fits[[w]] = lm(f ~ d, data = dt[ d >= L[w]])
      mps[[w]] = (U[w] + L[w] ) /2
      sds [[w]] = (U[w] - L[w] ) /alpha
      }
    if(w > 1 & w < W){
      fits[[w]] = lm(f ~ d, data = dt[d <= U[w] & d >= L[w]])
      mps[[w]] = ( U[w] + L[w] ) /2
      sds [[w]] = ( U[w] - L[w] ) /alpha
      }
  }
  
  # adjust break points based on intersection of fits and claculate new sd and mps and U/L 
  
  # because we have low resolution data we are going to create T not based on our repeted measurements
  # instead we are going to make a 0.5 m profile
  T = seq(ceiling(min(depth)), floor(max(depth)), 0.1)
  
  # predict values and rates of change for all local fits over T
  mean_fits = matrix(NA,nrow = length(T), ncol = W) # mean fit
  int_fits = matrix(NA,nrow = length(T), ncol = W)
  rate_fits = matrix(NA,nrow = length(T), ncol = W) # slope fits
  wind_dist = matrix(NA,nrow = length(T), ncol = W) # weights middle values higher in the fit
  norm_wind_dist = matrix(NA,nrow = length(T), ncol = W)
  wind_dist_mask = matrix(0,nrow = length(T), ncol = W)
  for(w in 1:W){
    int = predict(fits[[w]], data.frame(d = T), interval = "prediction")
    mean_fits[,w] = int[,1]
    # because we have sparse data dont worry about CI?
    #int_fits[[w]] = int[,1] - int[,2]
    rate_fits[,w] = fits[[w]]$coefficients[2]
    
    # window influence
    # using a flat prior
    # and modelling distribution of regresstion measurements
    # b_w
    # repeated measures distribution for local fits
    #
    wind_dist[,w] = dnorm(T, mps[[w]], sds[[w]]) 
    idx = c(max(1,w-1):min(W,w+1))
    wind_dist_mask[T <= U[w] & T >= L[w],idx] = 1
    rm(idx)
    
    # b_w[t]
  }
  wind_dist_masked = wind_dist * wind_dist_mask
  for(w in 1:W){
  # posterior distribution for fits 
  norm_wind_dist[,w] = wind_dist_masked[,w] / rowSums(wind_dist_masked)
  }
  # set up matricies
  
  # calculate final fit and U/L and rate of change
  pbsm = rowSums(norm_wind_dist * mean_fits * wind_dist_mask)
  pbsm_deriv = rowSums(norm_wind_dist * rate_fits * wind_dist_mask)
  
  return(list(depth = T,fluor_fit = pbsm, rate_fit = pbsm_deriv))
}


# Now we fit a probabalistic broken-stick model to produce a smoothed fit with contiuous rate of change


breaks = bps$break_points$depth
U = breaks[1]
L = breaks[length(breaks)]
# fit linear regressions -> thetas

#test
pf = "33LG20060101_CTD_576_1"
depth = prof_data$DEPTH
fluor = prof_data$CTDFLUOR
max_depth = 200
plot( prof_data$CTDFLUOR,prof_data$DEPTH, col = "green",ylim = c(100,0), type = "l",xaxt = "n", ylab = "Depth", xlab = "")
lines( bps$fluor,bps$depth,col = "black",pch = 4)
lines(pbsm, T, col  = "red")
for(w in 1:W){
  lines(mean_fits[,w],T, col = "blue")
}


# get the part of teh profile that we are interested in and remove clear local outliers
# order profile by pressure and cut at max_depth
depth = depth[order(depth)]
fluor = fluor[order(depth)]
idx = depth < max_depth
depth = depth[idx]
fluor = fluor[idx]
rm(idx)
# remove spikes/outliers that may sway the analysis (+/- 3sd around a 5pt MA mean)
# may need to consider not doing this for low resolution profiles?
fluor = rm_out(fluor)
idx = complete.cases(depth,fluor)
depth = depth[idx]
fluor = fluor[idx]
o_depth = depth
o_fluor = fluor
rm(idx)


