
# A probabalistic broken-stick model to produce a smoothed fit with contiuous rate of change


alpha = 4
# local fitting
# calculate the regression for each window
# need to decide if we want an overlap or not - I dont think so as there ae sharp increases and we dont want to contaminate segments "too much"
# choose an alpha that corresponds to d = 6 sd

## adjust so only influenced by neighboring segments
# triDiag(1,1,1,W)
PBSM_fit = function(result, depth, fluor, alpha = 4){
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

