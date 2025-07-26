
NPQ_BCPmaxfun = function(depth, fluor, bbp, temp, dens, MLD, SFM, TML, FWL){
  
  if(min(depth,na.rm = T) > 10){
    DENS_10 = dens[!is.na(dens)][1]
    T_10 = temp[!is.na(temp)][1]}else{DENS_10 = approx(x = depth, y = dens, xout = 10)$y
    T_10 = approx(x = depth, y = temp, xout = 10)$y}
  
  r = fluor/bbp
  
  if(MLD > SFM){
    cond = dens < (DENS_10 + 0.004)
    max_r = which.max(r[cond]) 
    NPQ_depth = depth[max_r]
    max_r = r[max_r]
    }
  if(MLD < SFM){    
    cond = dens < (DENS_10 + 0.060)
    max_r = which.max(r[cond]) 
    NPQ_depth = depth[max_r]
    max_r = r[max_r]
  }
  
  if(!is.na(FWL)){
    if(FWL){
      cond = dens < (DENS_10 + 0.002)
      max_r = which.max(r[cond]) 
      NPQ_depth = depth[max_r]
      max_r = r[max_r]
    }}
  if(!is.na(TML)){
    if(TML){
      cond = dens < (DENS_10 + 0.24)
      max_r = which.max(r[cond]) 
      NPQ_depth = depth[max_r]
      max_r = r[max_r]
    }}
  
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(df$depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}

NPQ_BCPmax_light = function(depth, fluor, bbp, temp, dens, MLD, SFM, TML, FWL){
  
  df = data.frame(depth,fluor,bbp)
  df = df[complete.cases(df),]
  kd = 0.0232+0.074*(fluor^0.674)
  ddiff = df$depth - c(0,depth[1:(length(fluor)-1)])
  kddiff = kd - c(0,kd[1:(length(kd)-1)])
  integral_tmp = ddiff*(kd + kddiff/2) # area under curve using linear approximation
  integral = cumsum(integral_tmp)
  ltmp = exp(-integral)
  if(min(depth,na.rm = T) > 10){
    DENS_10 = dens[!is.na(dens)][1]
    T_10 = temp[!is.na(temp)][1]}else{DENS_10 = approx(x = depth, y = dens, xout = 10)$y
    T_10 = approx(x = depth, y = temp, xout = 10)$y}
  
  r = fluor/bbp
  
  if(MLD > SFM){
    cond = ltmp > 0.14
    max_r = which.max(r[cond]) 
    NPQ_depth = depth[max_r]
    max_r = r[max_r]
  }
  if(MLD < SFM){    
    cond = ltmp > 0.11
    max_r = which.max(r[cond]) 
    NPQ_depth = depth[max_r]
    max_r = r[max_r]
  }
  
  if(!is.na(FWL)){
    if(FWL){
      cond = temp < (T_10 + 0.004)
      max_r = which.max(r[cond]) 
      NPQ_depth = depth[max_r]
      max_r = r[max_r]
    }}
  if(!is.na(TML)){
    if(TML){
      cond = ltmp > 0.083
      max_r = which.max(r[cond]) 
      NPQ_depth = depth[max_r]
      max_r = r[max_r]
    }}  
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(df$depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}
