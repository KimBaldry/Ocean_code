NPQ_X12maxopt = function(depth, fluor, NPQ_depth){
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = df$fluor[which.closest(df$depth,NPQ_depth)]
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
   return(corr_fluor)
  
}

NPQ_X12maxfun = function(depth, fluor, temp, dens, MLD, SFM, TML, FWL){
  
  if(min(depth,na.rm = T) > 10){
    DENS_10 = dens[!is.na(dens)][1]
    T_10 = temp[!is.na(temp)][1]}else{DENS_10 = approx(x = depth, y = dens, xout = 10)$y
    T_10 = approx(x = depth, y = temp, xout = 10)$y}
  if(MLD > SFM){NPQ_depth = depth[which.max(fluor[dens < (DENS_10 + 0.010)])]}# NPQ MR2 only - 0.016
  if(MLD < SFM){NPQ_depth = depth[which.max(fluor[dens < (DENS_10 + 0.076)])]}# NPQ MR1 only - 0.11
  
  if(!is.na(FWL)){ 
  if(FWL){NPQ_depth = depth[which.max(fluor[dens < (DENS_10 + 0.008)])]}} #MR 4
  if(!is.na(TML)){
  if(TML){NPQ_depth = depth[which.max(fluor[temp < (DENS_10 + 0.21)])]}} #MR3
  
  
  
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = df$fluor[which.closest(df$depth,NPQ_depth)]
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}

NPQ_X12max_light = function(depth, fluor, temp, dens, MLD, SFM, TML, FWL){
  
  df = data.frame(depth,fluor)
  df = df[complete.cases(df),]
  kd = 0.0232+0.074*(fluor^0.674)
  ddiff = depth - c(0,df$depth[1:(length(fluor)-1)])
  kddiff = kd - c(0,kd[1:(length(kd)-1)])
  integral_tmp = ddiff*(kd + kddiff/2) # area under curve using linear approximation
  integral = cumsum(integral_tmp)
  ltmp = exp(-integral)
  
  if(min(depth,na.rm = T) > 10){
    DENS_10 = dens[!is.na(dens)][1]
    T_10 = temp[!is.na(temp)][1]}else{DENS_10 = approx(x = depth, y = dens, xout = 10)$y
    T_10 = approx(x = depth, y = temp, xout = 10)$y}
  
  if(MLD > SFM){NPQ_depth = depth[which.max(fluor[ltmp > 0.16])]}# 0.065 for NPQ only
  if(MLD < SFM){NPQ_depth = depth[which.max(fluor[ltmp > 0.084])]} # 0.037 for NPQ only
  
  if(!is.na(FWL)){
    if(FWL){NPQ_depth = depth[which.max(fluor[temp < (T_10 + 0.004)])]}}
  if(!is.na(TML)){
    if(TML){NPQ_depth = depth[which.max(fluor[ltmp > 0.083])]}}

  
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = df$fluor[which.closest(df$depth,NPQ_depth)]
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}
