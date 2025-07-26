NPQ_X12_surf = function(depth, fluor, mld,SFM, lthresh){
  # This NPQ correction is new, adjusted from Xing et al. (2012)
  #
  #
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  df$fluor[df$fluor < 0] = 0
  # identify maximum f in upper MLD
  Z = Zeu(df$depth, df$fluor,lthresh)
    depth_max = Z
    NPQ_depth = df$depth[which.max(df$fluor[df$depth < depth_max])]
  
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max(fluor[depth <= NPQ_depth], na.rm = T)
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  if(is.empty(NPQ_depth)){NPQ_depth = 0}
  if(!is.na(NPQ_depth)){if(NPQ_depth == min(df$depth)){NPQ_depth = 0}}
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}

