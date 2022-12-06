NPQ_S08_surf = function(depth, fluor, bbp, lthresh){
  # This NPQ correction is new, adjusted from Sackmann et al. (2008)
  #
  #
  df = data.frame(depth, fluor,bbp)
  df = df[complete.cases(df),]
  df$fluor[df$fluor < 0] = 0
  r = df$fluor/df$bbp # ratio used by Sackmann et al. (2008)
  Z = Zeu(df$depth, df$fluor,lthresh)
  # identify maximum f in upper MLD
    depth_max = min(Z,df$depth[which.max(df$fluor)], na.rm = T)
    max_r = max(r[df$depth < depth_max], na.rm = T) # identify maximum r in upper MLD
    NPQ_depth = df$depth[which(r == max_r)[1]]

  
  corr_fluor = fluor 
  if(length(which(df$depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  if(is.empty(NPQ_depth)){NPQ_depth = 0}
  if(!is.na(NPQ_depth)){if(NPQ_depth == min(df$depth)){NPQ_depth = 0}}
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}


