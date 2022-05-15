NPQ_X12 = function(depth, fluor, mld){
  # This NPQ correction is from Xing et al. (2012)
  # 
  # ref: Xing, X., Claustre, H., Blain, S., D'Ortenzio, F., Antoine, D., Ras, J., & Guinet, C. (2012). 
  # Quenching correction for in vivo chlorophyll fluorescence acquired by autonomous platforms: 
  # A case study with instrumented elephant seals in the Kerguelen region (Southern Ocean). 
  # Limnology and Oceanography: Methods, 10(7), 483-495.
  #
  #
  #
  # origionally uses a 0.03 critereon for MLD
  #
  #
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  # identify maximum f in upper MLD
  NPQ_depth = df$depth[which.max(df$fluor[df$depth < 0.9*mld])]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max(fluor[depth <= NPQ_depth], na.rm = T)
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  if(is.empty(NPQ_depth)){NPQ_depth = NA}
  if(!is.na(NPQ_depth)){if(NPQ_depth == min(df$depth)){NPQ_depth = NA}}
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}

