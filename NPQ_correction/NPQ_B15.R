NPQ_B15 <- function(depth, fluor, ed){
  # This NPQ correction is from Biermann et al. (2015)
  #
  # ref: Biermann, L., Guinet, C., Bester, M., Brierley, A. S., & Boehme, L. (2015). 
  # An alternative method for correcting fluorescence quenching. 
  # Ocean Science.
  #
  #
  # 
  #
  #
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  max_f = max(df$fluor[df$depth < ed], na.rm = T) # identify maximum f in upper MLD
  NPQ_depth_idx = which(df$fluor[df$depth < ed] == max_f)
  NPQ_depth_idx = NPQ_depth_idx[length(NPQ_depth_idx)] 
  NPQ_depth = df$depth[NPQ_depth_idx] # depth at which maximum f ocurrs
  corr_fluor = fluor 
  corr_fluor[depth < NPQ_depth] = max_f # correct above NPQ_depth
  if(is.empty(NPQ_depth)){NPQ_depth = NA}
  if(!is.na(NPQ_depth)){if(NPQ_depth == min(df$depth)){NPQ_depth = NA}}
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
}
