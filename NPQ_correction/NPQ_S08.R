NPQ_S08  <- function(depth, fluor, bbp,mld){
  # this function calculates the NPQ_depth and a NPQ corrected fluoresence profile based on Sackmann et al. (2008) 
  # ref:  Sackmann, B., Perry, M., and Eriksen, C.: 
  #       Seaglider observations of variability in daytime fluorescence quenching of chlorophyll-a in 
  #       Northeastern Pacific coastal waters, 
  #       Biogeosciences Discussions, 5, 2839-2865, 2008.
  #
  #
  # this paper used MLD = 0.125 kg/m^3
  # this method originally includes a visual classification into
  # 1. quenched above MLD
  # 2. quenched below MLD, but uniform ratio exissts below
  # 3. quenched below MLD but no uniform region of ratio exists -> these profiles cannot be corrected
  #
  #
  # Note that units for NPQ depth are inherrited from pres. So, m for Depth and dbar for pressure measurements
  #
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  r = df$fluor/df$bbp # ratio used by Sackmann et al. (2008)
  max_r = max(r[df$depth < mld], na.rm = T) # identify maximum r in upper MLD
  NPQ_depth = df$depth[which(r == max_r)[1]]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  rm(NPQ_depth,NPQ_depth_idx,corr_fluor,r,max_r)
}