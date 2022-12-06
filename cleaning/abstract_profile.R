abstract_fluor = function(depth, fluor){
chl_d = CHL_20(depth, fluor) + 75
if(length(which(depth >= chl_d)) > 30){
  fluor_mod = fit_bp_segmented(depth[which(depth < chl_d)], fluor[which(depth < chl_d)], CD_thresh = 0)
  deep = fluor[which(depth >= chl_d & depth < 500)]
  if( all( ((deep - mean(deep, na.rm = T))/diff(range(fluor, na.rm = T) )) < 0.05, na.rm = T)){
    fluor_mod_down = list("fluor.out" = rep(mean(deep, na.rm = T),length(deep)),"bp" = NULL)
  }else{
    fluor_mod_down = fit_bp_segmented(depth[which(depth >= chl_d & depth < 500)], fluor[which(depth >= chl_d & depth < 500)], CD_thresh = 0)
  }
  fluor_mod$fluor.out = c(fluor_mod$fluor.out, fluor_mod_down$fluor.out)
  fluor_mod$bp = c(fluor_mod$bp, chl_d, fluor_mod_down$bp)}else{
    fluor_mod = fit_bp_segmented(depth[which(depth < 500)], fluor[which(depth < 500)], CD_thresh = 0)
  }

f_off = optic_dark(depth, fluor,  fluor_mod$bp[length(fluor_mod$bp)])
fluor_adj = fluor_mod$fluor.out
fluor_adj[which(fluor_adj > mean(fluor_adj[depth < chl_d], na.rm = T)+6*sd(fluor_adj[depth < chl_d], na.rm = T))] = NA
fluor_adj = (fluor_adj - f_off)
return(c(fluor_adj,(fluor[which(depth >= 500)] - f_off)))

}

abstract_bbp = function(depth,bbp,fluor){
  d = depth[is.finite(bbp)]
  if(mean(d[-length(d)] - d[-1], na.rm = T) < 3){
    bbp_adj = runmed(bbp,5, na.action = "na.omit")}else{bbp_adj = bbp}
  
  chl_d = CHL_20(depth, fluor) + 75
  if(length(which(depth >= chl_d)) > 30){
    fluor_mod = fit_bp_segmented(depth[which(depth < chl_d)], bbp_adj[which(depth < chl_d)], CD_thresh = 0)
    deep = bbp_adj[which(depth >= chl_d & depth < 500)]
    if( all( ((deep - mean(deep, na.rm = T))/diff(range(bbp_adj, na.rm = T) )) < 0.05, na.rm = T)){
      fluor_mod_down = list("fluor.out" = rep(mean(deep, na.rm = T),length(deep)),"bp" = NULL)
    }else{
      fluor_mod_down = fit_bp_segmented(depth[which(depth >= chl_d & depth < 500)], bbp_adj[which(depth >= chl_d & depth < 500)], CD_thresh = 0)
    }
    fluor_mod$fluor.out = c(fluor_mod$fluor.out, fluor_mod_down$fluor.out)
    fluor_mod$bp = c(fluor_mod$bp, chl_d, fluor_mod_down$bp)}else{
      fluor_mod = fit_bp_segmented(depth[which(depth < 500)], bbp_adj[which(depth < 500)], CD_thresh = 0)
    }
  
  bbp_off = optic_dark(depth, bbp_adj, fluor_mod$bp[length(fluor_mod$bp)])  
  bbp_adj = fluor_mod$fluor.out - bbp_off
  return(c(bbp_adj,bbp[which(depth >= 500)] - bbp_off))
}