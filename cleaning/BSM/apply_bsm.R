apply_bsm = function(depth, fluor){
  chl_d = CHL_20(depth, fluor) + 75
  if(length(which(depth >= chl_d)) > 30){
    fluor_mod = fit_bp_segmented(depth[which(depth < chl_d)], fluor[which(depth < chl_d)], CD_thresh = 0)
    deep = fluor[which(depth >= chl_d & depth< 500)]
    if( all( ((deep - mean(deep, na.rm = T))/diff(range(fluor, na.rm = T) )) < 0.05, na.rm = T)){
      fluor_mod_down = list("fluor.out" = rep(mean(deep, na.rm = T),length(deep)),"bp" = NULL)
    }else{
      fluor_mod_down = fit_bp_segmented(depth[which(depth >= chl_d & depth < 500)], fluor[which(depth >= chl_d & depth < 500)], CD_thresh = 0)
    }
    fluor_mod$fluor.out = c(fluor_mod$fluor.out, fluor_mod_down$fluor.out)
    fluor_mod$bp = c(fluor_mod$bp, chl_d, fluor_mod_down$bp)}else{
      fluor_mod = fit_bp_segmented(depth[which(depth < 500)], fluor[which(depth < 500)], CD_thresh = 0)
    }
  fluor_mod
}