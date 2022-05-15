NPQ_T18 <- function(pres, fluor, bbp, mean_night_f, mean_night_f_bbp){
  #
  # ref:Thomalla, S. J., Moutier, W., Ryan‐Keogh, T. J., Gregor, L., & Schütt, J. (2018). 
  # An optimized method for correcting fluorescence quenching using optical backscattering on autonomous platforms.
  # Limnology and Oceanography: Methods, 16(2), 132-144.  #
  #
  #
  # despiking and smoothing was employed in this paper: a lot of it
  #
  #
  diff = mean_night_f - fluor
  signs = diff/abs(diff)
  idx_five_abs = order[abs(diff)][1:5]
  idx_sign_change = which((signs[-length(signs)]+signs[-1]) == 0) + 1 # it looks like they are grabbing the measurement after the sign change
  NPQ_depth_idx = min(idx_five_abs,idx_sign_change)
  NPQ_depth = pres[NPQ_depth_idx]
  corr_fluor = fluor 
  for(sp in 1:NPQ_depth_idx)
  {
    if(fluor[sp] < bbp[1:NPQ_depth_idx]*mean_night_f_bbp){
        corr_fluor[sp] = bbp[1:NPQ_depth_idx]*mean_night_f_bbp # correct above NPQ_depth
    }
  }
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  rm(diff, signs, idx_sign_change, idx_five_abs, NPQ_depth, NPQ_depth_idx)
}

