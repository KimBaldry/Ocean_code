sbm_adjust = function(depth,fluor,bbp){
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  sfm = df$depth[which.max(df$fluor)]
  
  local_m = localMaxima(df$bbp)
  peak_vals = df$bbp[local_m]
  
  # local bbp/bcp maxima
  local_sbm = which.closest(df$depth[local_m],sfm)
  sbm = length(local_sbm) > 0
  if(sbm){
    sh = sfm - df$depth[local_m[local_sbm]]
  }else{sh = NA}
  if(abs(sh) > 0 & abs(sh) <= 5){
  depth_bbp = df$depth + sh
  new_bbp = approx(x = depth_bbp,y = df$bbp, xout = depth)$y}else{new_bbp = bbp}
  return(list(bbp = new_bbp,shift = sh,SBM_test = sbm))
}
