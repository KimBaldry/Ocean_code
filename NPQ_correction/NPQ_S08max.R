NPQ_S08maxopt  <- function(depth, fluor, bbp, NPQ_depth, max_r){
  df = data.frame(depth, fluor, bbp)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max_r*(bbp[depth <= NPQ_depth])
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
  
  return(corr_fluor)
  rm(NPQ_depth,NPQ_depth_idx,corr_fluor,r,max_r)
}