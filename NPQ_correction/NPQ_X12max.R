NPQ_X12maxopt = function(depth, fluor, NPQ_depth){
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  corr_fluor = fluor 
  if(length(which(depth <= NPQ_depth)) < 1){corr_fluor = fluor
  NPQ_depth = 0}else{
    fluor_NPQ = max(fluor[depth <= NPQ_depth], na.rm = T)
    corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
  }
   return(corr_fluor)
  
}

