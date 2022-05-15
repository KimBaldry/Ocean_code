NPQ_B19 <- function(depth, fluor, bbp){
  # This is a new NPQ correction method based off of the maximum angle method for defining MLD
  # The method works by first calculating tan theta of the fluor:bbp profile 
  # then identifying where the double derivative of tan theta is zero to find a maximum rate of change 
  # where tan theta is also positive
  # This method requires de-spinking or a very smooth profile
  # 
  #
  #
  #
  # 
  
  # calculate eco_mld
  emld = Eco_MLD(depth,fluor)$EMLD
   df = data.frame(depth, fluor, bbp)
   df = df[order(depth),]
   df = df[complete.cases(depth,fluor,bbp),]

   ratio_fb = df$bbp/df$fluor
   max_depth_idx = which.closest(df$depth, emld) + 10
   av_resolution = mean(df$depth[2:max_depth_idx] - df$depth[1:(max_depth_idx-1)], na.rm = T)
   window = 2*round((20/av_resolution)/2) + 1
   window_half = round((20/av_resolution)/2)   
   if(window < 7){window = 7
   window_half = 3}
   ratio_fb_dir =  diff(ratio_fb[1:length(df$fluor)]) < 0
   if(max_depth_idx-window < 10){NPQ_depth = NA
   corr_fluor = fluor}else{
     FLUOR = ratio_fb[1:max_depth_idx]
     DEPTH = df$depth[1:max_depth_idx]
     DEPTH_optic = DEPTH[window_half:(max_depth_idx-window_half)]
     
     tan_theta = vector("list",max_depth_idx-(2*window_half))
     for(sp in window_half:(length(FLUOR)-window_half))
     {
       if(sp < window){
         m1 = lm(FLUOR[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
         m2 = lm(FLUOR[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
         tan_theta[[sp]] = (m2 - m1)/(1+(m2*m1))
       }else{
         m1 = lm(FLUOR[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
         m2 = lm(FLUOR[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
         
         tan_theta[[sp+1 - window_half]] = (m2 - m1)/(1+(m2*m1))}
     }
     tan_theta = unlist(tan_theta)
     ddtan_theta = diff(diff(tan_theta)/diff(depth[1:(length(tan_theta))]))/diff(depth[1:(length(tan_theta))])
 
     corr_fluor = fluor
     if(length(which(ratio_fb_dir[2:max_depth_idx] == T))==0 | length(which(diff(sign(ddtan_theta)) != 0 & tan_theta[3:(length(tan_theta)-3)] > 0 & diff(tan_theta[2:(length(tan_theta)-2)]) < 0)) == 0){
        NPQ_idx = NA
        NPQ_depth = NA
     }else{
        NPQ_idx = which(diff(sign(ddtan_theta)) != 0 & tan_theta[3:(length(tan_theta)-3)] > 0 & diff(tan_theta[2:(length(tan_theta)-2)]) < 0)[1]
        NPQ_depth = DEPTH_optic[NPQ_idx+3]
        
        if(length(which(depth <= NPQ_depth)) < 1){
        NPQ_depth = 0}else{
          fluor_NPQ = (bbp[depth <= NPQ_depth])/ratio_fb[which.closest(depth,NPQ_depth)]
          corr_fluor[depth <= NPQ_depth] = fluor_NPQ # correct above NPQ_depth
        }
        
        }
     
     
   }
   return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  }