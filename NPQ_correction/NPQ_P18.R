NPQ_P18 <- function(depth, fluor, mld){
  #for this correction fluor must be in chla or a scale factor (sf) provided
  #
  # ref: Xiaogang Xing, Nathan Briggs, Emmanuel Boss, and Hervé Claustre, 
  # "Improved correction for non-photochemical quenching of in situ chlorophyll fluorescence based on a synchronous irradiance profile," 
  # Opt. Express 26, 24734-24751 (2018)
  #
  # This calculates zeu from [Chla], but we dont know the sensitivity to this.. and the [Chla] will be less with quenching... 
  # I suppose that this means that the the ez will just be deeper..
  # maybe a sensitivity test should be done
  # The test will be different for BGC-Argo 
  #
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  # kd = 0.0232+0.074*(df$fluor^0.674)
  # ddiff = df$depth - c(0,df$depth[1:(length(df$fluor)-1)])
  # kddiff = kd - c(0,kd[1:(length(kd)-1)])
  # integral_tmp = ddiff*(kd + kddiff/2) # area under curve using linear approximation
  # integral = cumsum(integral_tmp)
  # tmp = exp(-integral)
  # zeu = approx(tmp,df$depth, 0.01)$y
  zeu = Zeu(df$depth, df$fluor)
  NPQ_depth = df$depth[which.max(df$fluor[df$depth <= min(zeu,mld)])]
  corr_fluor = fluor
  if(length(which(depth <= NPQ_depth)) < 1){
  NPQ_depth = 0}else{
  fluor_zp18 = max(fluor[depth <= NPQ_depth], na.rm = T)
  corr_fluor[depth <= NPQ_depth] = fluor_zp18 # correct above NPQ_depth
  }
  if(is.empty(NPQ_depth)){NPQ_depth = 0}
  if(!is.na(NPQ_depth)){if(NPQ_depth == min(df$depth)){NPQ_depth = 0}}
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))

}
