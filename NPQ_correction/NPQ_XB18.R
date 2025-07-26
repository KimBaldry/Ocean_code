NPQ_XB18 <- function(pres, fluor, PAR){
  #
  # Xiaogang Xing, Nathan Briggs, Emmanuel Boss, and Hervé Claustre, 
  # "Improved correction for non-photochemical quenching of in situ chlorophyll fluorescence based on a synchronous irradiance profile," 
  # Opt. Express 26, 24734-24751 (2018)
  #
  #
  #
  #
  #
  corr_fluor = fluor
  corr_fluor[pres < 10] = fluor[pres < 10]/(0.092 + 0.908/(1+(PAR[which.closest(pres,10)[length(which.closest(pres,10))]]/261)^2.2))
  corr_fluor[pres >= 10] = fluor[pres >= 10]/(0.092 + 0.908/(1+(PAR[pres >= 10]/261)^2.2)) # correct above NPQ_depth
  NPQ_depth = max(pres[(corr_fluor-fluor)/fluor > 0.1]) # retrospectively calculate NPQ_depth 
  return(list("NPQdepth" = NPQ_depth, "corr_fluor" = corr_fluor))
  rm(corr_fluor, NPQ_depth)
}