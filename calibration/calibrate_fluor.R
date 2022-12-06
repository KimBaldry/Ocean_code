

calibrate_fluor = function(depth, fluor, chl, chl_depth){
  cor = unlist(lapply(1:length(chl), FUN = function(z){mean_5m(depth,fluor, depth_out = chl_depth[z], window = 2)}))
  if(is.na(cor[1])){cor[1] = fluor[which(is.finite(fluor))[1]]}
  sc_data = data.frame(y = chl, fy = cor, depth = chl_depth)[which(chl > 0.2*max(chl, na.rm = T) & cor>0),]
  sc = optim(par = 1, fn = function(par){sum((sc_data$y/(max(sc_data$y, na.rm = T))^2)*(sc_data$y-par*sc_data$fy)^2,na.rm = T)},method =  "L-BFGS-B", lower = 0, upper = Inf)
  fluor*sc$par
}
