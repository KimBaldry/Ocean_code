### This calculates different MLD

## MLD
# use linear interpolation for now but look into using a smoothing fit
# Assuptions
MLD <- function(press, sal, temp, lat,dens_thresh)
{
  # flag 1 = good
  # flag 2 = MLD cant be calculated due to not enough water column sampled
  # flag 3 = reference density taken from deeper than 10 m

  # create ctd data frame
  ctd_data = data.frame(cbind(press, sal, temp))
  colnames(ctd_data) = c("press","sal", "temp")
  ctd_data = ctd_data[complete.cases(ctd_data),]

  # calculate density
  DENS = as.numeric(rho(S = ctd_data$sal, T = ctd_data$temp, P = 0))

  # remove spikes or outliers (3 standard dev away from neighboring 5 measurements)
  DENS_s = rm_out(DENS)

  # calculate depth
  DEPTH_gws = swDepth(ctd_data$press, latitude = lat)

  # remove missing values
  dens_prof = data.frame("depth" = DEPTH_gws,"density" = DENS_s)
  dens_prof = dens_prof[complete.cases(dens_prof),]


  max_d = max(dens_prof$depth, na.rm = T)
  min_d = min(dens_prof$depth, na.rm = T)
if(max_d <= 10){mld = NA
flag = 2}else{
  
  if(min_d > 10){
  DENS_10 = dens_prof$dens[which(dens_prof$depth == min_d)]
  if(max(dens_prof$density,na.rm = T) < DENS_10 + dens_thresh){flag = 2
  mld = NA}else{
    if(length(dens_prof$dens) >4){
     mld = approx(x = dens_prof$dens, y = dens_prof$depth, xout = DENS_10 + dens_thresh)$y
      flag = 3 
    }else{flag = 2}
      }

  }else{
    if(length(dens_prof$dens) > 4){
    DENS_10 = approx(x = dens_prof$depth, y = dens_prof$dens, xout = 10)$y
    if(max(dens_prof$density,na.rm = T) < DENS_10 + dens_thresh | length(which(dens_prof$depth >= 10)) < 3){flag = 2
    mld = NA}else{
      mld = approx(x = dens_prof$dens[which(dens_prof$depth >= 10)], y = dens_prof$depth[which(dens_prof$depth >= 10)], xout = DENS_10 + dens_thresh)$y
      flag = 1}
    }else{flag = 2
  mld = NA}}}
  out = list()
  out$MLD = mld
  out$FLAG = flag
  out
}

## calculate fit on increasing intervals
sigmoid_fit <- function(depth, val, mld, max_d){9

  dat = data.frame(x = depth, y = val)
  mod = optim(par = c(0,1,0), fn = SS_res,method =  "BFGS", data = dat, lower = c(0,0,0), upper = c(max_d, Inf, max_d- (1.5*mld)))

  r2 = 1 - (SS_res(dat, par, surf_val)/SS_tot(data))

  z_0.5 = mod$par[1]
  s = mod$par[2]
  z_max = par[3]
  val_res = with(data[1:(mld*1.5 + par[3]),], y - surf_val/(1 + (s*exp(z_0.5 - x))))

}

SS_res <- function(data, par, surf_val, mld) {
  with(data[1:(mld*1.5 + par[3]),], sum((surf_val/(1 + (par[2]*exp(par[1] - x))) - y)^2))
}

SS_tot <- function(data){
  with(data, sum((y-mean(data$y))^2))
}
Eco_MLD <- function(depth, fluor){
  depth = depth[order(depth)]
  fluor = fluor[order(depth)]
  cdx = complete.cases(depth, fluor)
  depth = depth[cdx]
  fluor = fluor[cdx]
  min_depth_idx = max(which.max(fluor))
  max_depth_idx = length(depth)
  av_resolution = mean(depth[2:max_depth_idx] - depth[1:(max_depth_idx-1)], na.rm = T)
  window = 2*round((20/av_resolution)/2) + 1
  window_half = round((20/av_resolution)/2)
  if(window < 7){window = 7
  window_half = 3}
  if(max_depth_idx-window < 10){Bloom_depth = NA
  QI_bd = NA}else{
    FLUOR = fluor[1:max_depth_idx]
    DEPTH = depth[1:max_depth_idx]
    DEPTH_optic = DEPTH[window_half:(max_depth_idx-window_half)]
    tan_theta = vector("list",max_depth_idx-(2*window_half))
    for(sp in window_half:(length(FLUOR)-window_half))
    {
      if(sp < window){
        #if(sd(FLUOR[sp:(sp-window_half)]))
        m1 = lm(FLUOR[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
        m2 = lm(FLUOR[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
        tan_theta[[sp]] = (m2 - m1)/(1+(m2*m1))
      }else{
        m1 = lm(FLUOR[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
        m2 = lm(FLUOR[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
        
        tan_theta[[sp+1 - window_half]] = (m2 - m1)/(1+(m2*m1))}
    }
    tan_theta = unlist(tan_theta)
    
    local_m = localMaxima(tan_theta)
    sign = diff(FLUOR,lag = 1)[local_m+window_half]
    local_m = local_m[which(sign < 0)]
    local_m = local_m[which((local_m+window_half) >= min_depth_idx)]
    if(length(local_m) == 0){ Bloom_depth = NA
    QI_bd = NA}else{
      peak_vals = tan_theta[local_m]
      bloom_idx = local_m[which.max(peak_vals)[1]]
      Bloom_depth = DEPTH_optic[bloom_idx]
      
      idx_u = which.closest(DEPTH_optic, 1.5*Bloom_depth)
      if(idx_u > max_depth_idx){idx_u = max_depth_idx}
      idx_l = 1
      # if(idx_l < 1){idx_l =1}
      if(idx_l == bloom_idx){QI_bd = 0}else{
        QI_bd = 1 - sd(FLUOR[bloom_idx:idx_u]-mean(FLUOR[bloom_idx:idx_u]))/
          sd(FLUOR[idx_l:idx_u]-mean(FLUOR[idx_l:idx_u]))}
    }
  }
  return(list("EMLD" = Bloom_depth, "QI" = QI_bd))
}



CHL_50 <- function(pres, fluor){
  pres = pres[order(pres)]
  fluor = fluor[order(pres)]
  fluor = rm_out(fluor)
  val = (range(fluor,na.rm = T)[2]-range(fluor,na.rm = T)[1])/2 + range(fluor,na.rm = T)[1]
  pres[max(which(fluor > val & pres < 300),na.rm = T)]
}

CHL_20 <- function(pres, fluor){
  pres = pres[order(pres)]
  fluor = fluor[order(pres)]
  val = (range(fluor,na.rm = T)[2]-range(fluor,na.rm = T)[1])*0.2 + range(fluor,na.rm = T)[1]
  pres[max(which(fluor > val & pres < 300),na.rm = T)]
}

Zeu = function(depth,fluor,l = 0.01){
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  df$fluor[df$fluor < 0] = 0
  kd = 0.0232+0.074*(df$fluor^0.674)
  ddiff = df$depth - c(0,df$depth[1:(length(df$fluor)-1)])
  kddiff = kd - c(0,kd[1:(length(kd)-1)])
  integral_tmp = ddiff*(kd + kddiff/2) # area under curve using linear approximation
  integral = cumsum(integral_tmp)
  tmp = exp(-integral)
  
  zeu = approx(tmp,df$depth, l)$y
  if(is.na(zeu)){zeu = min(df$depth)}
  zeu
  
}


MLD_mam <- function(depth, dens, max_d){
  depth = depth[order(depth)]
  dens = dens[order(depth)]
  dx = complete.cases(depth, dens)
  depth = depth[dx]
  dens = dens[dx]
  max_depth_idx = max(which(depth < max_d), na.rm = T)
  av_resolution = mean(depth[2:max_depth_idx] - depth[1:(max_depth_idx-1)], na.rm = T)
  window = 2*round((20/av_resolution)/2) + 1
  window_half = round((20/av_resolution)/2)
  if(window < 7){window = 7
  window_half = 3}
  if(max_depth_idx-window < 10){Bloom_depth = NA
  QI_bd = NA}else{
    DENS = dens[1:max_depth_idx]
    DEPTH = depth[1:max_depth_idx]
    DEPTH_optic = DEPTH[window_half:(max_depth_idx-window_half)]
    tan_theta = vector("list",max_depth_idx-(2*window_half))
    for(sp in window_half:(length(DENS)-window_half))
    {
      if(sp < window){
        #if(sd(DENS[sp:(sp-window_half)]))
        m1 = lm(DENS[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
        m2 = lm(DENS[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
        tan_theta[[sp]] = (m2 - m1)/(1+(m2*m1))
      }else{
        m1 = lm(DENS[sp:(sp-window_half)]~DEPTH[sp:(sp-window_half)])$coefficients[2]
        m2 = lm(DENS[(sp):(sp+window_half)]~DEPTH[(sp):(sp+window_half)])$coefficients[2]
        
        tan_theta[[sp+1 - window_half]] = (m2 - m1)/(1+(m2*m1))}
    }
    tan_theta = unlist(tan_theta)
    ## first local max
    local_m = localMaxima(tan_theta)
    peak_vals = tan_theta[local_m]
    bloom_idx = local_m[which.max(peak_vals)[1]]
    
    ## second local max
    peak_vals = peak_vals[-(which(local_m == bloom_idx))]
    local_m = local_m[-(which(local_m == bloom_idx))]
    bloom_idx2 = local_m[which.max(peak_vals)[1]]
    # depths
    Bloom_depth = DEPTH_optic[bloom_idx]
    Bloom_depth2 = DEPTH_optic[bloom_idx2]
    
    # qualities
    idx_u = which.closest(DEPTH_optic, 1.5*Bloom_depth)
    QI_bd1 = 1 - sd(DENS[1:bloom_idx]-mean(DENS[1:bloom_idx]))/
        sd(DENS[1:idx_u]-mean(DENS[1:idx_u]))
    if(bloom_idx > bloom_idx2){
      idx_u = ifelse(bloom_idx > 6, bloom_idx - 5, bloom_idx)
      idx_l = 1}else{
      idx_l = bloom_idx + 5
      idx_u = which.closest(DEPTH_optic, 1.5*Bloom_depth2)}
    QI_bd2 = 1 - sd(DENS[idx_l:bloom_idx2]-mean(DENS[idx_l:bloom_idx2]))/
      sd(DENS[idx_l:idx_u]-mean(DENS[idx_l:idx_u]))
  }
  return(list("mam1" = Bloom_depth, "QI1" = QI_bd1, "mam2" = Bloom_depth2, "QI2" = QI_bd2))
}

complex_fluor = function(depth, fluor, emld){
  if(is.na(emld)){emld = 500}
  df = data.frame(depth, fluor)
  df = df[complete.cases(df),]
  df = df[df$depth < emld,]
  
  local_m = localMaxima(df$fluor)
  peak_vals = df$fluor[local_m]
  local_m = local_m[order(peak_vals, decreasing = T)]
  peak_vals = peak_vals[order(peak_vals, decreasing = T)]
  local_m = local_m[which(local_m != 1)]
  
  if(length(local_m) > 1){  
  for(sp in 2:length(local_m)){
    idx = sp
    if(abs(df$depth[local_m[1]] - df$depth[local_m[sp]]) > 20){sp = length(local_m)}
  }
  
  

  if(peak_vals[idx] > 1.2*mean(df$fluor[1:3])){
    r = TRUE
  }else{
    r = FALSE
  }}else{r = FALSE}
  r
}
