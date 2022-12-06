library(segmented)

fit_bp_segmented = function(depth, fluor, CD_thresh = 0){
  #CHL50 = CHL_50(depth[is.finite(fluor)],fluor[is.finite(fluor)])
  # EMLD = Eco_MLD(depth[is.finite(fluor)],prof_data$CTDFLUOR[is.finite(fluor)])
  # if(EMLD$QI > 0.3){
  #   depth = depth[which(depth < (EMLD$EMLD + 50))]
  #   fluor = fluor[which(depth < (EMLD$EMLD + 50))]
  # }
  test = new.env()
  assign("br",F,env = test)
  sp = 1
  fluor.in = fluor
  depth.in = depth
  na.idx = is.finite(fluor.in) & is.finite(depth.in)
  fluor = fluor[na.idx]
  depth = depth[na.idx]
  fit_lm = lm(fluor~depth)
  tryCatch({mod = segmented.lm(fit_lm, seg.Z = ~depth, npsi = sp)},error = function(e){
    print("cannot fit segmented profile with one break point")
    assign("br",T,env = test)})
  if(get("br",test)){return(NA)}else{
  BIC_new = BIC(mod)
  TSS = tot_var(fluor)
  # calculate first RSS
  RSS_new = RSS(fluor, mod$fitted.values)
  # coefficient of determination or r2 
  # Note - can have less than 0 as the first fit can be worse than the mean
  CD_new = (TSS-RSS_new)/TSS
  BIC_old = BIC_new
  CD_old = 0
  RSS_old = RSS_new
  assign("br",F,env = test)
  
  while( (RSS_old >= RSS_new) & ((CD_new - CD_old) > CD_thresh) & (BIC_old >= BIC_new) & !get("br",test)){
    
    sp = sp+1
    BIC_old = BIC_new
    CD_old = CD_new
    RSS_old = RSS_new
    tryCatch({mod = segmented.lm(fit_lm, seg.Z = ~depth, npsi = sp)
    BIC_new = BIC(mod)
    RSS_new = RSS(fluor, mod$fitted.values)
    CD_new = (TSS-RSS_new)/TSS
    }, warning=function(w){ 
      assign("br",T,env = test)}, error = function(e){assign("br",T,env = test)})
 
  }
  
  final_sp = sp - 1
  #Two chances for better fit
  assign("br",F,env = test)
  chance = 1
  rm.out = T
  final_out = F
  while(chance <= 5){
    sp = sp+1
    tryCatch({mod1 = segmented.lm(fit_lm, seg.Z = ~depth, npsi = sp)
    BIC_new = BIC(mod1)
    RSS_new = RSS(fluor, mod1$fitted.values)
    CD_new = (TSS-RSS_new)/TSS
    }, warning=function(w){
      assign("br",T,env = test)}, error = function(e){assign("br",T,env = test)})
    
    
    if((RSS_old >= RSS_new) & ((CD_new - CD_old) > CD_thresh) & (BIC_old >= BIC_new) & !get("br",test)){
      while( (RSS_old >= RSS_new) & ((CD_new - CD_old) > CD_thresh) & (BIC_old >= BIC_new) & !get("br",test)){
        sp= sp+1
        # test this
        # reset iteration
        chance = 1 
        BIC_old = BIC_new
        CD_old = CD_new
        RSS_old = RSS_new
        tryCatch({mod = segmented.lm(fit_lm, seg.Z = ~depth, npsi = sp)
        BIC_new = BIC(mod)
        RSS_new = RSS(fluor, mod$fitted.values)
        CD_new = (TSS-RSS_new)/TSS
        }, warning=function(w){
          assign("br",T,env = test)}, error = function(e){assign("br",T,env = test)})
       # print(paste("sp = ",sp))
        
      }
      final_sp = sp-1
      
    }
    assign("br",F,env = test)
    # print(paste("ch =",chance))
    chance = chance + 1
    
  if(chance == 5 & rm.out == T){
   # print(paste("final sp = ",final_sp))
  final_mod = segmented.lm(fit_lm, seg.Z = ~depth, npsi = final_sp)
  final_bps = final_mod$psi[,2]
  bps = final_mod$psi[,2]
  # check for outliers, remove and re-fit
  out_idx = NULL
  sp = final_sp
  for(bp in 1:(length(final_bps)+1)){
    if(bp == 1){idx = which(depth < final_bps[bp])
    idx1 = 1}else{
      if(bp == length(final_bps)+1){idx1 = which(depth >= final_bps[bp-1])[1]
      idx2 = length(depth)}else{
      idx1 = which(depth >= final_bps[bp-1])[1]
    idx2 = rev(which(depth < final_bps[bp]))[1]}
    idx = idx1:idx2}
    s = sd(final_mod$residuals[idx])
    m = mean(final_mod$residuals[idx])
    out_idx = c(out_idx, idx1 -1 + which(final_mod$residuals[idx] < (m - 3*s) | final_mod$residuals[idx] > (m + 3*s)))
  }
  
  if(length(out_idx) > 0){
    final_out = T
    back_up_mod = final_mod
    fluor_o = fluor[-out_idx]
    depth_o = depth[-out_idx]
    fit_lm = lm(fluor_o~depth_o)
    tryCatch({final_mod = segmented.lm(fit_lm, seg.Z = ~depth_o, fixed.psi = bps)
    final_bps = bps},error = function(e){assign("final_mod", back_up_mod,env = .GlobalEnv)
    assign("final_bps", back_up_mod$psi[,2],env = .GlobalEnv)
    assign("out_idx", NULL,env = .GlobalEnv)
    })
    # sometimes removing outliers screws it up.. so we need a bacckup where we dont remove outliers in this case
    if(!any(grepl(pattern = "psi",names(final_mod)))){final_mod = back_up_mod
    final_bps = back_up_mod$psi[,2]
    out_idx = NULL}else{
      if(!is.null(out_idx)){
      fluor = fluor_o
      depth = depth_o
      rm.out_sp = final_sp}}
    BIC_new = BIC(final_mod)
    RSS_new = RSS(fluor, final_mod$fitted.values)
    CD_new = (TSS-RSS_new)/TSS
    if(is.null(out_idx)){final_out = F}
  }
  rm.out = F
  #reset and resume iteration
  chance = 1
  }
    
    if(exists("rm.out_sp")){
      if(final_sp > (rm.out_sp + 5) ){
        final_mod = segmented.lm(fit_lm, seg.Z = ~depth, npsi = final_sp)
        final_bps = final_mod$psi[,2]
        
      }
    }
    
    
  }
  rm(test)
  fluor.out = fluor.in
  if(final_out){
    fluor.out[which(na.idx)[out_idx]] = NA}
  fluor.out[!na.idx] = NA
  if(length(out_idx) != 0 | final_out){
  fluor.out[-c(which(na.idx)[out_idx],which(!na.idx))] = final_mod$fitted.values
  }else{
    fluor.out[which(na.idx)] = final_mod$fitted.values}
  return(list("mod" = mod, "BIC" = BIC_new, "RSS" = RSS_new, "fluor.out" = fluor.out, "bp" = final_bps, "n_out" = length(out_idx), "idx_out" = out_idx))
  }
}

