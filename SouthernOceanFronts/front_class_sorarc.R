# this function classifies profiles into front zones

front_class_sorarc = function(press_prof, temp_prof, sal_prof){
  prof_array = data.frame(cbind(press_prof, sal_prof,temp_prof))
  prof_array = prof_array[complete.cases(prof_array), ]
  colnames(prof_array) = c("press","sal","temp")
  if(max(press_prof,na.rm = T) < 510){return("max press < 500 dbar")}else{
  # Start code
  ## STZ
  # Orsi1995
  stz_p  = 105; 
  stz_t_lb = 11.5; 
  stz_s_lb = 35.05;
  
  # P-S conditions
  fun_stz_sp = prof_array$press <= stz_p & prof_array$sal >= stz_s_lb
  # P-T conditions
  fun_stz_tp =prof_array$press <= stz_p & prof_array$temp >= stz_t_lb
  # both conditions
  fun_stz = apply(cbind(fun_stz_sp,fun_stz_tp),1,function(x){x[1] == T & x[2] == T})
  #classify
  isstz = any(fun_stz)

  
  ## SAZ characterisation
  # S < 34.55 at 100 dbar
  # T > 6.85 at 400 dbar
  sazs_ub = 34.6;
  sazs_ub_pub = 105;
  
  sazt_lb = 6.85;
  sazt_lb_pub = 405;
  
  # P-S condition
  fun_saz_s = any(prof_array$sal < sazs_ub & prof_array$press <= sazs_ub_pub)
  
  # P-T condition
  fun_saz_t = any(prof_array$temp > sazt_lb & prof_array$press <= sazt_lb_pub)
  
  # both
  fun_saz = (fun_saz_s == 1 & fun_saz_t == 1)
  issaz = fun_saz
  
  ## PZ characterisation
  # T < 2.63 at 400 dbar
  # T < 2 at 200 dbar
  pzt_ub = 2.63;
  pzt_ub_pub = 405;
  pzt_ub_plb = 395;
  
  pzt_lb = 2;
  pzt_lb_pub = 205;
  pzt_lb_plb = 195;
  
  # P-T condition 1
  fun_pz_ub = any(prof_array$temp < pzt_ub & (prof_array$press >= pzt_ub_plb & prof_array$press <= pzt_ub_pub))
  
  # P-T condition 2
  fun_pz_lb = any(prof_array$temp > pzt_lb & (prof_array$press >= pzt_lb_plb & prof_array$press <= pzt_lb_pub))
  
  # both
  fun_pz = (fun_pz_ub == 1 & fun_pz_lb == 1)
  ispz = fun_pz
  
  ## AZ characterisation
  # T > 1.8 at 500 dbar
  # T < 2 at 200 dbar
  
  az_lb = 1.8; 
  az_lb_plb = 490; 
  az_lb_pub = 510; 
  
  az_ub = 2; 
  az_ub_plb = 190; 
  az_ub_pub = 210; 
  
  # P-T condition 1
  fun_az_lb = any(prof_array$temp > az_lb & (prof_array$press >= az_lb_plb & prof_array$press <= az_lb_pub))
  
  # P-T condition 2
  fun_az_ub = any(prof_array$temp < az_ub & (prof_array$press >= az_ub_plb & prof_array$press <= az_ub_pub))
  
  # both
  fun_az = (fun_az_ub == 1 & fun_az_lb == 1)
  isaz = fun_az 
  
  ## SZ characterisation
  # T > -0.66 at Tmin
  # T < ~1.8 at 500 dbar
  szt_lb = -0.66;
  szt_ub = 1.75;
  szt_ub_plb = 490;
  szt_ub_pub = 510;
  
  # P-T condition 1
  fun_szt_lb = any(min(prof_array$temp) >= szt_lb)
  
  # P-T condition 2
  fun_szt_ub  = any(prof_array$temp <= szt_ub & (prof_array$press >= szt_ub_plb & prof_array$press <= szt_ub_pub))
  
  # both
  fun_sz = (fun_szt_ub == 1 &  fun_szt_lb == 1)
  issz = fun_sz 
  
  ## Subpolar region
  # Tmin <= -1.24 C
  spr_t_lb = -1.24;
  
  fun_spr = any(prof_array$temp < spr_t_lb)
  isspr = fun_spr
  
  #combine result of classification tests
  result = c(isstz, issaz, ispz, isaz, issz, isspr)
  
  isunclass = sum(result) > 1
  isnoclass = sum(result) == 0
  if(isunclass){
    if(sum(result) == 2 & abs(which(result == 1)[1] - which(result == 1)[2] ) == 1){
    return(sample(c("stz","saz","pz","az","sz","spr")[which(result == 1)],1))
  }else{return("unclassified")}}
  if(isnoclass){return("no class")}
  if(sum(result) == 1){return(c("stz","saz","pz","az","sz","spr")[which(result == 1)])}
  #return(paste(c("stz","saz","pz","az","sz","spr")[which(result == 1)],sep = ","))
  }
}

