front_class_insitu = function(press_prof, temp_prof, sal_prof,lat,lon){
  prof_array = data.frame(cbind(press_prof, swTheta(salinity = sal_prof, temperature = temp_prof, pressure = press_prof,latitude = lat, longitude = lon), temp_prof, sal_prof))
  prof_array = prof_array[complete.cases(prof_array), ]
  colnames(prof_array) = c("press","temp", "temp2", "sal")
  if(max(press_prof,na.rm = T) < 500){return("max press < 500 dbar")}else{
    
    t_min  = min(prof_array$temp, na.rm = T)
    t_max = max(prof_array$temp, na.rm = T)
    
    ub_plb = 390;
    ub_pub = 410;
    
    ## AZ
    # Tmin 
    t3_ub = 0.98
    
    # Tmax
    t3_ub_max = 2.11
    
    tmin_az = any(t_min < t3_ub)
    tmax_az = any(t_max < t3_ub_max)
    fun_az = tmin_az == 1 & tmax_az == 1
    isaz = fun_az
    
    ## PFZ
    pfz_ub = 2.78
    
    tmin_pfz = any(t_min > t3_ub)
    tmax_pfz = any(t_max > t3_ub_max)
    t_400  = any(prof_array$temp < pfz_ub & (prof_array$press >= ub_plb & prof_array$press <= ub_pub))
    fun_pfz = tmax_pfz == 1 & t_400 == 1
    
    ispfz = fun_pfz
    
    ## SAFZ
    saf_ub = 6.06
    saf_lb = 2.78
    
    fun_safz_ub  = any(prof_array$temp > saf_lb & prof_array$temp < saf_ub & (prof_array$press >= ub_plb & prof_array$press <= ub_pub))
    issafz = fun_safz_ub
    
    ## STZ/SAZ
    stz_lb = 6.06;
    stz_ub_plb = 70;
    stz_ub_pub = 130;
    
    fun_saz_ub  = any(prof_array$temp > stz_lb & (prof_array$press >= ub_plb & prof_array$press <= ub_pub))
    
    fun_stz_ub  = any(prof_array$temp2 > 11 & prof_array$sal >35 & (prof_array$press >= stz_ub_plb & prof_array$press <= stz_ub_pub))
    isstz = fun_stz_ub ==1 & fun_saz_ub == 1
    issaz = fun_stz_ub == 0 & fun_saz_ub == 1
    
    #combine result of classification tests
    result = c(isstz,issaz, issafz, ispfz, isaz)
    
    isunclass = sum(result) > 1
    isnoclass = sum(result) == 0
    if(isunclass){
      if(sum(result) == 2 & abs(which(result == 1)[1] - which(result == 1)[2] ) == 1){
        return(sample(c("STZ","SAZ","SAFZ","PFZ","AZ")[which(result == 1)],1))
      }else{return("unclassified")}
    }
    if(isnoclass){return("no class")}
    if(sum(result) == 1){return(c("STZ","SAZ","SAFZ","PFZ","AZ")[which(result == 1)])}
  }
  
}

